module WeaklyCompressibleSubcycling

using LinearAlgebra
using Printf
using TOML
using SHA
using Serialization

export Config, State, bubble, capillary_dt, stable_dt, step!, diagnostics,
       run_bubble, write_vtk, source_digest

function save_checkpoint(path,s,c,mode,factor,digest)
    open(path*".tmp","w") do io
        serialize(io,(state=s,config=c,mode=mode,factor=factor,
            source_sha256=digest,julia_version=string(VERSION)))
    end
    mv(path*".tmp",path;force=true)
end

function load_checkpoint(path)
    saved=open(deserialize,path)
    saved.julia_version==string(VERSION) || error("Checkpoint Julia version differs")
    saved.source_sha256==source_digest() || error("Checkpoint solver revision differs")
    saved
end

"""SI units; phi=0 gas, phi=1 liquid. Unvalidated research implementation."""
Base.@kwdef struct Config
    nx::Int = 256
    ny::Int = 512
    lx::Float64 = 0.010
    ly::Float64 = 0.020
    rho_g::Float64 = 1.2
    rho_l::Float64 = 998.0
    mu_g::Float64 = 1.8e-5
    mu_l::Float64 = 1.0e-3
    sigma::Float64 = 0.072
    gravity::Float64 = -9.8
    radius::Float64 = 0.00125
    center_x::Float64 = 0.005
    center_y::Float64 = 0.005
    epp_iterations::Int = 10
    sound_cfl::Float64 = 0.4
    poisson_rtol::Float64 = 1e-9
    poisson_atol::Float64 = 1e-10
    poisson_maxiter::Int = 10000
    weno_epsilon::Float64 = 1e-12
    pressure_backend::Symbol = :jacobi
end

dx(c::Config) = c.lx/c.nx
dy(c::Config) = c.ly/c.ny
h(c::Config) = min(dx(c), dy(c))
density(phi, c) = c.rho_g .+ (c.rho_l-c.rho_g).*phi
viscosity(phi, c) = c.mu_g .+ (c.mu_l-c.mu_g).*phi

mutable struct State
    phi::Matrix{Float64}
    u::Matrix{Float64}
    v::Matrix{Float64}
    p::Matrix{Float64}
    uf::Matrix{Float64}
    vf::Matrix{Float64}
    t::Float64
    steps::Int
    substeps::Int
    poisson_solves::Int
    last_residual::Float64
end

function bubble(c::Config)
    c.nx >= 4 && c.ny >= 4 || error("At least four cells per direction required")
    isapprox(dx(c), dy(c); rtol=1e-12) || error("This implementation requires square cells")
    c.radius > 0 || error("Bubble radius must be positive")
    min(c.rho_g,c.rho_l) > 0 || error("Density must be positive")
    min(c.mu_g,c.mu_l) >= 0 || error("Viscosity must be nonnegative")
    c.sigma >= 0 || error("Surface tension must be nonnegative")
    0 < c.sound_cfl <= 0.4 || error("EPP sound CFL must lie in (0,0.4]")
    c.epp_iterations >= 1 || error("EPP iteration count must be positive")
    phi = [0.5*(1+tanh((hypot((i-0.5)*dx(c)-c.center_x,
              (j-0.5)*dy(c)-c.center_y)-c.radius)/(2h(c))))
           for i in 1:c.nx, j in 1:c.ny]
    State(phi, zeros(c.nx,c.ny), zeros(c.nx,c.ny), zeros(c.nx,c.ny),
          zeros(c.nx+1,c.ny), zeros(c.nx,c.ny+1), 0.0, 0, 0, 0, 0.0)
end

# Reflection supplies the WENO halo. Scalars: zero normal derivative.
# u: odd on all walls. v: even on x walls (free slip), odd on y walls.
@inline function sample(q, i, j, kind::Symbol=:scalar)
    nx,ny = size(q)
    sign = 1.0
    if i < 1
        i = 1-i
        kind in (:u,:nx) && (sign = -sign)
    elseif i > nx
        i = 2nx+1-i
        kind in (:u,:nx) && (sign = -sign)
    end
    if j < 1
        j = 1-j
        kind in (:u,:v,:ny) && (sign = -sign)
    elseif j > ny
        j = 2ny+1-j
        kind in (:u,:v,:ny) && (sign = -sign)
    end
    sign*q[i,j]
end

function interface_geometry(phi,c)
    epsilon = h(c)
    psi = @. epsilon*log((clamp(phi,0,1)+1e-100)/(1-clamp(phi,0,1)+1e-100))
    nx = similar(phi); ny = similar(phi); kappa = similar(phi)
    for j in 1:c.ny, i in 1:c.nx
        a = (sample(psi,i+1,j)-sample(psi,i-1,j))/(2dx(c))
        b = (sample(psi,i,j+1)-sample(psi,i,j-1))/(2dy(c))
        d = hypot(a,b)+1e-100
        nx[i,j] = a/d; ny[i,j] = b/d
    end
    for j in 1:c.ny, i in 1:c.nx
        kappa[i,j] = -(sample(nx,i+1,j,:nx)-sample(nx,i-1,j,:nx))/(2dx(c)) -
                     (sample(ny,i,j+1,:ny)-sample(ny,i,j-1,:ny))/(2dy(c))
    end
    psi,nx,ny,kappa
end

# A conservative bound on the face speed magnitude; exact norm reconstruction
# in the authors' implementation is not specified. Documented implementation choice.
function speed(s)
    nx,ny=size(s.phi);largest=0.0
    # Reconstruct the tangential component at each normal-velocity face.
    @inbounds for j in 1:ny,i in 1:nx+1
        il=max(i-1,1);ir=min(i,nx)
        tangential=(s.vf[il,j]+s.vf[ir,j]+s.vf[il,j+1]+s.vf[ir,j+1])/4
        largest=max(largest,hypot(s.uf[i,j],tangential))
    end
    @inbounds for j in 1:ny+1,i in 1:nx
        tangential=if j==1 || j==ny+1
            0.0 # no-slip horizontal walls
        else
            (s.uf[i,j-1]+s.uf[i+1,j-1]+s.uf[i,j]+s.uf[i+1,j])/4
        end
        largest=max(largest,hypot(tangential,s.vf[i,j]))
    end
    largest
end
capillary_dt(c::Config) = c.sigma == 0 ? Inf : sqrt((c.rho_g+c.rho_l)*h(c)^3/(4pi*c.sigma))
function stable_dt(s::State,c::Config)
    umax = speed(s)
    adv = umax == 0 ? Inf : h(c)/umax
    nu = max(c.mu_g/c.rho_g,c.mu_l/c.rho_l)
    visc = nu == 0 ? Inf : h(c)^2/(4nu)
    phase = umax == 0 ? Inf : h(c)/(4umax)
    min(0.95adv,0.5visc,capillary_dt(c),0.5phase)
end

@inline function weno_left(a,b,d,epsw)
    beta0 = (b-a)^2; beta1 = (d-b)^2
    w0 = (1/3)/(epsw+beta0)^2
    w1 = (2/3)/(epsw+beta1)^2
    (w0*(1.5b-0.5a)+w1*(0.5b+0.5d))/(w0+w1)
end
@inline function face_value(q,i,j,axis,velocity,kind,epsw)
    if axis == 1
        a=sample(q,i-2,j,kind); b=sample(q,i-1,j,kind)
        d=sample(q,i,j,kind); e=sample(q,i+1,j,kind)
    else
        a=sample(q,i,j-2,kind); b=sample(q,i,j-1,kind)
        d=sample(q,i,j,kind); e=sample(q,i,j+1,kind)
    end
    velocity >= 0 ? weno_left(a,b,d,epsw) : weno_left(e,d,b,epsw)
end

function fluxes(s::State,c::Config)
    phi,u,v = s.phi,s.u,s.v
    psi,nx,ny,_ = interface_geometry(phi,c)
    mu = viscosity(phi,c)
    gamma = speed(s)
    fx=zeros(c.nx+1,c.ny); fy=zeros(c.nx,c.ny+1)
    ux=similar(fx); vx=similar(fx); uy=similar(fy); vy=similar(fy)
    for j in 1:c.ny, i in 1:c.nx+1
        vel=s.uf[i,j]
        ph=(sample(phi,i-1,j)+sample(phi,i,j))/2
        ps=(sample(psi,i-1,j)+sample(psi,i,j))/2
        normal=(sample(nx,i-1,j,:nx)+sample(nx,i,j,:nx))/2
        fx[i,j] = -vel*ph + gamma*(h(c)*(sample(phi,i,j)-sample(phi,i-1,j))/dx(c) -
                   0.25*(1-tanh(ps/(2h(c)))^2)*normal)
        # Closed boundary: zero phase/mass flux, including regularization.
        (i==1 || i==c.nx+1) && (fx[i,j]=0)
        mass=-c.rho_g*vel+(c.rho_l-c.rho_g)*fx[i,j]
        muf=(sample(mu,i-1,j)+sample(mu,i,j))/2
        dudx=(sample(u,i,j,:u)-sample(u,i-1,j,:u))/dx(c)
        dvdx=(sample(v,i,j,:v)-sample(v,i-1,j,:v))/dx(c)
        dudy=(sample(u,i,j+1,:u)+sample(u,i-1,j+1,:u)-
              sample(u,i,j-1,:u)-sample(u,i-1,j-1,:u))/(4dy(c))
        ux[i,j]=mass*face_value(u,i,j,1,vel,:u,c.weno_epsilon)+2muf*dudx
        vx[i,j]=mass*face_value(v,i,j,1,vel,:v,c.weno_epsilon)+muf*(dudy+dvdx)
    end
    for j in 1:c.ny+1, i in 1:c.nx
        vel=s.vf[i,j]
        ph=(sample(phi,i,j-1)+sample(phi,i,j))/2
        ps=(sample(psi,i,j-1)+sample(psi,i,j))/2
        normal=(sample(ny,i,j-1,:ny)+sample(ny,i,j,:ny))/2
        fy[i,j] = -vel*ph + gamma*(h(c)*(sample(phi,i,j)-sample(phi,i,j-1))/dy(c) -
                   0.25*(1-tanh(ps/(2h(c)))^2)*normal)
        (j==1 || j==c.ny+1) && (fy[i,j]=0)
        mass=-c.rho_g*vel+(c.rho_l-c.rho_g)*fy[i,j]
        muf=(sample(mu,i,j-1)+sample(mu,i,j))/2
        dudy=(sample(u,i,j,:u)-sample(u,i,j-1,:u))/dy(c)
        dvdy=(sample(v,i,j,:v)-sample(v,i,j-1,:v))/dy(c)
        dvdx=(sample(v,i+1,j,:v)+sample(v,i+1,j-1,:v)-
              sample(v,i-1,j,:v)-sample(v,i-1,j-1,:v))/(4dx(c))
        uy[i,j]=mass*face_value(u,i,j,2,vel,:u,c.weno_epsilon)+muf*(dudy+dvdx)
        vy[i,j]=mass*face_value(v,i,j,2,vel,:v,c.weno_epsilon)+2muf*dvdy
    end
    (;fx,fy,ux,uy,vx,vy)
end

function divergence(fx,fy,c)
    [(fx[i+1,j]-fx[i,j])/dx(c)+(fy[i,j+1]-fy[i,j])/dy(c)
     for i in 1:c.nx,j in 1:c.ny]
end
function face_density(phi,c)
    rho=density(phi,c)
    rx=[(sample(rho,i-1,j)+sample(rho,i,j))/2 for i in 1:c.nx+1,j in 1:c.ny]
    ry=[(sample(rho,i,j-1)+sample(rho,i,j))/2 for i in 1:c.nx,j in 1:c.ny+1]
    rx,ry
end

function forces(phi,c)
    _,_,_,k=interface_geometry(phi,c)
    rx,ry=face_density(phi,c)
    fx=zeros(c.nx+1,c.ny); fy=zeros(c.nx,c.ny+1)
    for j in 1:c.ny, i in 1:c.nx+1
        ph=(sample(phi,i-1,j)+sample(phi,i,j))/2
        kap=(sample(k,i-1,j)+sample(k,i,j))/2
        fx[i,j]=6ph*(1-ph)*c.sigma*kap*(sample(phi,i,j)-sample(phi,i-1,j))/dx(c)
    end
    for j in 1:c.ny+1, i in 1:c.nx
        ph=(sample(phi,i,j-1)+sample(phi,i,j))/2
        kap=(sample(k,i,j-1)+sample(k,i,j))/2
        fy[i,j]=6ph*(1-ph)*c.sigma*kap*(sample(phi,i,j)-sample(phi,i,j-1))/dy(c) +
                ry[i,j]*c.gravity
    end
    fx,fy
end

function add_face_increment!(s,ax,ay,c)
    s.uf .+= ax; s.vf .+= ay
    for j in 1:c.ny,i in 1:c.nx
        s.u[i,j] += (ax[i,j]+ax[i+1,j])/2
        s.v[i,j] += (ay[i,j]+ay[i,j+1])/2
    end
end

function provisional!(s,old,iphi,imu,imv,ifx,ify,c)
    s.phi .= old.phi .+ iphi
    rho=density(s.phi,c)
    minimum(rho)>0 || error("Nonpositive density; refusing to clip the transported phase")
    minimum(viscosity(s.phi,c))>=0 || error("Negative viscosity")
    rho0=density(old.phi,c)
    s.u .= (rho0.*old.u .+ imu)./rho
    s.v .= (rho0.*old.v .+ imv)./rho
    for j in 1:c.ny,i in 1:c.nx+1
        s.uf[i,j]=(sample(s.u,i-1,j,:u)+sample(s.u,i,j,:u))/2
    end
    for j in 1:c.ny+1,i in 1:c.nx
        s.vf[i,j]=(sample(s.v,i,j-1,:v)+sample(s.v,i,j,:v))/2
    end
    rx,ry=face_density(s.phi,c)
    add_face_increment!(s,ifx./rx,ify./ry,c)
end

# Positive semidefinite Neumann operator A = -div((1/rho_f) grad).
# The pressure constant nullspace is projected out, never fixed by a penalty.
function pressure_operator!(out,p,rx,ry,c)
    fill!(out,0)
    for j in 1:c.ny, i in 2:c.nx
        f=(p[i,j]-p[i-1,j])/(rx[i,j]*dx(c)^2)
        out[i,j]+=f; out[i-1,j]-=f
    end
    for j in 2:c.ny, i in 1:c.nx
        f=(p[i,j]-p[i,j-1])/(ry[i,j]*dy(c)^2)
        out[i,j]+=f; out[i,j-1]-=f
    end
    out
end
zero_mean!(q) = (q .-= sum(q)/length(q); q)

include("Multigrid.jl")

function poisson!(p,b,rx,ry,c)
    zero_mean!(b); zero_mean!(p)
    diag=zeros(size(p))
    for j in 1:c.ny,i in 2:c.nx
        a=1/(rx[i,j]*dx(c)^2); diag[i,j]+=a; diag[i-1,j]+=a
    end
    for j in 2:c.ny,i in 1:c.nx
        a=1/(ry[i,j]*dy(c)^2); diag[i,j]+=a; diag[i,j-1]+=a
    end
    ap=similar(p); pressure_operator!(ap,p,rx,ry,c)
    r=b.-ap; zero_mean!(r)
    tol=max(c.poisson_atol,c.poisson_rtol*norm(b))
    norm(r)<=tol && return (0,norm(r))
    c.pressure_backend in (:jacobi,:multigrid) || error("Unknown pressure backend")
    mg=c.pressure_backend==:multigrid ? mg_setup(rx,ry,c) : nothing
    z=zero_mean!(r./diag)
    mg!==nothing && mg_precondition!(z,r,mg)
    d=copy(z); rz=dot(r,z)
    for iteration in 1:c.poisson_maxiter
        pressure_operator!(ap,d,rx,ry,c)
        denominator=dot(d,ap)
        denominator>0 || error("PCG lost positive definiteness")
        alpha=rz/denominator
        p .+= alpha.*d; r .-= alpha.*ap
        # Verify against the true residual before reporting convergence.
        if norm(r)<=tol || iteration%50==0
            pressure_operator!(ap,p,rx,ry,c)
            r .= b.-ap; zero_mean!(r)
            if norm(r)<=tol
                zero_mean!(p)
                return iteration,norm(r)
            end
        end
        if mg===nothing
            z .= r./diag; zero_mean!(z)
        else
            mg_precondition!(z,r,mg)
        end
        rznew=dot(r,z)
        d .= z .+ (rznew/rz).*d
        rz=rznew
    end
    error("Pressure Poisson solver did not converge; residual=$(norm(r)), tolerance=$tol")
end

function project!(s,dt,c; weak=false)
    rx,ry=face_density(s.phi,c)
    # Lift prescribed wall fluxes before using the homogeneous Neumann operator.
    # Retain original wall velocity for the matching cell-center correction.
    wallx=copy(s.uf); wally=copy(s.vf)
    s.uf[1,:].=0; s.uf[end,:].=0
    s.vf[:,1].=0; s.vf[:,end].=0
    div=divergence(s.uf,s.vf,c)
    if weak
        rho=density(s.phi,c)
        cs=c.sound_cfl*h(c)/dt
        ap=similar(s.p)
        for _ in 1:c.epp_iterations
            pressure_operator!(ap,s.p,rx,ry,c)
            s.p .-= (dt*cs^2).*rho.*(div .+ dt.*ap)
            zero_mean!(s.p)
        end
    else
        _,residual=poisson!(s.p,-div./dt,rx,ry,c)
        s.last_residual=residual
        s.poisson_solves+=1
    end
    ax=zeros(size(s.uf)); ay=zeros(size(s.vf))
    for j in 1:c.ny,i in 2:c.nx
        ax[i,j]=-dt*(s.p[i,j]-s.p[i-1,j])/(dx(c)*rx[i,j])
    end
    for j in 2:c.ny,i in 1:c.nx
        ay[i,j]=-dt*(s.p[i,j]-s.p[i,j-1])/(dy(c)*ry[i,j])
    end
    # Restore uncorrected wall values so a single face/cell operation enforces BCs.
    s.uf[1,:].=wallx[1,:]; s.uf[end,:].=wallx[end,:]
    s.vf[:,1].=wally[:,1]; s.vf[:,end].=wally[:,end]
    ax[1,:].=-wallx[1,:]; ax[end,:].=-wallx[end,:]
    ay[:,1].=-wally[:,1]; ay[:,end].=-wally[:,end]
    add_face_increment!(s,ax,ay,c)
    nothing
end

function increments(s,dt,c)
    f=fluxes(s,c)
    iphi=dt.*divergence(f.fx,f.fy,c)
    imu=dt.*divergence(f.ux,f.uy,c)
    imv=dt.*divergence(f.vx,f.vy,c)
    fx,fy=forces(s.phi.+iphi,c)
    iphi,imu,imv,dt.*fx,dt.*fy
end

"""Advance one main step; accumulated quantities are flux/force integrals, not averaged velocities."""
function step!(s::State,dt,c::Config; mode::Symbol=:subcycling)
    mode in (:standard,:weak,:subcycling) || error("Unknown solver mode $mode")
    isfinite(dt) && dt>0 || error("Time step must be finite and positive")
    old=deepcopy(s)
    if mode==:subcycling
        work=deepcopy(s)
        sums=(zeros(size(s.phi)),zeros(size(s.u)),zeros(size(s.v)),
              zeros(size(s.uf)),zeros(size(s.vf)))
        elapsed=0.0; nsub=0
        while elapsed < dt
            remaining=dt-elapsed
            limit=stable_dt(work,c)
            # Absorb a rounding-sized excess into the last physical substep.
            # Otherwise an integer multiplier can produce a ~1e-19 s extra
            # EPP step, whose artificial sound speed is proportional to 1/dt.
            subdt=remaining<=limit*(1+64eps(Float64)) ? remaining : limit
            elapsed+subdt>elapsed || error("Substep underflow")
            inc=increments(work,subdt,c)
            for k in 1:5
                sums[k].+=inc[k]
            end
            prev=deepcopy(work)
            provisional!(work,prev,inc...,c)
            project!(work,subdt,c;weak=true)
            elapsed+=subdt; nsub+=1
            all(isfinite,work.phi) && all(isfinite,work.u) && all(isfinite,work.v) &&
                all(isfinite,work.p) || error("Nonfinite substep state")
        end
        provisional!(s,old,sums...,c)
        # Pressure is solved afresh from the accumulated RHS, initialized with p^n.
        project!(s,dt,c)
        s.substeps+=nsub
    else
        dt<=stable_dt(s,c)*(1+1e-12) || error("Baseline step exceeds explicit stability limit")
        inc=increments(s,dt,c)
        provisional!(s,old,inc...,c)
        project!(s,dt,c;weak=(mode==:weak))
        s.substeps+=1
    end
    s.t=old.t+dt; s.steps+=1
    all(isfinite,s.phi) && all(isfinite,s.u) && all(isfinite,s.v) && all(isfinite,s.p) ||
        error("Nonfinite state")
    s
end

function diagnostics(s,c)
    gas=1 .- s.phi
    area=sum(gas)*dx(c)*dy(c)
    area>0 || error("Nonpositive bubble area")
    y=sum(gas[i,j]*(j-0.5)*dy(c) for i in 1:c.nx,j in 1:c.ny)/sum(gas)
    rise=sum(gas.*s.v)/sum(gas)
    div=divergence(s.uf,s.vf,c)
    (;time=s.t, rise_velocity=rise, centroid_y=y, gas_area=area,
      phi_min=minimum(s.phi),phi_max=maximum(s.phi),
      divergence_linf=maximum(abs,div),
      kinetic_energy=sum(density(s.phi,c).*(s.u.^2+s.v.^2))*dx(c)*dy(c)/2,
      steps=s.steps,substeps=s.substeps,poisson_solves=s.poisson_solves,
      poisson_residual=s.last_residual)
end

function write_vtk(path,s,c)
    open(path,"w") do io
        println(io,"# vtk DataFile Version 3.0\nUnvalidated ACDI calculation\nASCII\nDATASET STRUCTURED_POINTS")
        println(io,"DIMENSIONS $(c.nx+1) $(c.ny+1) 1\nORIGIN 0 0 0\nSPACING $(dx(c)) $(dy(c)) 1")
        println(io,"CELL_DATA $(c.nx*c.ny)")
        for (name,q) in (("phi_liquid",s.phi),("pressure_Pa",s.p))
            println(io,"SCALARS $name double 1\nLOOKUP_TABLE default")
            for value in q
                println(io,value)
            end
        end
        println(io,"VECTORS velocity_m_per_s double")
        for j in 1:c.ny,i in 1:c.nx
            println(io,"$(s.u[i,j]) $(s.v[i,j]) 0")
        end
    end
end

function source_digest(source_dir=@__DIR__)
    root=dirname(source_dir)
    files=sort(filter(p->endswith(p,".jl"),[joinpath(d,f) for (d,_,fs) in walkdir(source_dir) for f in fs]))
    io=IOBuffer()
    for path in files
        write(io,relpath(path,root),UInt8(0),read(path),UInt8(0))
    end
    bytes2hex(sha256(take!(io)))
end

function run_bubble(c::Config; mode=:subcycling,factor=30.0,t_end=0.001,
                    output="results/bubble",snapshot_interval=0.001,resume=nothing)
    t_end>0 && isfinite(t_end) || error("Positive finite end time required")
    factor>0 && isfinite(factor) || error("Positive finite capillary multiplier required")
    snapshot_interval>0 || error("Positive output interval required")
    ispath(output) && error("Output already exists; choose a new directory: $output")
    mkpath(output)
    cp(@__DIR__,joinpath(output,"src"))
    cp(joinpath(dirname(@__DIR__),"Project.toml"),joinpath(output,"Project.toml"))
    s=if resume===nothing
        bubble(c)
    else
        saved=load_checkpoint(resume)
        saved.config==c && saved.mode==mode && saved.factor==factor ||
            error("Checkpoint configuration/mode/factor differs")
        deepcopy(saved.state)
    end
    s.t<t_end || error("End time must follow checkpoint time")
    meta=Dict{String,Any}("status"=>"running_unvalidated","julia_version"=>string(VERSION),
        "source_sha256"=>source_digest(),"mode"=>string(mode),"factor"=>factor,
        "requested_t_end"=>t_end,"capillary_dt"=>capillary_dt(c),
        "snapshot_interval"=>snapshot_interval,
        "initial_time"=>s.t,"resumed_from"=>(resume===nothing ? "" : abspath(resume)),
        "config"=>Dict(string(f)=>(getfield(c,f) isa Symbol ? string(getfield(c,f)) : getfield(c,f)) for f in fieldnames(Config)))
    function metadata!()
        open(joinpath(output,"run.toml"),"w") do io
            TOML.print(io,meta)
        end
    end
    metadata!()
    start=time_ns(); next_snapshot=s.t+snapshot_interval; snapshot=0
    try
        open(joinpath(output,"history.csv"),"w") do io
            snapshots=joinpath(output,"snapshots.csv")
            write(snapshots,"file,time\nfield_000000.vtk,$(s.t)\n")
            d=diagnostics(s,c); println(io,join(string.(keys(d)),',')); println(io,join(values(d),','))
            write_vtk(joinpath(output,"field_000000.vtk"),s,c)
            while s.t<t_end
                dt=mode==:subcycling ? factor*capillary_dt(c) : min(factor*capillary_dt(c),stable_dt(s,c))
                dt=min(dt,t_end-s.t)
                step!(s,dt,c;mode=mode)
                println(io,join(values(diagnostics(s,c)),','))
                if s.t>=next_snapshot || s.t>=t_end
                    snapshot+=1
                    write_vtk(joinpath(output,@sprintf("field_%06d.vtk",snapshot)),s,c)
                    open(snapshots,"a") do si
                        println(si,@sprintf("field_%06d.vtk,%.17g",snapshot,s.t))
                    end
                    save_checkpoint(joinpath(output,"checkpoint.jls"),s,c,mode,factor,meta["source_sha256"])
                    next_snapshot+=snapshot_interval
                    flush(io)
                    @printf("t=%.8g step=%d substeps=%d\n",s.t,s.steps,s.substeps)
                end
            end
        end
        meta["status"]="completed_unvalidated"
    catch err
        meta["status"]="failed"
        meta["error"]=sprint(showerror,err)
        rethrow()
    finally
        meta["elapsed_seconds"]=(time_ns()-start)/1e9
        meta["actual_t_end"]=s.t
        metadata!()
    end
    s
end

end
