# Conservative ACDI and consistent advective momentum transport. Fluxes have
# the same negative-advection sign as the reviewed bubble solver. Solid ghosts
# are scratch values; they never replace the conserved cut-volume quantities.
Base.@kwdef struct TransportParameters
    rho_l::Float64=998.0
    rho_g::Float64=1.2
    theta::Float64=150.0
    inlet_phi::Float64=1.0
    inlet_u::Float64=.003
    inlet_v::Float64=0.0
    weno_epsilon::Float64=1e-12
end

@inline function transport_sample(q,i,j,kind,c)
    nx,ny=size(q)
    if i<1
        value=kind==:phi ? c.inlet_phi : kind==:u ? c.inlet_u : c.inlet_v
        return 2value-q[min(1-i,nx),clamp(j,1,ny)]
    end
    q[clamp(i,1,nx),clamp(j,1,ny)]
end
@inline function transport_weno(a,b,d,e,vel,epsw)
    aa,bb,dd=vel>=0 ? (a,b,d) : (e,d,b)
    w0=(1/3)/(epsw+(bb-aa)^2)^2;w1=(2/3)/(epsw+(dd-bb)^2)^2
    (w0*(1.5bb-.5aa)+w1*(.5bb+.5dd))/(w0+w1)
end
@inline function transported_velocity(q,i,j,axis,vel,kind,c)
    di,dj=axis==1 ? (1,0) : (0,1)
    a=transport_sample(q,i-2di,j-2dj,kind,c)
    b=transport_sample(q,i-di,j-dj,kind,c)
    d=transport_sample(q,i,j,kind,c)
    e=transport_sample(q,i+di,j+dj,kind,c)
    transport_weno(a,b,d,e,vel,c.weno_epsilon)
end

function transport_fluxes(phi,u,v,uf,vf,g,c::TransportParameters;
                          gamma=hypot(maximum(abs,uf),maximum(abs,vf)),wall_ghosts=true)
    gamma>=0 || error("Negative ACDI regularization speed")
    h=min(g.dx,g.dy);q=copy(phi);uu=copy(u);vv=copy(v)
    if wall_ghosts
        extend_contact!(q,g.distance,g.normal_x,g.normal_y,g.dx,g.dy,c.theta)
        impose_ghost_velocity!(uu,vv,g)
    end
    psi=@. h*log((clamp(q,0,1)+1e-100)/(1-clamp(q,0,1)+1e-100))
    nnx=zeros(size(q));nny=similar(nnx)
    # Phase halos follow the imposed inlet and zero-gradient outlet. The
    # regularization flux is explicitly zero on all exterior faces.
    for j in 1:g.ny,i in 1:g.nx
        gx=(psi[min(i+1,g.nx),j]-psi[max(i-1,1),j])/(2g.dx)
        gy=(psi[i,min(j+1,g.ny)]-psi[i,max(j-1,1)])/(2g.dy)
        r=hypot(gx,gy)
        # psi has length units, hence this threshold is dimensionless. Do not
        # amplify roundoff in a constant field into an O(1) sharpening flux.
        nnx[i,j]=r>1e-12 ? gx/r : 0.;nny[i,j]=r>1e-12 ? gy/r : 0.
    end
    fx=zeros(g.nx+1,g.ny);fy=zeros(g.nx,g.ny+1)
    mx=similar(fx);my=similar(fy);ux=similar(fx);uy=similar(fy);vx=similar(fx);vy=similar(fy)
    for (axis,flux,mass,momx,momy,velocities,aperture,normal,spacing) in
        ((1,fx,mx,ux,vx,uf,g.aperture_x,nnx,g.dx),(2,fy,my,uy,vy,vf,g.aperture_y,nny,g.dy))
        di,dj=axis==1 ? (1,0) : (0,1)
        for j in axes(flux,2),i in axes(flux,1)
            vel=velocities[i,j]
            if aperture[i,j]==0 || (axis==2 && (j==1 || j==g.ny+1))
                flux[i,j]=mass[i,j]=momx[i,j]=momy[i,j]=0;continue
            end
            if axis==1 && (i==1 || i==g.nx+1)
                # Outlet zero-gradient phase; backflow needs a specified
                # reservoir and is rejected instead of silently inventing one.
                i==g.nx+1 && vel< -1e-12 && error("Outlet backflow requires a phase boundary condition")
                ph=i==1 ? c.inlet_phi : phi[end,j]
                f=-vel*ph
                uface=i==1 ? c.inlet_u : u[end,j]
                vface=i==1 ? c.inlet_v : v[end,j]
            else
                a,b=i-di,j-dj
                ph=(q[a,b]+q[i,j])/2;ps=(psi[a,b]+psi[i,j])/2
                nr=(normal[a,b]+normal[i,j])/2
                f=-vel*ph+gamma*(h*(q[i,j]-q[a,b])/spacing-.25*(1-tanh(ps/(2h))^2)*nr)
                uface=transported_velocity(uu,i,j,axis,vel,:u,c)
                vface=transported_velocity(vv,i,j,axis,vel,:v,c)
            end
            m=-c.rho_g*vel+(c.rho_l-c.rho_g)*f
            flux[i,j]=f;mass[i,j]=m;momx[i,j]=m*uface;momy[i,j]=m*vface
        end
    end
    (;fx,fy,mx,my,ux,uy,vx,vy)
end

function volume_integrals(q,g,cv)
    out=zeros(length(cv.volumes))
    for k in eachindex(q)
        a=cv.owner[k];a==0 && continue
        out[a]+=g.volume[k]*g.dx*g.dy*q[k]
    end
    out
end
function scatter_control_volumes!(q,values,cv)
    for k in eachindex(q)
        a=cv.owner[k];a==0 && continue;q[k]=values[a]
    end
    q
end

"""Explicit conservative phase and advective momentum predictor only.
No viscosity, CSF, pressure correction or velocity boundary enforcement is
included here. Returns the exact numerical flux budget used by this step.
The caller supplies face velocities; this is not a two-phase CFD integrator.
"""
function transport_step!(phi,u,v,uf,vf,dt,g,cv,c::TransportParameters;kwargs...)
    dt>0 || error("Positive transport timestep required")
    f=transport_fluxes(phi,u,v,uf,vf,g,c;kwargs...)
    rp=integrated_divergence(f.fx,f.fy,g,cv)
    ru=integrated_divergence(f.ux,f.uy,g,cv);rv=integrated_divergence(f.vx,f.vy,g,cv)
    rho=@. c.rho_g+(c.rho_l-c.rho_g)*phi
    phin=(volume_integrals(phi,g,cv).+dt.*rp)./cv.volumes
    rhon=@. c.rho_g+(c.rho_l-c.rho_g)*phin
    all(isfinite,phin) && minimum(rhon)>0 || error("Transport generated invalid phase/density; no clipping applied")
    un=(volume_integrals(rho.*u,g,cv).+dt.*ru)./(rhon.*cv.volumes)
    vn=(volume_integrals(rho.*v,g,cv).+dt.*rv)./(rhon.*cv.volumes)
    scatter_control_volumes!(phi,phin,cv);scatter_control_volumes!(u,un,cv);scatter_control_volumes!(v,vn,cv)
    (;phase_rate=sum(rp),momentum_x_rate=sum(ru),momentum_y_rate=sum(rv),fluxes=f)
end
