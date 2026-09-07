using Test
using LinearAlgebra
using WeaklyCompressibleSubcycling
const W=WeaklyCompressibleSubcycling
include("multigrid.jl")

@testset "Checkpoint preserves all state and rejects changed configuration" begin
    c=Config(nx=16,ny=32,pressure_backend=:multigrid)
    s=bubble(c);dt=capillary_dt(c)/2
    step!(s,dt,c;mode=:standard)
    mktempdir() do dir
        p=joinpath(dir,"state.jls")
        W.save_checkpoint(p,s,c,:standard,0.5,source_digest())
        saved=W.load_checkpoint(p)
        for f in fieldnames(State)
            @test getfield(s,f)==getfield(saved.state,f)
        end
        step!(s,dt,c;mode=:standard)
        resumed=run_bubble(c;mode=:standard,factor=0.5,t_end=2dt,
            output=joinpath(dir,"resumed"),resume=p)
        @test s.phi==resumed.phi
        @test s.uf==resumed.uf
        @test s.p==resumed.p
        @test_throws ErrorException run_bubble(c;mode=:weak,factor=0.5,t_end=2dt,
            output=joinpath(dir,"wrong"),resume=p)
    end
end

@testset "No roundoff-sized extra EPP substep" begin
    c=Config(nx=32,ny=64,gravity=0.0,pressure_backend=:multigrid)
    for factor in (15,30)
        s=bubble(c)
        step!(s,factor*capillary_dt(c),c;mode=:subcycling)
        @test s.substeps==factor
        @test s.t==factor*capillary_dt(c)
        @test all(isfinite,s.p)
    end
end

@testset "Paper capillary timestep, SI units" begin
    c=Config()
    @test isapprox(capillary_dt(c),8.1e-6;rtol=0.01)
    @test isapprox(capillary_dt(Config(nx=128,ny=256))/capillary_dt(c),2^1.5;rtol=1e-14)
end

@testset "WENO reproduces constant and affine data" begin
    @test W.weno_left(3.0,3.0,3.0,1e-12)≈3.0
    @test W.weno_left(1.0,2.0,3.0,1e-12)≈2.5
end

@testset "Conservative divergence and phase flux" begin
    c=Config(nx=16,ny=32)
    fx=[sin(i+2j) for i in 1:c.nx+1,j in 1:c.ny]
    fy=[cos(2i-j) for i in 1:c.nx,j in 1:c.ny+1]
    fx[1,:].=0;fx[end,:].=0;fy[:,1].=0;fy[:,end].=0
    @test abs(sum(W.divergence(fx,fy,c))*W.dx(c)*W.dy(c))<1e-12
    s=bubble(c)
    s.uf.=fx;s.vf.=fy
    f=W.fluxes(s,c)
    @test abs(sum(W.divergence(f.fx,f.fy,c))*W.dx(c)*W.dy(c))<1e-12
end

@testset "Variable-density Neumann pressure solve" begin
    c=Config(nx=16,ny=32,poisson_rtol=1e-10)
    phi=[0.5+0.45sin(2pi*(i-0.5)/c.nx)*sin(2pi*(j-0.5)/c.ny)
          for i in 1:c.nx,j in 1:c.ny]
    rx,ry=W.face_density(phi,c)
    exact=[cos(pi*(i-0.5)/c.nx)*cos(2pi*(j-0.5)/c.ny)
           for i in 1:c.nx,j in 1:c.ny]
    b=similar(exact);W.pressure_operator!(b,exact,rx,ry,c)
    p=zeros(size(exact));W.poisson!(p,copy(b),rx,ry,c)
    @test norm(p-exact)/norm(exact)<1e-7
    s=bubble(c);s.phi.=phi
    for j in 1:c.ny,i in 2:c.nx
        s.uf[i,j]=sin(0.3i+0.2j)
    end
    for j in 2:c.ny,i in 1:c.nx
        s.vf[i,j]=cos(0.4i-0.1j)
    end
    before=norm(W.divergence(s.uf,s.vf,c))
    W.project!(s,1e-5,c)
    @test norm(W.divergence(s.uf,s.vf,c))/before<1e-8
end

@testset "Pressure operator analytical Neumann eigenmode" begin
    c=Config(nx=16,ny=32)
    phi=ones(c.nx,c.ny)
    rx,ry=W.face_density(phi,c)
    p=[cos(pi*(i-0.5)/c.nx)*cos(2pi*(j-0.5)/c.ny)
       for i in 1:c.nx,j in 1:c.ny]
    eigenvalue=(4sin(pi/(2c.nx))^2/W.dx(c)^2 +
                4sin(pi/c.ny)^2/W.dy(c)^2)/c.rho_l
    ap=similar(p);W.pressure_operator!(ap,p,rx,ry,c)
    @test norm(ap-eigenvalue*p)/norm(ap)<1e-12
end

@testset "Single-phase hydrostatic balance, including wall lifting" begin
    c=Config(nx=16,ny=32)
    s=bubble(c);fill!(s.phi,1.0)
    dt=0.1capillary_dt(c)
    step!(s,dt,c;mode=:standard)
    @test maximum(abs,s.uf)<1e-9
    @test maximum(abs,s.vf)<1e-9
    @test maximum(abs,s.u)<1e-9
    @test maximum(abs,s.v)<1e-9
    expected=c.rho_l*c.gravity*W.dy(c)
    @test maximum(abs,diff(s.p;dims=2).-expected)<1e-5
end

@testset "One-substep identity and bubble mass conservation" begin
    c=Config(nx=32,ny=64)
    a=bubble(c); b=deepcopy(a)
    dt=0.25stable_dt(a,c)
    mass=sum(a.phi)
    step!(a,dt,c;mode=:standard)
    step!(b,dt,c;mode=:subcycling)
    @test b.substeps==1
    @test a.phi≈b.phi
    @test a.u≈b.u atol=1e-12
    @test a.v≈b.v atol=1e-12
    @test abs(sum(a.phi)-mass)<1e-10
    @test all(isfinite,a.p)
end

@testset "Outward bubble normal gives negative signed curvature" begin
    c=Config(nx=128,ny=256)
    s=bubble(c)
    _,_,_,k=W.interface_geometry(s.phi,c)
    mask=(s.phi.>0.4).&(s.phi.<0.6)
    @test abs(sum(k[mask])/count(mask)+1/c.radius)/(1/c.radius)<0.05
end
