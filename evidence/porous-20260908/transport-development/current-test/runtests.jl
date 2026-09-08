using Test,LinearAlgebra
include("../src/PorousMedia.jl")
using .PorousMedia
const P=PorousMedia
@testset "Exact disk cut geometry" begin
    c=Circle(.53,.47,.21)
    for n in (20,40,80)
        g=cut_geometry(n,n,1.0,1.0,[c])
        @test abs(sum(1 .-g.volume)/n^2-pi*c.radius^2)<1e-12
        @test minimum(g.volume)>=0 && maximum(g.volume)<=1
        @test all((g.volume.==0).|(g.volume.==1))==false
        # Exact line intersections, not a mask of face centres.
        x=.5;expected=1-2sqrt(c.radius^2-(x-c.x)^2)
        @test sum(g.aperture_x[n÷2+1,:])/n≈expected atol=1e-13
        fx=[sin(i+.7j) for i in 1:n+1,j in 1:n];fy=[cos(.8i-j) for i in 1:n,j in 1:n+1]
        div=cut_divergence(fx,fy,g)
        net=(sum(g.aperture_x[end,:].*fx[end,:])-sum(g.aperture_x[1,:].*fx[1,:]))/n+
            (sum(g.aperture_y[:,end].*fy[:,end])-sum(g.aperture_y[:,1].*fy[:,1]))/n
        @test sum(g.volume.*div)/n^2≈net atol=1e-12
    end
    g=cut_geometry(40,40,1.,1.,[Circle(.5,0.,.2)])
    @test sum(1 .-g.volume)/40^2≈pi*.2^2/2 atol=1e-12
    @test length(paper_circles())==18
    @test paper_circles()[3].y==.00202
end

@testset "Cut-cell explicit pressure has a bounded relaxation spectrum" begin
    g=cut_geometry(100,60,.005,.003,paper_circles());cv=P.control_volumes(g)
    rho=[(i+j)%3==0 ? 998. : 1.2 for i in 1:g.nx,j in 1:g.ny]
    P.agglomerate!(rho,g,cv);a=P.pressure_graph(g,cv,rho)
    exact=[sin(.71k) for k in eachindex(cv.volumes)];b=similar(exact);P.graph_apply!(b,exact,a)
    density=P.volume_integrals(rho,g,cv)./cv.volumes;weight=cv.volumes./density
    p=zeros(length(exact));old=sum(weight.*(p.-exact).^2)
    report=P.explicit_pressure!(p,b,a,g,cv,rho,1e-5;iterations=1,sound_cfl=4.)
    @test report.sound_speed_scale<1
    @test report.relaxation_spectral_upper_bound<=1.8*(1+1e-14)
    @test sum(weight.*(p.-exact).^2)<old
    # Check contraction at each iteration, not just bounded scalar parameters.
    contraction=true
    for _ in 1:100
        old=sum(weight.*(p.-exact).^2)
        P.explicit_pressure!(p,b,a,g,cv,rho,1e-5;iterations=1)
        contraction &= sum(weight.*(p.-exact).^2)<=old*(1+1e-13)
    end
    @test contraction
    uf=zeros(g.nx+1,g.ny);vf=zeros(g.nx,g.ny+1);fill!(p,0)
    r=P.project_cut_weak!(uf,vf,p,1e-5,g,cv,rho)
    d=P.integrated_divergence(uf,vf,g,cv)
    @test norm(d)≈1e-5*r.residual rtol=1e-12
    @test uf[1,:]==fill(.003,g.ny)
end

@testset "Cut-volume transport conserves phase and consistent momentum" begin
    g=cut_geometry(100,60,.005,.003,paper_circles());cv=P.control_volumes(g)
    c=P.TransportParameters(inlet_u=.003,inlet_v=-.001,no_slip_outer_walls=false)
    phi=[.5*(1-tanh(((i-.5)*g.dx-.0005)/(2g.dx))) for i in 1:g.nx,j in 1:g.ny]
    P.agglomerate!(phi,g,cv)
    uf=zeros(g.nx+1,g.ny);vf=zeros(g.nx,g.ny+1);pressure=zeros(length(cv.volumes))
    P.project_cut!(uf,vf,pressure,1e-5,g,cv,fill(998.,size(phi));backend=:multigrid)
    # Manufactured uniform transported velocity is independent of the
    # divergence-free advecting field. No physical no-slip claim in this test.
    u=fill(c.inlet_u,size(phi));v=fill(c.inlet_v,size(phi));dt=1e-5
    rho=@. c.rho_g+(c.rho_l-c.rho_g)*phi
    m0=sum(P.volume_integrals(phi,g,cv));px0=sum(P.volume_integrals(rho.*u,g,cv));py0=sum(P.volume_integrals(rho.*v,g,cv))
    r=P.transport_step!(phi,u,v,uf,vf,dt,g,cv,c;wall_ghosts=false)
    function exterior(fx,fy)
        g.dy*(sum(g.aperture_x[end,:].*fx[end,:])-sum(g.aperture_x[1,:].*fx[1,:]))+
        g.dx*(sum(g.aperture_y[:,end].*fy[:,end])-sum(g.aperture_y[:,1].*fy[:,1]))
    end
    rho=@. c.rho_g+(c.rho_l-c.rho_g)*phi
    @test sum(P.volume_integrals(phi,g,cv))-m0≈dt*exterior(r.fluxes.fx,r.fluxes.fy) atol=1e-20
    @test sum(P.volume_integrals(rho.*u,g,cv))-px0≈dt*exterior(r.fluxes.ux,r.fluxes.uy) atol=1e-20
    @test sum(P.volume_integrals(rho.*v,g,cv))-py0≈dt*exterior(r.fluxes.vx,r.fluxes.vy) atol=1e-20
    @test maximum(abs.(u[cv.owner.>0].-c.inlet_u))<1e-11
    @test maximum(abs.(v[cv.owner.>0].-c.inlet_v))<1e-11
    # Boundary extension and momentum ghosts must never mutate conserved data.
    old=copy(phi);u.=0;v.=0;P.transport_fluxes(phi,u,v,uf,vf,g,c)
    @test phi==old && all(iszero,u) && all(iszero,v)
    c0=P.TransportParameters(inlet_phi=.37)
    phi.=.37;u.=0;v.=0
    P.transport_step!(phi,u,v,uf,vf,dt,g,cv,c0)
    @test maximum(abs.(phi[cv.owner.>0].-.37))<1e-11
end

@testset "Viscous force includes open faces and circle traction" begin
    g=cut_geometry(60,40,1.,1.,Circle[]);phi=ones(size(g.volume))
    u=[((j-.5)*g.dy)*(1-(j-.5)*g.dy) for i in 1:g.nx,j in 1:g.ny];v=zeros(size(u))
    c=P.TransportParameters(mu_l=.7,inlet_u=0.)
    f=P.viscous_force(phi,u,v,g,P.WallArc[],c)
    @test maximum(abs.(f.fx[3:end-2,3:end-2]./(g.dx*g.dy).+1.4))<1e-10
    @test maximum(abs,f.fy[3:end-2,3:end-2])<1e-12
    circles=[Circle(.51,.49,.21)];g=cut_geometry(80,80,1.,1.,circles);cv=P.control_volumes(g)
    phi=ones(size(g.volume));u=[sin(.1i)*cos(.2j) for i in 1:g.nx,j in 1:g.ny];v=.3 .*u
    f=P.viscous_force(phi,u,v,g,P.wall_arcs(g,circles),c)
    exterior(fx,fy)=g.dy*(sum(g.aperture_x[end,:].*fx[end,:])-sum(g.aperture_x[1,:].*fx[1,:]))+
        g.dx*(sum(g.aperture_y[:,end].*fy[:,end])-sum(g.aperture_y[:,1].*fy[:,1]))
    @test sum(f.fx)≈exterior(f.ux,f.uy)+sum(f.wallx) atol=1e-11
    @test sum(f.fy)≈exterior(f.vx,f.vy)+sum(f.wally) atol=1e-11
    # Exercise the conservative predictor's actual force accumulation path.
    uf=zeros(g.nx+1,g.ny);vf=zeros(g.nx,g.ny+1)
    mass0=sum(P.volume_integrals(998 .*u,g,cv));dt=1e-7
    P.transport_step!(phi,u,v,uf,vf,dt,g,cv,c;integrated_force=(f.fx,f.fy),gamma=0.)
    @test sum(P.volume_integrals(998 .*u,g,cv))-mass0≈dt*sum(f.fx) atol=1e-11
end

@testset "Smooth advection rate on an unobstructed grid converges" begin
    errors=Float64[]
    for n in (40,80,160)
        g=cut_geometry(n,8,1.,.2,Circle[]);cv=P.control_volumes(g)
        phi=[.5+.2sin(2pi*(i-.5)/n) for i in 1:n,j in 1:8]
        u=fill(.1,size(phi));v=zeros(size(phi));uf=fill(.1,n+1,8);vf=zeros(n,9)
        f=P.transport_fluxes(phi,u,v,uf,vf,g,P.TransportParameters();gamma=0.)
        rate=P.integrated_divergence(f.fx,f.fy,g,cv)./cv.volumes
        expected=[-.04pi*cos(2pi*(i-.5)/n) for i in 3:n-2,j in 1:8]
        numerical=[rate[cv.owner[i,j]] for i in 3:n-2,j in 1:8]
        push!(errors,norm(numerical-expected)/sqrt(length(expected)))
    end
    @test errors[2]<.27errors[1]
    @test errors[3]<.27errors[2]
end
@testset "Stationary immersed wall, manufactured normal-linear velocity" begin
    g=cut_geometry(80,80,1.,1.,[Circle(.5,.5,.23)])
    u=copy(g.distance);v=2 .*g.distance
    u[g.distance.<0].=0;v[g.distance.<0].=0
    its,res=impose_ghost_velocity!(u,v,g)
    @test res<1e-12
    ids=(-2g.dx .<g.distance).&(g.distance.<0)
    @test maximum(abs.(u[ids].-g.distance[ids]))<2g.dx^2/.23
    @test maximum(abs.(v[ids].-2g.distance[ids]))<4g.dx^2/.23
    @test u[g.distance.>=0]==g.distance[g.distance.>=0]
end
@testset "Liquid-side contact angle, planar affine interface extension" begin
    n=80;h=1/n;distance=[(j-.5)*h-.5 for i in 1:n,j in 1:n]
    wx=zeros(n,n);wy=ones(n,n)
    for theta in (30.,90.,150.)
        # phi need not be bounded for this affine manufactured extension test.
        exact=[sinpi(theta/180)*((i-.5)*h-.5)-cospi(theta/180)*distance[i,j] for i in 1:n,j in 1:n]
        q=copy(exact)
        for j in 1:n,i in 1:n
            distance[i,j]<0 && (q[i,j]=exact[i,n÷2+1])
        end
        fluid=copy(q[distance.>=0]);extend_contact!(q,distance,wx,wy,h,h,theta;iterations=160,band=4h)
        j=n÷2
        angles=[liquid_contact_angle((q[i+1,j]-q[i-1,j])/(2h),(q[i,j+1]-q[i,j-1])/(2h),0.,1.) for i in 15:65]
        @test maximum(abs.(angles.-theta))<.1
        @test q[distance.>=0]==fluid
    end
end
@testset "Conservative tiny-cell agglomeration and open-boundary projection" begin
    g=cut_geometry(60,36,.005,.003,paper_circles());cv=P.control_volumes(g)
    @test sum(cv.volumes)≈sum(g.volume)*g.dx*g.dy atol=1e-18
    q=[sin(.3i+.7j) for i in 1:g.nx,j in 1:g.ny];mass=sum(q.*g.volume)
    P.agglomerate!(q,g,cv);@test sum(q.*g.volume)≈mass atol=1e-11
    q.=1;P.agglomerate!(q,g,cv);@test maximum(abs.(q.-1))<1e-14
    rho=fill(998.,size(g.volume));a=P.pressure_graph(g,cv,rho)
    p=[sin(k*.1) for k in eachindex(cv.volumes)];b=similar(p);P.graph_apply!(b,p,a)
    sol=zeros(length(p));P.graph_pcg!(sol,b,a)
    @test norm(sol-p)/norm(p)<1e-7
    uf=zeros(g.nx+1,g.ny);vf=zeros(g.nx,g.ny+1);fill!(sol,0)
    P.project_cut!(uf,vf,sol,1e-5,g,cv,rho)
    d=P.integrated_divergence(uf,vf,g,cv)
    @test maximum(abs.(d./cv.volumes))<1e-5
    @test sum(g.aperture_x[end,:].*uf[end,:])≈sum(g.aperture_x[1,:].*uf[1,:]) rtol=1e-8
    @test all(iszero,uf[g.aperture_x.==0]) && all(iszero,vf[g.aperture_y.==0])
end
@testset "Paper grid has no phantom fluid in wholly solid cells" begin
    coarse=cut_geometry(250,150,.005,.003,paper_circles())
    @test coarse.volume[46,12]==0
    @test minimum(P.control_volumes(coarse).volumes)>=.5coarse.dx*coarse.dy*(1-1e-10)
    g=cut_geometry(1250,750,.005,.003,paper_circles())
    @test g.volume[123,5]==0
    @test sum(g.volume)*g.dx*g.dy≈.005*.003-15pi*.0004^2 atol=1e-15
    cv=P.control_volumes(g)
    @test minimum(cv.volumes)>=.5g.dx*g.dy*(1-1e-10)
end
@testset "Cut-cell graph multigrid is symmetric positive and solves the same equation" begin
    g=cut_geometry(60,36,.005,.003,paper_circles());cv=P.control_volumes(g)
    rho=[1.2+996.8*(1+sin(.1i+.15j))/2 for i in 1:g.nx,j in 1:g.ny]
    a=P.pressure_graph(g,cv,rho);mg=P.graph_mg(a,cv,g);n=length(cv.volumes)
    x=[sin(.17k) for k in 1:n];y=[cos(.23k) for k in 1:n];zx=zeros(n);zy=zeros(n)
    P.graph_precondition!(zx,x,mg);P.graph_precondition!(zy,y,mg)
    @test dot(x,zx)>0
    @test dot(x,zy)≈dot(y,zx) rtol=1e-10 atol=1e-8
    b=zeros(n);P.graph_apply!(b,x,a);sol=zeros(n)
    its,res=P.graph_pcg!(sol,b,a;mg=mg)
    @test norm(sol-x)/norm(x)<1e-7
    @test norm(P.graph_apply!(similar(b),sol,a)-b)<=1e-10norm(b)
end
@testset "Cut boundary normals close every fluid control volume" begin
    g=cut_geometry(80,80,1.,1.,[Circle(.51,.49,.213)])
    arcs=P.wall_arcs(g,[Circle(.51,.49,.213)])
    xx=(g.aperture_x[2:end,:]-g.aperture_x[1:end-1,:])*g.dy
    yy=(g.aperture_y[:,2:end]-g.aperture_y[:,1:end-1])*g.dx
    for a in arcs
        x,y=P.wall_normal_integral(a);xx[a.i,a.j]+=x;yy[a.i,a.j]+=y
    end
    @test maximum(abs,xx)<1e-12
    @test maximum(abs,yy)<1e-12
    @test sum(a.radius*(a.angle1-a.angle0) for a in arcs)≈2pi*.213 atol=1e-12
end
@testset "Fluid-side curved contact curvature excludes arbitrary solid values" begin
    errors=Float64[]
    for n in (200,400)
        rs=.23;rl=.12;sep=sqrt(rs^2+rl^2-2rs*rl*cosd(150.));sx=.4;sy=.5;cx=sx+sep
        g=cut_geometry(n,n,1.,1.,[Circle(sx,sy,rs)]);h=1/n
        psi=[rl-hypot((i-.5)*h-cx,(j-.5)*h-sy) for i in 1:n,j in 1:n]
        fits=P.boundary_fits(g);k=zeros(n,n);P.boundary_curvature!(k,psi,g,fits)
        xc=(rs^2-rl^2+sep^2)/(2sep);yc=sqrt(rs^2-xc^2)
        value=P.interpolate(k,sx+xc,sy+yc,h,h);push!(errors,abs(value-1/rl)/(1/rl))
        old=copy(k);psi[g.distance.<0].=12345.
        P.boundary_curvature!(k,psi,g,fits)
        @test k==old
    end
    @test errors[2]<.01
    @test errors[2]<.4errors[1]
end
