@testset "Symmetric multigrid preconditioner and unchanged fine solve" begin
    c=Config(nx=32,ny=64,pressure_backend=:multigrid,poisson_rtol=1e-10)
    s=bubble(c);rx,ry=W.face_density(s.phi,c)
    mg=W.mg_setup(rx,ry,c)
    x=W.zero_mean!([sin(i+2j) for i in 1:c.nx,j in 1:c.ny])
    y=W.zero_mean!([cos(2i-j) for i in 1:c.nx,j in 1:c.ny])
    mx=similar(x);my=similar(y)
    W.mg_precondition!(mx,x,mg);W.mg_precondition!(my,y,mg)
    @test dot(x,mx)>0
    @test isapprox(dot(x,my),dot(mx,y);rtol=1e-10,atol=1e-18)
    b=similar(x);W.pressure_operator!(b,x,rx,ry,c)
    p=zeros(size(x));W.poisson!(p,copy(b),rx,ry,c)
    @test norm(p-x)/norm(x)<1e-7
    a=bubble(c);b=deepcopy(a);base=Config(nx=32,ny=64,poisson_rtol=1e-10)
    for _ in 1:3
        dt=min(stable_dt(a,c),stable_dt(b,base))
        step!(a,dt,c;mode=:standard);step!(b,dt,base;mode=:standard)
    end
    @test maximum(abs,a.phi-b.phi)<1e-9
    @test maximum(abs,a.v-b.v)<1e-8
    @test maximum(abs,a.p-b.p)<1e-5
    odd=Config(nx=19,ny=38,pressure_backend=:multigrid)
    ox,oy=W.face_density(bubble(odd).phi,odd)
    @test_throws ErrorException W.mg_setup(ox,oy,odd)
end
