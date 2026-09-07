using WeaklyCompressibleSubcycling,LinearAlgebra
const W=WeaklyCompressibleSubcycling
BLAS.set_num_threads(1)
println("nx,backend,iterations,residual,elapsed_seconds,allocated_bytes")
for n in (32,128,256), backend in (:jacobi,:multigrid)
    c=Config(nx=n,ny=2n,pressure_backend=backend)
    s=bubble(c);rx,ry=W.face_density(s.phi,c)
    inc=W.increments(s,capillary_dt(c),c)
    W.provisional!(s,deepcopy(s),inc...,c)
    s.uf[1,:].=0;s.uf[end,:].=0;s.vf[:,1].=0;s.vf[:,end].=0
    b=-W.divergence(s.uf,s.vf,c)./capillary_dt(c)
    W.poisson!(s.p,copy(b),rx,ry,c)
    fill!(s.p,0)
    measurement=@timed W.poisson!(s.p,copy(b),rx,ry,c)
    println(join((n,backend,measurement.value...,measurement.time,measurement.bytes),','));flush(stdout)
end
