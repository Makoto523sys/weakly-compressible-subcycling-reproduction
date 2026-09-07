# Diagnostic initial-force projection only; NOT a modified paper-case simulation.
using WeaklyCompressibleSubcycling, LinearAlgebra
const W=WeaklyCompressibleSubcycling
BLAS.set_num_threads(1)
println("nx,epsilon_over_dx,epsilon_m,projected_acceleration_bound,pressure_jump_Pa,pressure_residual")
for (n,m) in ((128,1.0),(256,1.0),(512,1.0),(256,2.0),(512,4.0))
    c=Config(nx=n,ny=2n,gravity=0.0,pressure_backend=:multigrid)
    s=bubble(c); epsilon=m*W.h(c)
    for j in 1:c.ny,i in 1:c.nx
        r=hypot((i-.5)*W.dx(c)-c.center_x,(j-.5)*W.dy(c)-c.center_y)
        s.phi[i,j]=.5*(1+tanh((r-c.radius)/(2epsilon)))
    end
    fx,fy=W.forces(s.phi,c);rx,ry=W.face_density(s.phi,c);dt=capillary_dt(c)
    W.add_face_increment!(s,dt.*fx./rx,dt.*fy./ry,c);W.project!(s,dt,c)
    inside=[hypot((i-.5)*W.dx(c)-c.center_x,(j-.5)*W.dy(c)-c.center_y)<c.radius/2 for i in 1:c.nx,j in 1:c.ny]
    outside=[hypot((i-.5)*W.dx(c)-c.center_x,(j-.5)*W.dy(c)-c.center_y)>2c.radius for i in 1:c.nx,j in 1:c.ny]
    jump=sum(s.p[inside])/count(inside)-sum(s.p[outside])/count(outside)
    println(join((n,m,epsilon,W.speed(s)/dt,jump,s.last_residual),','));flush(stdout)
end
