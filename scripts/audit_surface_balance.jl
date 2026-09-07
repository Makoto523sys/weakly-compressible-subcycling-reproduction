using WeaklyCompressibleSubcycling, LinearAlgebra
const W=WeaklyCompressibleSubcycling
BLAS.set_num_threads(1)
# Diagnostic alternatives only: do not replace the paper's CSF in production.
println("nx,force_model,max_projected_face_acceleration,pressure_jump")
for n in (64,128,256), model in (:paper,:constant_curvature,:discrete_potential)
    c=Config(nx=n,ny=2n,gravity=0.0);s=bubble(c)
    fx,fy=W.forces(s.phi,c);rx,ry=W.face_density(s.phi,c)
    if model != :paper
        potential=@. -c.sigma/c.radius*(3s.phi^2-2s.phi^3)
        for j in 1:c.ny,i in 2:c.nx
            a=s.phi[i-1,j];b=s.phi[i,j];ph=(a+b)/2
            fx[i,j]=model==:constant_curvature ? -6ph*(1-ph)*c.sigma/c.radius*(b-a)/W.dx(c) :
                (potential[i,j]-potential[i-1,j])/W.dx(c)
        end
        for j in 2:c.ny,i in 1:c.nx
            a=s.phi[i,j-1];b=s.phi[i,j];ph=(a+b)/2
            fy[i,j]=model==:constant_curvature ? -6ph*(1-ph)*c.sigma/c.radius*(b-a)/W.dy(c) :
                (potential[i,j]-potential[i,j-1])/W.dy(c)
        end
    end
    dt=capillary_dt(c)
    W.add_face_increment!(s,dt.*fx./rx,dt.*fy./ry,c)
    W.project!(s,dt,c)
    println(join((n,model,W.speed(s)/dt,maximum(s.p)-minimum(s.p)),','));flush(stdout)
end
