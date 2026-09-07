using WeaklyCompressibleSubcycling
using LinearAlgebra, TOML
const W = WeaklyCompressibleSubcycling
BLAS.set_num_threads(1)

function main(args)
    output = isempty(args) ? "results/static-local" : args[1]
    ispath(output) && error("Output exists: $output")
    mkpath(output)
    rows = Any[]
    # Same fluids, radius and box as the rising bubble; gravity alone is disabled
    # for this separate Laplace verification problem. Start from p=u=0.
    for nx in (64,128,256)
        c=Config(nx=nx,ny=2nx,gravity=0.0)
        s=bubble(c); area0=diagnostics(s,c).gas_area
        inside=[hypot((i-0.5)*W.dx(c)-c.center_x,
                      (j-0.5)*W.dy(c)-c.center_y)<c.radius/2
                for i in 1:c.nx,j in 1:c.ny]
        outside=[hypot((i-0.5)*W.dx(c)-c.center_x,
                       (j-0.5)*W.dy(c)-c.center_y)>2c.radius
                 for i in 1:c.nx,j in 1:c.ny]
        start=time_ns(); peak=0.0
        open(joinpath(output,"history-$nx.csv"),"w") do io
            println(io,"time,pressure_jump_Pa,relative_jump_error,max_face_speed,gas_area_drift,divergence_linf,phi_min,phi_max")
            while s.t < 0.001
                step!(s,min(stable_dt(s,c),0.001-s.t),c;mode=:standard)
                d=diagnostics(s,c)
                jump=sum(s.p[inside])/count(inside)-sum(s.p[outside])/count(outside)
                vmax=W.speed(s); peak=max(peak,vmax)
                println(io,join((s.t,jump,abs(jump-c.sigma/c.radius)/(c.sigma/c.radius),vmax,
                    (d.gas_area-area0)/area0,d.divergence_linf,d.phi_min,d.phi_max),','))
            end
        end
        elapsed=(time_ns()-start)/1e9
        row=Dict("nx"=>nx,"steps"=>s.steps,"elapsed_seconds"=>elapsed,
                 "peak_face_speed_bound"=>peak,"source_sha256"=>source_digest())
        push!(rows,row); println(row); flush(stdout)
        W.write_vtk(joinpath(output,"final-$nx.vtk"),s,c)
        open(joinpath(output,"summary.toml"),"w") do io
            TOML.print(io,Dict("status"=>"measured_requires_review","runs"=>rows,
                "julia_version"=>string(VERSION),"t_end"=>0.001,
                "expected_pressure_jump_Pa"=>c.sigma/c.radius))
        end
    end
end
main(ARGS)
