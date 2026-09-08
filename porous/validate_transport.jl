using LinearAlgebra,TOML,SHA
include("src/PorousMedia.jl")
using .PorousMedia
const P=PorousMedia
BLAS.set_num_threads(1)
function main()
    out=length(ARGS)>0 ? ARGS[1] : "results/porous-prescribed-transport"
    nx=length(ARGS)>1 ? parse(Int,ARGS[2]) : 250
    tend=length(ARGS)>2 ? parse(Float64,ARGS[3]) : .02
    factor=length(ARGS)>3 ? parse(Float64,ARGS[4]) : 1.
    nx%5==0 && factor>0 || error("Invalid grid or timestep factor")
    ispath(out) && error("Output exists")
    mkpath(out);cp("porous/src",joinpath(out,"src"));cp(@__FILE__,joinpath(out,"validate_transport.jl"))
    hashes=Dict(f=>bytes2hex(sha256(read(joinpath("porous/src",f)))) for f in readdir("porous/src"))
    g=cut_geometry(nx,3nx÷5,.005,.003,paper_circles());cv=P.control_volumes(g)
    phi=[.5*(1-tanh(((i-.5)*g.dx-.0005)/(2g.dx))) for i in 1:g.nx,j in 1:g.ny]
    P.agglomerate!(phi,g,cv);u=zeros(size(phi));v=zeros(size(phi))
    uf=zeros(g.nx+1,g.ny);vf=zeros(g.nx,g.ny+1);p=zeros(length(cv.volumes))
    P.project_cut!(uf,vf,p,1e-5,g,cv,fill(998.,size(phi));backend=:multigrid)
    c=P.TransportParameters();gamma=hypot(maximum(abs,uf),maximum(abs,vf))
    # Conservative provisional restriction, not a proved nonlinear cut-cell
    # stability bound. Test halving and raw phase extrema, never clip phi.
    dt=factor*.0625min(g.dx,g.dy)/gamma
    initial=sum(P.volume_integrals(phi,g,cv));budget=0.;maxerror=0.;lo=minimum(phi[cv.owner.>0]);hi=maximum(phi[cv.owner.>0])
    t=0.;steps=0;elapsed=0.;allocated=0
    open(joinpath(out,"history.csv"),"w") do io
        println(io,"step,time_s,liquid_area_m2,boundary_integral_m2,balance_error_m2,phi_min,phi_max")
        println(io,"0,0,$initial,0,0,$lo,$hi")
        while t<tend
            delta=min(dt,tend-t)
            measurement=@timed P.transport_step!(phi,u,v,uf,vf,delta,g,cv,c)
            elapsed+=measurement.time;allocated+=measurement.bytes
            # u,v here are passive transported momentum and do not drive uf,vf.
            budget+=delta*measurement.value.phase_rate;t+=delta;steps+=1
            area=sum(P.volume_integrals(phi,g,cv));err=area-initial-budget
            low=minimum(phi[cv.owner.>0]);high=maximum(phi[cv.owner.>0])
            maxerror=max(maxerror,abs(err));lo=min(lo,low);hi=max(hi,high)
            println(io,"$steps,$t,$area,$budget,$err,$low,$high");flush(io)
        end
    end
    summary=Dict("status"=>"prescribed_velocity_transport_only_not_CFD","source_sha256"=>hashes,
        "nx"=>g.nx,"ny"=>g.ny,"time_s"=>t,"steps"=>steps,"nominal_dt_s"=>dt,
        "timestep_factor"=>factor,"speed_bound_m_s"=>gamma,"contact_angle_deg"=>c.theta,
        "liquid_area_initial_m2"=>initial,"liquid_area_final_m2"=>sum(P.volume_integrals(phi,g,cv)),
        "boundary_phase_integral_m2"=>budget,"max_balance_error_m2"=>maxerror,
        "max_relative_balance_error"=>maxerror/initial,"phi_min_all_steps"=>lo,"phi_max_all_steps"=>hi,
        "transport_seconds_including_first_JIT"=>elapsed,"transport_allocated_bytes"=>allocated,
        "group_divergence_linf_per_s"=>maximum(abs.(P.integrated_divergence(uf,vf,g,cv)./cv.volumes)),
        "julia_version"=>string(VERSION),"note"=>"Fixed single-density projected face velocity; no viscosity, CSF, two-phase pressure coupling or no-slip flow solution. Geometry is Figure 10, grid is diagnostic. No pressure or invasion reproduction claim.")
    open(joinpath(out,"summary.toml"),"w") do io;TOML.print(io,summary);end
    open(joinpath(out,"phase.csv"),"w") do io
        println(io,"i,j,fluid_fraction,phi")
        for j in 1:g.ny,i in 1:g.nx
            println(io,"$i,$j,$(g.volume[i,j]),$(phi[i,j])")
        end
    end
    println(summary)
end
main()
