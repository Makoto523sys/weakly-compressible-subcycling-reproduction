using LinearAlgebra,TOML,SHA
include("src/PorousMedia.jl")
using .PorousMedia
const P=PorousMedia
BLAS.set_num_threads(1)
function main()
    output=isempty(ARGS) ? "results/porous-geometry-projection" : ARGS[1]
    ispath(output) && error("Output exists")
    mkpath(output);cp("porous/src",joinpath(output,"src"))
    hashes=Dict(f=>bytes2hex(sha256(read(joinpath("porous/src",f)))) for f in readdir("porous/src"))
    measurement=@timed cut_geometry(1250,750,.005,.003,paper_circles())
    g=measurement.value
    println("Geometry seconds: ",measurement.time," allocated bytes: ",measurement.bytes);flush(stdout)
    cvmeasure=@timed P.control_volumes(g);cv=cvmeasure.value
    println("Control volumes: ",length(cv.volumes)," seconds: ",cvmeasure.time);flush(stdout)
    vol=sum(g.volume)*g.dx*g.dy;expected=.005*.003-15pi*.0004^2
    abs(vol-expected)<1e-15 || error("Paper geometry area mismatch")
    uf=zeros(g.nx+1,g.ny);vf=zeros(g.nx,g.ny+1);p=zeros(length(cv.volumes));rho=fill(998.,size(g.volume))
    # Single-phase projection of the imposed inflow only; NOT a Navier-Stokes run.
    projection=@timed P.project_cut!(uf,vf,p,1e-5,g,cv,rho;backend=:multigrid)
    d=P.integrated_divergence(uf,vf,g,cv)
    inlet=sum(g.aperture_x[1,:].*uf[1,:])*g.dy;outlet=sum(g.aperture_x[end,:].*uf[end,:])*g.dy
    summary=Dict("pressure_backend"=>"Galerkin aggregation multigrid PCG","status"=>"geometry_and_single_projection_only","julia_version"=>string(VERSION),
        "source_sha256"=>hashes,"nx"=>g.nx,"ny"=>g.ny,"dx_m"=>g.dx,
        "fluid_area_m2"=>vol,"exact_fluid_area_m2"=>expected,"porosity"=>vol/(.005*.003),
        "minimum_raw_cut_fraction"=>minimum(filter(>(0),g.volume)),
        "control_volumes"=>length(cv.volumes),"minimum_merged_fraction"=>minimum(cv.volumes)/(g.dx*g.dy),
        "geometry_seconds"=>measurement.time,"geometry_allocated_bytes"=>measurement.bytes,
        "projection_seconds"=>projection.time,"projection_allocated_bytes"=>projection.bytes,
        "pressure_iterations"=>projection.value[1],"absolute_pressure_residual"=>projection.value[2],
        "group_divergence_linf_per_s"=>maximum(abs.(d./cv.volumes)),
        "inlet_area_flux_m2_s"=>inlet,"outlet_area_flux_m2_s"=>outlet,
        "relative_flow_mismatch"=>abs(outlet-inlet)/inlet,
        "note"=>"Constant liquid density, one pressure projection; no viscous flow, two-phase transport or contact-line simulation")
    open(joinpath(output,"summary.toml"),"w") do io;TOML.print(io,summary);end
    println(summary)
end
main()
