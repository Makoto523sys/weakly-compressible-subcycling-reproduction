using TOML,SHA,LinearAlgebra
length(ARGS)==2 || error("Usage: baseline_source_directory output_directory")
const baseline=abspath(ARGS[1]);const out=abspath(ARGS[2])
ispath(out) && error("Output exists")
mkpath(out)
module Before
    include(joinpath(Main.baseline,"PorousMedia.jl"))
end
include("src/PorousMedia.jl")
const B=Before.PorousMedia;const P=PorousMedia
BLAS.set_num_threads(1)
function main()
    g=P.cut_geometry(1250,750,.005,.003,P.paper_circles())
    q=[.5*(1-tanh(((i-.5)*g.dx-.0005)/(2g.dx))) for i in 1:g.nx,j in 1:g.ny]
    # Compile both paths before timing; use identical independent initial data.
    args=(g.distance,g.normal_x,g.normal_y,g.dx,g.dy,150.)
    before=copy(q);after=copy(q)
    B.extend_contact!(before,args...);P.extend_contact!(after,args...)
    before==after || error("Optimization changed the extension result")
    bt=Float64[];at=Float64[];bb=Int[];ab=Int[]
    for _ in 1:5
        before.=q;after.=q
        b=@timed B.extend_contact!(before,args...)
        a=@timed P.extend_contact!(after,args...)
        before==after || error("Result mismatch")
        push!(bt,b.time);push!(at,a.time);push!(bb,b.bytes);push!(ab,a.bytes)
    end
    summary=Dict("nx"=>g.nx,"ny"=>g.ny,"bitwise_equal"=>true,"before_seconds"=>bt,"after_seconds"=>at,
        "before_allocated_bytes"=>bb,"after_allocated_bytes"=>ab,"julia_version"=>string(VERSION),
        "baseline_source"=>baseline,"baseline_sha256"=>bytes2hex(sha256(read(joinpath(baseline,"PorousMedia.jl")))),
        "current_sha256"=>bytes2hex(sha256(read(joinpath(@__DIR__,"src/PorousMedia.jl")))),
        "scope"=>"80-iteration contact extension alone, one paper-grid phase input; not full CFD speedup")
    open(joinpath(out,"summary.toml"),"w") do io;TOML.print(io,summary);end
    cp(joinpath(@__DIR__,"src"),joinpath(out,"src"));cp(@__FILE__,joinpath(out,"benchmark_contact_extension.jl"))
    println(summary)
end
main()
