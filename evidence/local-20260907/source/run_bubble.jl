using WeaklyCompressibleSubcycling

function main(args)
    if isempty(args) || args[1] in ("-h","--help")
        println("Usage: julia --project=. scripts/run_bubble.jl MODE NX FACTOR T_END OUTPUT")
        println("MODE = standard | weak | subcycling. NY = 2*NX. All units SI.")
        println("A short run is only a smoke test, not paper validation.")
        return
    end
    length(args)==5 || error("Expected five arguments; use --help")
    mode=Symbol(args[1]); nx=parse(Int,args[2]); factor=parse(Float64,args[3])
    tend=parse(Float64,args[4]); output=args[5]
    mode in (:standard,:weak,:subcycling) || error("Unknown mode")
    c=Config(nx=nx,ny=2nx)
    run_bubble(c;mode=mode,factor=factor,t_end=tend,output=output,
               snapshot_interval=max(tend/20,1e-12))
end

main(ARGS)
