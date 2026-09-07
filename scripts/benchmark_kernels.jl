using LinearAlgebra
BLAS.set_num_threads(1)
root=dirname(@__DIR__)
old=Module(:Baseline);candidate=Module(:Candidate)
Base.include(old,joinpath(root,"evidence/local-20260907/pre-kernel-optimization/src/WeaklyCompressibleSubcycling.jl"))
Base.include(candidate,joinpath(root,"evidence/local-20260907/kernel-candidate/src/WeaklyCompressibleSubcycling.jl"))
O=old.WeaklyCompressibleSubcycling;C=candidate.WeaklyCompressibleSubcycling
oc=O.Config(nx=256,ny=512,pressure_backend=:multigrid)
cc=C.Config(nx=256,ny=512,pressure_backend=:multigrid)
os=O.bubble(oc);cs=C.bubble(cc)
for (name,mod,s,c) in (("baseline",O,os,oc),("candidate",C,cs,cc))
    mod.fluxes(s,c)
    m=@timed for _ in 1:20;mod.fluxes(s,c);end
    println((name=name,flux_seconds=m.time/20,bytes=m.bytes/20));flush(stdout)
end
for k in 1:3
    O.step!(os,O.capillary_dt(oc),oc;mode=:standard)
    C.step!(cs,C.capillary_dt(cc),cc;mode=:standard)
    for field in fieldnames(O.State)
        a=getfield(os,field);b=getfield(cs,field)
        a==b || error("Kernel change differs in $field at step $k")
    end
end
println("All state fields bitwise identical after 3 paper-grid steps")
