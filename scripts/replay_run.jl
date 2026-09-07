using WeaklyCompressibleSubcycling, TOML

length(ARGS)==2 || error("Usage: julia --project=ARCHIVED_RUN scripts/replay_run.jl ORIGINAL_RUN_TOML NEW_OUTPUT")
meta=TOML.parsefile(ARGS[1])
source_digest()==meta["source_sha256"] || error("Activate the archived run's project; solver digest differs")
parameters=Pair{Symbol,Any}[]
for f in fieldnames(Config)
    value=meta["config"][string(f)]
    fieldtype(Config,f)==Symbol && (value=Symbol(value))
    push!(parameters,f=>value)
end
c=Config(;parameters...)
run_bubble(c;mode=Symbol(meta["mode"]),factor=meta["factor"],t_end=meta["requested_t_end"],
    output=ARGS[2],snapshot_interval=get(meta,"snapshot_interval",meta["requested_t_end"]/20))
