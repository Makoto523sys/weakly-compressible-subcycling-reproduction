"""Generate auditable plots/metrics for completed or partial runs; never certify them."""
import csv
import hashlib
import json
import pathlib
import tomllib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "evidence/local-20260907"
OUT.mkdir(parents=True, exist_ok=True)

def read_csv(path):
    with path.open() as f:
        reader = csv.DictReader(f)
        rows = [r for r in reader if all(r.get(k) is not None for k in reader.fieldnames)]
    return {k: np.array([float(r[k]) for r in rows]) for k in (reader.fieldnames or ["time"])}

report = {"status": "incomplete_requires_review", "bubble": {}, "static": {}}
fig, axes = plt.subplots(3, 1, figsize=(9, 10), sharex=True)
for ax, mode, curve in zip(axes, ["standard", "weak", "proposed"], ["standard", "weak", "proposed"]):
    run_name = "bubble-standard-mg" if mode == "standard" else f"bubble-{mode}-initial"
    path = ROOT / "results" / run_name
    reference = ROOT / f"reference/bubble_{curve}.csv"
    ref = read_csv(reference)
    ax.plot(ref["time"], ref["rise_velocity"], color="black", lw=1, label=f"Figure 5 {curve} (digitized)")
    if path.exists():
        sim = read_csv(path / "history.csv")
        meta = tomllib.loads((path / "run.toml").read_text())
        if len(sim["time"]) > 1:
            end = min(sim["time"][-1], ref["time"][-1], 0.07)
            t = np.linspace(max(sim["time"][0], ref["time"][0]), end, 1001)
            err = np.interp(t, sim["time"], sim["rise_velocity"]) - np.interp(t, ref["time"], ref["rise_velocity"])
            scale = max(abs(ref["rise_velocity"]))
            report["bubble"][mode] = {
                "run_status": meta["status"], "comparison_end_s": float(end),
                "full_interval": bool(end >= 0.07 - 1e-12),
                "nrmse": float(np.sqrt(np.sum((err[:-1]**2+err[1:]**2)*np.diff(t)/2)/(t[-1]-t[0]))/scale),
                "max_error_m_s": float(max(abs(err))),
                "area_drift": float(max(abs(sim["gas_area"]/sim["gas_area"][0]-1))),
                "phi_min": float(min(sim["phi_min"])), "phi_max": float(max(sim["phi_max"])),
                "max_divergence_per_s": float(max(sim["divergence_linf"])),
                "source_sha256": meta["source_sha256"],
                "history_sha256": hashlib.sha256((path/"history.csv").read_bytes()).hexdigest(),
                "reference_sha256": hashlib.sha256(reference.read_bytes()).hexdigest(),
            }
            ax.plot(sim["time"], sim["rise_velocity"], label=f"Julia {mode}: {end:.5f} s")
    ax.set_ylabel("Rise velocity [m/s]"); ax.legend(); ax.grid(alpha=.3)
axes[-1].set_xlabel("Time [s]")
fig.suptitle("Diagnostic comparison — reproduction NOT established")
fig.tight_layout(); fig.savefig(OUT/"bubble-overlay.png", dpi=160); plt.close(fig)

fig, axes = plt.subplots(1, 2, figsize=(10, 4))
for path in sorted((ROOT/"results").glob("static-*/history-*.csv")):
    if (path.parent/"PROVENANCE_WARNING.md").exists():
        continue
    d = read_csv(path)
    if not len(d["time"]):
        continue
    label = f"{path.parent.name}, nx={path.stem.split('-')[-1]}"
    report["static"][label] = {k: float(v[-1]) for k, v in d.items()}
    axes[0].plot(d["time"], d["pressure_jump_Pa"], label=label)
    axes[1].plot(d["time"], d["max_face_speed"], label=label)
axes[0].axhline(57.6, ls="--", color="black", label="Laplace 57.6 Pa")
axes[0].set_ylabel("Pressure jump [Pa]")
axes[1].set_ylabel("Face speed upper bound [m/s]")
for ax in axes:
    ax.set_xlabel("Time [s]"); ax.grid(alpha=.3)
axes[0].legend(fontsize=6)
fig.tight_layout(); fig.savefig(OUT/"static-balance.png", dpi=160); plt.close(fig)
(OUT/"metrics.json").write_text(json.dumps(report, indent=2)+"\n")
print(json.dumps(report, indent=2))
