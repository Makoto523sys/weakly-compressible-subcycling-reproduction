"""Compare two prescribed-velocity transport runs, never label them CFD."""
import hashlib
import json
from pathlib import Path
import shutil
import sys
import tomllib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

left, right, output = map(Path, sys.argv[1:])
output.mkdir(parents=True, exist_ok=False)
summaries = [tomllib.loads((p / "summary.toml").read_text()) for p in (left, right)]
assert all(s["status"] == "prescribed_velocity_transport_only_not_CFD" for s in summaries)
assert all(summaries[0][k] == summaries[1][k] for k in ("nx", "ny", "time_s"))
data = [np.genfromtxt(p / "phase.csv", delimiter=",", names=True) for p in (left, right)]
assert np.array_equal(data[0]["fluid_fraction"], data[1]["fluid_fraction"])
w = data[0]["fluid_fraction"]
error = data[1]["phi"] - data[0]["phi"]
report = {
    "scope": "Timestep comparison of prescribed-velocity phase transport, not two-phase CFD",
    "volume_weighted_l1_phase_difference": float(np.sum(w * abs(error)) / sum(w)),
    "volume_weighted_rms_phase_difference": float(np.sqrt(np.sum(w * error**2) / sum(w))),
    "fluid_linf_phase_difference": float(np.max(abs(error[w > 0]))),
    "inputs_sha256": {
        str(p / "phase.csv"): hashlib.sha256((p / "phase.csv").read_bytes()).hexdigest()
        for p in (left, right)
    },
}
(output / "comparison.json").write_text(json.dumps(report, indent=2) + "\n")
ny, nx = summaries[0]["ny"], summaries[0]["nx"]
fig, axes = plt.subplots(1, 3, figsize=(12, 3.3), constrained_layout=True)
for ax, values, title, cmap, limits in zip(
    axes, [data[0]["phi"], data[1]["phi"], error],
    ["Nominal timestep", "Half timestep", "Half minus nominal"],
    ["Blues", "Blues", "coolwarm"],
    [(0, 1), (0, 1), (-max(abs(error[w > 0])), max(abs(error[w > 0])))],
):
    field = np.ma.array(values.reshape(ny, nx), mask=(w == 0).reshape(ny, nx))
    palette = plt.get_cmap(cmap).copy()
    palette.set_bad("0.65")
    image = ax.imshow(field, origin="lower", extent=(0, 5, 0, 3), cmap=palette,
                      vmin=limits[0], vmax=limits[1], interpolation="nearest")
    ax.set(title=title, xlabel="x (mm)", ylabel="y (mm)")
    fig.colorbar(image, ax=ax, shrink=.7)
fig.suptitle("Prescribed-velocity transport, t = 0.02 s — NOT a two-phase CFD result")
fig.savefig(output / "phase_comparison.png", dpi=180)
shutil.copy2(__file__, output / Path(__file__).name)
print(json.dumps(report, indent=2))
