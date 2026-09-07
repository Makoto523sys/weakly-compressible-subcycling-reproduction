# External reference data: Figure 5

Do not populate this directory with synthetic or hand-guessed paper results.

The three CSV files were extracted from the visible SVG paths of Figure 5:
https://arxiv.org/html/2608.23110v1/rising_bubble_velocity.svg
They contain the plotted polyline vertices, not the authors' original simulation
time series. Column units are seconds and metres/second. Black, red and blue
paths correspond respectively to standard, proposed and weak solvers; the legend
and axes were visually checked. `figure5_paths.json` preserves the path coordinates
and axis calibration for audit. The TOML sidecars contain hashes and provenance.

The axis spans 0–0.07 s. SVG coordinates map to physical coordinates as
`t=(x-75.6)/5760` and `v=(y-87)/660`. The path's y coordinates precede the SVG
vertical flip. A last path vertex can lie just outside the plot box; preserve it
for interpolation but only compare inside 0–0.07 s. A conservative working
uncertainty of 0.0005 m/s accounts for plotted-curve simplification; it is an
analyst choice, not a published uncertainty bound. Original data would be better.

Sections 4.2 and 4.4 are a water/air-like, high-density-ratio problem.
Do not replace Section 4.2 with the Hysing benchmark: its properties and domain
are different. A match between two solvers in this package is an internal
consistency test, not independent evidence of reproduction.
