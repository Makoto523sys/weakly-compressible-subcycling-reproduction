2026-09-08: Section 4.2 review completed with documented limits; see [final review](BUBBLE_REVIEW_JA.md). Historical verification sequence follows.

# Validation sequence

Current result: **local execution underway; reproduction not established**.
See [2026-09-07 measurements and logs](LOCAL_VALIDATION_JA.md). The sequence and
initial acceptance thresholds below remain applicable.

## 1. Code verification

Run `julia --project=. test/runtests.jl` and retain the complete log.
The tests cover the capillary time scale, affine WENO reconstruction, closed-box
flux conservation, variable-density projection, a known Neumann eigenmode,
single-phase hydrostatic balance including walls, curvature sign/magnitude and
single-substep equivalence. The original 21 tests and 28 additional pressure,
substep-rounding and checkpoint assertions now pass locally.

Before the full bubble campaign, add and run a static two-phase bubble study:
check the 2D Laplace jump sigma/R, spurious currents and mesh sensitivity. The
existing curvature test alone is insufficient to verify discrete force balance.

## 2. Short bubble runs

Use a small grid and short duration only to detect crashes, nonfinite values,
pressure-solver stagnation and boundary/phase errors. This is a smoke test;
coarse cells substantially broaden the ACDI interface relative to bubble size.

## 3. Figure 5 comparison

Use the physical case in `configs/bubble-paper.toml`. Run standard, EPP-only and
subcycling (30 capillary steps) to 0.07 s at 256×512. All start at zero velocity
and zero pressure; do not hydrostatically preinitialize only one mode, since the
paper intentionally tests the acoustic response to sudden pressure adjustment.

Plot all three computed histories alongside their independently extracted
Figure 5 curves. The reference files contain simplified SVG polylines, not raw
author time histories. Record interpolation uncertainty and avoid extrapolation.
Check the weak-only acoustic oscillations, and whether the subcycled result
tracks the standard solver without comparable oscillations. Compare on the
same output-time grid; do not infer absence of ringing from sparsely saved data.

The comparison script applies the following **analyst-defined initial criteria**,
not tolerances claimed by the authors:

| Quantity | Initial criterion |
|---|---|
| Time-weighted RMS velocity error / reference maximum | ≤2%, plus reference extraction allowance |
| Maximum velocity error / reference maximum | ≤5%, plus reference extraction allowance |
| Relative gas-area drift | ≤1e-4 over the full history |
| Raw phase bounds | within [-1e-3, 1+1e-3] |
| Maximum face divergence × dx / reference speed | ≤1e-5 after projection |

Passing these gives `velocity_checks_passed_requires_review`, never an automatic
claim of reproduction. Inspect centroid history, shape, early-time oscillations,
phase tails, pressure residuals, and the choices listed in IMPLEMENTATION.md.
Do not tune physical properties or reference uncertainty to force a pass.

## 4. Solution verification

Refine grid and time step independently. Suggested grids are 128×256, 256×512,
and 512×1024, subject to measured compute cost. Because epsilon=dx, changing the
grid also changes interface thickness. Record that effect explicitly.
For temporal sensitivity at fixed 256×512, run standard factors 1 and 0.5, and
subcycling factors 30, 15 and 7.5. Require the differences to be comfortably
below the accepted paper comparison error, or report that grid/time uncertainty
prevents a conclusive match. Tighten pressure tolerance separately when needed.
Do not publish Richardson/GCI estimates without an observed asymptotic trend.

## 5. Section 4.4 — only after bubble acceptance

Implement fluid-volume/face-aperture geometry, cut-cell conservative fluxes,
an immersed-boundary velocity treatment and a geometrical liquid contact angle
of 150 degrees. Review the cited Meyer, Mittal and Yokoi methods before coding.
Explicitly resolve how tiny cut cells, boundary curvature and pressure coefficients
are handled; a binary staircase mask or omitted contact angle is not equivalent.

Check zero solid leakage, flow balance, static contact angle and capillary entry
pressure before the full porous case. Compare the invaded path at 0.18 s and
the inlet–outlet pressure history including the sudden drop. Preserve the five
small obstacle shifts; they select the widest route. Compare matched resolutions,
error levels and hardware before reporting a speedup. Runtime measured only to
0.005 s does not replace the full 0.18 s flow-pattern comparison.

## 6. Publish only supported claims

Keep exact source revision, Julia version, configurations, input/reference hashes,
logs, metrics and plots. Mark partial or failed tests explicitly. The user's
requested order is validation first, then upload to the specified repository.
Using a pre-validation GitHub branch as a compute venue is an alternate workflow
that needs agreement before publication.
