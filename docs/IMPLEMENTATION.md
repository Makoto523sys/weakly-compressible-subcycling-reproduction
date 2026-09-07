2026-09-08: Section 4.2 review completed with documented limits; see [final review](BUBBLE_REVIEW_JA.md). Historical verification sequence follows.

# Implementation audit

Local execution update (2026-09-07): load and 49 tests pass on Julia 1.12.7.
Full reproduction remains unestablished. See [local measurements](LOCAL_VALIDATION_JA.md).
Optional symmetric aggregation multigrid now preconditions PCG; the fine operator
and true-residual tolerance are unchanged. Rounding-sized final EPP substeps are
absorbed into the preceding substep within 64 machine epsilons of its limit.
New runs archive their source and save exact-state, source-checked checkpoints.

Source: [paper v1, Sections 2–3 and Appendix A](https://arxiv.org/html/2608.23110v1).
All arrays use Float64; all physical inputs use SI units. The implementation is
an initial CPU reference version with many allocations, not an optimized GPU port.

## Mapping

| Source | Code |
|---|---|
| Equations 13–14, ACDI phase flux | `fluxes` |
| Equations 15–20, consistent mass/momentum flux | `fluxes`, `face_value` |
| Equations 21–28, full Newtonian stress | `fluxes` |
| Phase update and material interpolation | `increments`, `provisional!` |
| Signed-distance normal and scalar curvature | `interface_geometry` |
| Equations 32–35, localized CSF | `forces`, `add_face_increment!` |
| Equation 8, EPP | `project!(...; weak=true)` |
| Equation 9, accumulated substeps | `step!(...; mode=:subcycling)` |
| Variable-density projection | `pressure_operator!`, `poisson!`, `project!` |

Summing the divergence of each flux times its substep is algebraically equivalent
to taking the divergence of its time integral on this fixed Cartesian grid.
Surface force impulses are summed **before division by the final face density**
in the main step. The final weakly compressible velocity is not used as the main
step's provisional velocity. For a single substep the method must recover the
standard solver's result to linear-solver accuracy.

## Interpretation of apparent typographical issues

The displayed Appendix A equation 36 omits division of the right-hand side by
the time step; Section 2.2 and the correction formula require it. This code solves
`div(grad(p)/rho_face) = div(u_star)/dt`. Otherwise units and the projection
identity fail. This is a dimensional interpretation, not an author-confirmed erratum.
The displayed curvature line separates derivative terms with a comma; the code
uses their scalar sum, as required by the divergence definition.

## Choices needing numerical verification

1. The pressure backend is a zero-mean projected PCG with Jacobi preconditioning,
   rather than the paper's diagonally scaled FlexGMRES/PFMG. The default relative
   tolerance is 1e-9 on the unscaled Euclidean residual, with a true-residual check.
   This is not the same stopping norm as the paper. Maximum iterations cause an
   explicit failure, never a silent continuation.
2. `hypot(max(abs(uf)), max(abs(vf)))` is a conservative velocity magnitude bound
   used for Gamma and substep selection. It can exceed the maximum collocated
   magnitude and therefore changes ACDI regularization. Confirm the authors'
   exact face-vector norm convention or assess sensitivity before reproduction.
3. The initial liquid field is `0.5*(1+tanh((r-R)/(2*epsilon)))`, with epsilon=dx.
   This equilibrium diffuse-interface initialization is inferred from ACDI;
   Section 4.2 does not spell out its discrete initialization.
4. Momentum uses finite-volume WENO-JS3, ideal weights 1/3 and 2/3, smoothness
   exponent 2 and epsilon=1e-12 in squared velocity units. The paper identifies
   WENO3 but does not give all these choices. Upwinding uses the face velocity.
5. Closed-wall phase flux is zero. Velocity halos are reflected to impose side
   free slip and top/bottom no slip. The provisional boundary force is lifted
   into the pressure boundary condition so the cell-center and face corrections
   remain consistent. The single-phase hydrostatic test checks this construction.
6. Gravity is included as a face force `rho_face*g` and accumulated like CSF.
   The general algorithm omits gravity for simplicity; this is an extension for
   Section 4.2 and must be reviewed through hydrostatic and bubble tests.
7. Only the log used to reconstruct psi is clipped, as in the paper. Transported
   phi is never clipped. Raw overshoots/undershoots are output; nonpositive density
   or negative viscosity terminates the run.

## Practical limits

Julia loading, unit tests, static-bubble and rising-bubble calculations have now
run locally. Static surface-force imbalance remains; grid/time verification and
the 512-grid refinement review are ongoing. All three paper-grid runs completed.
Execution is not a claim of reproduction.
Published H100 speedups are not performance targets for this CPU/PCG code.

No cut-cell, immersed-boundary or contact-angle solver is included yet.
`configs/porous-paper.toml` records dimensions read from Figure 10 and lattice
centers inferred from the stated 60-degree arrangement, including the five
signed shifts. That inference needs a geometry overlay review before meshing.
