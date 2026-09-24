# Case: independent validation of the rigorous 3D Lehmann–Goerisch lower bounds

**Case ID:** lg3d_validation_20260609 · **Date:** 2026-06-09 · **Status:** finding
identified, awaiting decision on the fix.

## Problem setting

Independently validate the rigorous (INTLAB) lower eigenvalue bounds for the
Dirichlet Laplacian produced by the archived drivers
`run_regular_tet_one_element_p12_lg.m` and `run_regular_tet_ref1_p12_lg.m`
(regular tetrahedron, P12) and `run_special_tet_p12_exactrho_lg.m` (fundamental
tetrahedron `T_F`, exact eigenvalues `mu_k = (pi^2/4)|k|^2`). Goal: re-derive the
certified Lehmann–Goerisch lower bounds through an independently written code path
that assembles **every** matrix entering the bound from interval element geometry,
and check them against the published numbers and (on `T_F`) the exact eigenvalues.

## What was checked

1. **Math audit (manuscript §4).** Bound formula, shift bookkeeping, index map
   (`nu_{q+1-j}`), and the B-positive-definiteness route
   (`rho_sh > Lambda_n ⇒ B PD`, manuscript Lemma) all reproduced correctly. The
   archived code certifies B PD only at the matrix **midpoint**
   (`min(eig(mid(B_lg)))`); the new `lg_certify.m` uses verified `isspd` + the
   Lemma margin.
2. **Interval-primitive smoke test.** In `INTERVAL_MODE=1` all VFEM primitives
   return genuine `intval`. The RT Piola transform `GetElementTransMat` has a
   genuine enclosure radius **~1e-8** at p=4 (midpoint matches double to ~1e-10),
   i.e. far larger than the ~1e-16 of a double-wrapped matrix.
3. **Rigorous re-implementation** (`../Code/`): interval `A0,A1`
   (`assemble_cg_lg_small_matrices_iv.m`), interval `A2` + verified `w2`
   (`interval_lg_rt_terms.m`, `build_interval_rt_matrices.m`), certified bounds
   (`lg_certify.m`), and certified **upper** bounds via `veig(A0−λ̂A1, A1)`.

## Result

The certified **upper** bounds reproduce the published Ritz values exactly
(one-element P12: 150.9727216782, 278.6982779029). The certified **lower** bounds
**could not be computed**: the interval Goerisch matrix `A2` blows up.

### Root cause (confirmed)

The approximate RT auxiliary field `w1`, obtained from the single-element RT
mixed saddle-point solve, has a **huge Euclidean norm** while its A_rt-norm is
tiny, because the high-order RT **mass** matrix `A_rt` is extremely
ill-conditioned (numerically singular at high p). Measured on the regular tet:

| p  | DimRT | ‖w1‖₂ | min eig(A_rt) | max eig(A_rt) |
|----|------:|------:|--------------:|--------------:|
| 6  | 280   | 4.07e2 | 3.24e-11 | 1.65e-2 |
| 10 | 924   | 2.10e4 | 1.99e-17 | 8.09e-2 |

The Goerisch term `A2[k,k] = w1ₖ'·A_rt·w1ₖ + λ̂·w2ₖ'·M_dg·w2ₖ` has a true value
~1e-3, but a **genuine** interval `A_rt` (entrywise enclosure radius ~1e-4 at p=6,
from the verified Piola transform) gives `rad(w1'·A_rt·w1) ≈ ‖w1‖₂²·rad(A_rt)`,
which swamps the value:

| p=6 | mid | rad | rad/\|mid\| |
|----|----:|----:|-----------:|
| k=1 | 1.238e-3 | 1.228e-2 | **9.92** |
| k=2 | 2.134e-3 | 7.240e-3 | 3.39 |
| k=3 | 2.826e-3 | 3.449e-3 | 1.22 |
| k=4 | 2.728e-3 | 2.751e-3 | 1.01 |

At p=12 (‖w1‖₂ ≳ 1e5) the constraint-residual / `A2` radius reaches ~8e9 and the
small EVP `B` cannot be certified PD → no lower bound.

### Why the archived code "worked"

The archived drivers form `A2` and the constraint with `intval(double(A_rt))`,
`intval(double(B_rt))` — **thin** wraps (radius ~1e-16) that do **not** enclose the
exact integrals on the irrational regular-tet geometry, and that hide the genuine
~1e-4–1e-8 enclosure width. Their saved `A2_iv` has radius only 1.1e-7 and the
published bounds (e.g. one-element λ₁ ≥ 150.9706556239) follow. This is exactly the
`intval(double(...))` construction the project note correctly forbids for `M_dg` —
but here applied to `A_rt` and `B_rt`, where (because ‖w1‖₂ is huge) it is **not**
a negligible shortcut: it materially affects whether a certified bound exists.

## Interpretation

The published one-/8-element 3D LG **lower** bounds are rigorous only under the
unstated assumption that `A_rt`, `B_rt` may be thin-wrapped. Genuine interval
assembly of those matrices (as rigor requires) makes `A2` uncomputable for the
current huge-norm `w1`. The certified **upper** bounds and the `M_dg`/`A0`/`A1`
treatment are sound.

The huge `‖w1‖₂` is a numerical artifact of the near-singular high-order RT
saddle point, **not** intrinsic: the auxiliary flux only needs `div(w1) ≈ −v`, and
its component in the small-eigenvalue subspace of `A_rt` does not change the true
`A2` value but poisons the interval. A small-A_rt-norm / small-Euclidean-norm `w1`
makes the rigorous interval `A2` tight.

## Options for the fix (decision needed)

- **(A) Small-norm `w1`.** Compute `w1` by minimizing its A_rt-norm subject to the
  divergence constraint (the Lehmann minimization the manuscript already states),
  e.g. via a regularized / orthonormalized-RT-basis solve, so `‖w1‖₂` is O(1).
  Then the rigorous interval `A2` is tight and the bound is genuinely certified.
- **(B) Certified thin-wrap correction.** Keep the efficient thin `A_rt`,`B_rt`
  but add an explicit, rigorously bounded correction for the
  `(double − exact)` matrix error contracted against `w1` — a perturbation term.
- **(C) Reformulate A2** in a well-conditioned (orthonormal) RT basis so `w1` is
  moderate by construction.

## Reproduce

```
matlab -batch "diag_qform(6)"      % ‖w1‖, A_rt spectrum, q-form value vs rad
matlab -batch "diag_w2"            % verifylss is fine on the real rhs (res 4.8e-10)
matlab -batch "run_regular_tet_one_element_lg_rigorous(12)"   % upper ok, lower NaN
```
Logs in `../Code/logs/`. Diagnostics `diag_qform.m`, `diag_w2.m`, `diag_brt.m`,
`smoke_test.m`. Archived `A2_iv` radius (1.1e-7) confirmed via
`../code_not_rigorous/diag_archived_a2.m`.

## Acceptance criteria

- [x] Rigorous interval assembly of all matrices entering the bound implemented.
- [x] Certified upper bounds reproduce published Ritz values.
- [ ] Certified lower bounds reproduced — **blocked** on the `w1` conditioning
  decision above.
- [ ] `T_F` exact-eigenvalue bracketing check — pending the same fix.
