# Rigorous 3D Lehmann–Goerisch eigenvalue bounds (validated re-implementation)

This folder is the **rigorous** re-implementation produced during the independent
validation of 2026-06-09. The previous code is archived in `../code_not_rigorous/`
(see `../code_not_rigorous/README_ARCHIVE.md`) and the plan in
`../code_not_rigorous/lg3d_validation_plan.html`.

## What changed (rigor fixes)

The archived LG drivers formed the small Lehmann–Goerisch matrices `A0, A1, A2`
from **floating-point** element matrices wrapped in `intval`. On domains with
irrational geometry (the regular tetrahedron) `intval(double(...))` encloses only
the *rounded* matrix, not the exact integral. The smoke test confirmed the RT
Piola transform (`GetElementTransMat`) carries a genuine enclosure radius ~1e-8,
i.e. far larger than the ~1e-16 of a double-wrapped matrix.

This re-implementation assembles **every** matrix entering the bound from interval
element geometry via the VFEM Bernstein primitives, element by element:

| Matrix | Archived (non-rigorous) | Here (rigorous) |
|---|---|---|
| `A0, A1` (CG energy/mass) | `intval(double(K,M))` | `assemble_cg_lg_small_matrices_iv.m` |
| `A2` RT term `w1'·A_rt·w1` | `w1'·intval(double(A_rt))·w1` | `interval_lg_rt_terms.m` |
| `A2` DG term, `w2` solve | interval `M_dg`, verified solve | kept (was already rigorous) |
| `B` positive definiteness | midpoint `eig` only | `isspd` + Lemma margin (`lg_certify.m`) |
| upper bounds | double Ritz | certified `veig(A0-λ̂A1, A1)` |

The approximate RT field `w1` is still produced by a double saddle-point solve —
the theorem permits `w1` to be any fixed field; only the constraint
`B_rt·w1 − λ̂·M_dg·w2 + M_dg·v = 0` must hold, and `w2` is solved rigorously so the
true (enclosed) `w2` satisfies it exactly.

## Files

Helpers:
- `assemble_cg_lg_small_matrices_iv.m` — interval `A0, A1` for conforming CG^p.
- `interval_lg_rt_terms.m` — interval `A2` (`= W1'·A_rt·W1 + λ̂·W2'·M_dg·W2`) and
  verified `w2`, element-wise (no global interval sparse matrix).
- `build_interval_dg_mass_matrix.m` — interval DG mass (kept, already rigorous).
- `lg_certify.m` — certified lower bounds for an eigenvalue group, with verified
  B-positive-definiteness.
- `local_register_cg_dof.m` — conforming CG^p DOF registration.
- `build_scalar_rt_matrices.m` — double RT matrices, used only to produce `w1`.
- `run_conforming_lg_rigorous.m` — shared CG+RT engine for conforming meshes.

Drivers:
- `run_regular_tet_one_element_lg_rigorous.m` (`p`) — regular tet, one bubble element.
- `run_regular_tet_ref1_lg_rigorous.m` (`p`) — regular tet, 8 tetrahedra.
- `run_special_tet_lg_rigorous.m` (`p`, `n_refine`) — fundamental tet `T_F`,
  bounds checked against the EXACT eigenvalues `mu_k = (pi^2/4)|k|^2`.

## Reproduce

```matlab
% from this folder, with MATLAB + INTLAB available
run_regular_tet_one_element_lg_rigorous(12)
run_regular_tet_ref1_lg_rigorous(12)
run_special_tet_lg_rigorous(12, 0)   % T_F exact-eigenvalue falsification check
```

Logs are written to `logs/`; case reports live in `reports/`.
