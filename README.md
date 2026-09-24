# LG-Method — Julia scripts for the two-stage Lehmann–Goerisch eigenvalue bounds

Companion code for

> Xuefeng Liu and Michael Plum,
> *A Two-Stage Finite Element Approach for High-precision Guaranteed Lower Eigenvalue Bounds*.

The paper proposes a two-stage algorithm for guaranteed eigenvalue bounds of
self-adjoint differential operators: a projection-based (Crouzeix–Raviart)
lower bound supplies the rough spectral separation, and the Lehmann–Goerisch
theorem with high-order conforming finite elements turns it into sharp certified
lower bounds.

This repository contains

* the **MATLAB/INTLAB implementation** with which the results in the paper were
  produced (`matlab/`, see `matlab/README.md`): the rigorous three-dimensional
  computations of Sections 5.3 and 5.4 with logs and saved results, the
  Crouzeix–Raviart prebound computation, and the 2D dumbbell example of
  Section 5.1;
* a **Julia version** of the three-dimensional computations (Sections 5.3 and
  5.4), built on the verified finite element library
  [`VFEM.jl`](https://github.com/xfliu/vfem) and the verified eigenvalue solver
  [`Veigs.jl`](https://github.com/xfliu/veigs).

The computations can also be run online at <https://ganjin.online/xfliu/LG-Method>.

## Contents

| Path | Purpose |
| --- | --- |
| `src/LGPaper.jl` | Package `LGPaper`: meshes of the two test tetrahedra, Stage 1 (`cr_projection_bounds`) and Stage 2 (`lg_bounds`, Algorithm 1 with the shift of Remark 4.7 and eigenvalue groups) on top of `VFEM.jl` |
| `scripts/run_special_tetrahedron.jl` | Section 5.3, Tables 5 and 6: the special tetrahedron `T_F` with analytic eigenvalues, `CG^p + RT^p` for `p = 2, 3, 6` |
| `scripts/run_regular_tetrahedron.jl` | Section 5.4, Table 7: the unit regular tetrahedron, `CG^8 + RT^8` on 8 tetrahedra; also recomputes the Stage-1 CR bounds on coarse refinements |
| `scripts/cr_efficiency_table.jl` | Section 5.4, Table 8: predicted cost of a direct CR computation |
| `scripts/special_tet_crosscheck.jl` | Independent cross-check of the CG Rayleigh–Ritz stage and the CR stage on `T_F` against the values used in the paper |
| `results/` | Output of the scripts as run on the authors' server (Julia 1.12.6) |
| `matlab/` | The MATLAB/INTLAB reference implementation (3D rigorous code, CR prebounds, 2D dumbbell example) — see `matlab/README.md` |

Each `run_*` script compares its output with the numbers reported in the paper
(columns `lb-paper`, `ub-paper`).

## Installation

`VFEM.jl` and `Veigs.jl` are not in the Julia General registry. `Manifest.toml`
pins the exact commits of both that were used for `results/`.

```bash
git clone https://github.com/xfliu/LG-Method.git
cd LG-Method
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

If instantiation from the manifest fails (e.g. on a newer Julia), add the two
dependencies explicitly:

```julia
using Pkg
Pkg.add(url = "https://github.com/xfliu/veigs", subdir = "VEIGS.jl")
Pkg.add(url = "https://github.com/xfliu/vfem")
Pkg.instantiate()
```

Requires Julia 1.10 or newer. All examples run on a laptop; the largest case
(`CG^6 + RT^6` on 64 tetrahedra, 15232 RT degrees of freedom) needs a few minutes.

## Usage

```bash
julia --project=. scripts/run_special_tetrahedron.jl          # Tables 5 and 6 (interval mode)
julia --project=. scripts/run_special_tetrahedron.jl 6 2      # one configuration (p = 6, level 2)
julia --project=. scripts/run_special_tetrahedron.jl --float  # floating point only
julia --project=. scripts/run_regular_tetrahedron.jl          # Table 7
julia --project=. scripts/cr_efficiency_table.jl              # Table 8
julia --project=. scripts/special_tet_crosscheck.jl           # CG/CR cross-check
```

From Julia:

```julia
using LGPaper
m = fundamental_tetrahedron_mesh(2)                    # T_F, 64 tetrahedra
ρ2, ρ4 = FUNDAMENTAL_TET_EXACT[2], FUNDAMENTAL_TET_EXACT[4]
r = lg_bounds(m, 3, 3, [1:1 => ρ2, 2:3 => ρ4]; lhat = 200.0, interval = true)
print_result(r; exact = FUNDAMENTAL_TET_EXACT[1:3])
```

`lg_bounds(m, p, neig, groups; lhat, RT_order, interval)` bounds the first
`neig` Dirichlet eigenvalues; each entry `idx => rho` of `groups` bounds the
eigenvalues with indices `idx` by one Lehmann–Goerisch pencil with the a-priori
separator `rho ≤ λ_{max(idx)+1}`.

## What is reproduced, and how

* **Floating-point mode** (`interval = false`) follows Algorithm 1 exactly:
  Rayleigh–Ritz eigenpairs in `CG^p`, auxiliary fields `(w¹, w²)` from the
  shifted mixed system in `RT^p × DG^p`, the Goerisch matrices
  `A0, A1, A2` of the shifted operator `−Δ + λ̂`, and the lower bound
  `λ_k ≥ ρ − ρ/(1 − ν)` with the shift removed afterwards.
* **Interval mode** (`interval = true`): `A0`, `A1` and the DG mass matrix are
  assembled by `VFEM.jl` in interval arithmetic, `w²` is enclosed by a verified
  element-wise solve of the constraint `div w¹ + λ̂ w² + v = 0` (Krawczyk
  iteration), the positive definiteness of `B` and the pencil eigenvalues `ν`
  are certified with `Veigs.jl`, and the upper bounds are certified enclosures
  of the projected Rayleigh–Ritz pencil.

### Comparison with the values reported in the paper

The certified values in the paper were produced with the authors' MATLAB/INTLAB
implementation, in which every matrix, including the RT mass and divergence
matrices, is assembled from interval element geometry. `VFEM.jl` currently
assembles the RT matrices in floating point only, so in interval mode they are
wrapped into thin intervals. The consequences, as recorded in `results/`:

| Quantity | Julia vs. paper |
| --- | --- |
| Upper bounds (all cases) | agree to `1e-10` – `4e-8` |
| Lower bounds, `CG^2`, `CG^3` on `T_F` | Julia larger by `1e-8` – `4e-6` |
| Lower bounds, `CG^6` on `T_F` | Julia larger by `4e-5` – `3e-2` |
| Lower bounds, `CG^8` on the regular tetrahedron | Julia larger by `2.4e-3` (λ₁) and `0.39` (λ₂–λ₄) |

The paper's lower bounds are the smaller (more conservative) ones because the
INTLAB computation carries the enclosure radius of the RT integration, which
grows quickly with the polynomial degree (the quantity `rad(A2)` discussed in
the paper). The Julia interval mode is therefore *verified except for the
rounding in the assembly of the RT matrices*; a reader who needs the full rigor
of the paper for the RT term should use the paper's values. In addition,
`IntervalArithmetic.jl` reports its `isguaranteed` flag as `false` for the
results, because `VFEM.jl` converts floating-point differences of node
coordinates into intervals during assembly (`print_result` shows this flag).

Two further conventions shared with the paper's computation:

* The vertices `(0.5, √3/2, 0)` and `(0.5, √3/6, √6/3)` of the regular
  tetrahedron are rounded to double precision.
* The Stage-1 separators of Section 5.4 (`ρ₂ = 264.5223894731`,
  `ρ₅ = 386.4419682205`) were obtained in the paper from a CR computation on a
  mesh with 127661 tetrahedra; `run_regular_tetrahedron.jl` records these
  values and recomputes the CR projection bounds on the uniform refinements of
  the regular tetrahedron (up to 4096 tetrahedra) to illustrate the stage.
  Note that `VFEM.jl`'s own `cr_liu_lower_bounds_3d` uses the constant
  `1/√10` for `C_h`; `LGPaper` uses `0.3804 h_max` as in the paper.

Run times on the authors' server (Julia 1.12.6, two Xeon Gold 5318Y): every
configuration of Tables 5–7 takes between one and 25 seconds in interval mode
after compilation.

## Not covered by the Julia code

The two-dimensional examples of Sections 5.1 and 5.2 (Dirichlet Laplacian on the
dumbbell domain, Steklov eigenvalues on the square and the L-shaped domain) were
computed in MATLAB; the dumbbell code is in `matlab/2d_dumbbell/`, while the
driver scripts of the Steklov runs are not included. `VFEM.jl` has no Steklov
eigenvalue support at present.

## Citation

```bibtex
@article{LiuPlum-two-stage,
  author  = {Liu, Xuefeng and Plum, Michael},
  title   = {A Two-Stage Finite Element Approach for High-precision Guaranteed Lower Eigenvalue Bounds},
  year    = {2026}
}
```

## License

MIT, see `LICENSE`.
