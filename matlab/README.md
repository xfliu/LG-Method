# MATLAB / INTLAB implementation

This directory holds the MATLAB code with which the results reported in the
paper were produced. It is the reference implementation; the Julia package in
`../src` re-implements the three-dimensional part on `VFEM.jl`.

| Directory | Paper | Content |
| --- | --- | --- |
| `3d_rigorous/` | Section 5.3 (Tables 5, 6), Section 5.4 (Table 7) | Rigorous Lehmann–Goerisch bounds on tetrahedral meshes with every matrix assembled in interval arithmetic (INTLAB). `README.md` inside describes the files; `logs/` holds the raw output of the runs behind the tables, `reports/` the case reports, `*_results.mat` the saved results (MATLAB v7.3). |
| `3d_cr_prebounds/` | Section 5.4 (separators ρ₂, ρ₅; Table 8) | The Crouzeix–Raviart computation on the level-6 regular-tetrahedron mesh (127661 tetrahedra) that supplies the separators, and the post-processing for the efficiency estimate of Table 8. |
| `2d_dumbbell/` | Section 5.1 | Self-contained 2D code (Lagrange `CG^p`, `RT^p` auxiliary problem, Lehmann–Goerisch bounds) for the Dirichlet Laplacian on the dumbbell domain, with the graded meshes (`mesh_data/*.xml`, DOLFIN XML format). |

## Requirements

* MATLAB (R2023b was used) with [INTLAB](https://www.tuhh.de/ti3/rump/intlab/)
  for interval arithmetic and the verified eigenvalue solver `veig`.
* The 3D drivers call the authors' MATLAB library **VFEM3D** (Bernstein-basis
  finite element primitives on tetrahedra, `my_intlab_mode_config`,
  `get_rt_funcmat`, ...). The drivers reference it through the variable
  `vfem_root` at the top of each `run_*.m` file; set it to the location of the
  library on your machine. The library is not distributed in this repository.
* The 2D code in `2d_dumbbell/` needs only MATLAB and INTLAB.

## Reproducing the 3D tables

From `3d_rigorous/` (the drivers write `logs/<name>_<date>.log` relative to the
working directory):

```matlab
run_tf_tables()                         % Section 5.3: all T_F configurations of Tables 5 and 6
run_regular_tet_ref1_lg_rigorous(8)     % Section 5.4: Table 7 (CG^8 + RT^8, 8 tetrahedra)
```

The separators `rho2 = 264.5223894731`, `rho5 = 386.4419682205` used in
`run_regular_tet_ref1_lg_rigorous.m` were computed by
`3d_cr_prebounds/run_cr_level6_tet.m` (log in the same directory); Table 8
follows from `3d_cr_prebounds/postprocess_regular_tet_cr_efficiency.m`.

## Reproducing the dumbbell example

From `2d_dumbbell/`, edit the mesh/`rho`/`neig` block at the top of
`main_lower_bound_laplace_eig.m` (the graded meshes of the paper are in
`mesh_data/`) and run it; set `INTERVAL_MODE = 1` for the verified computation.

## Not included

The Steklov computations of Section 5.2 (unit square and L-shaped domain) used
the Steklov module of the authors' 2D MATLAB library; the driver scripts for
those runs are not part of this repository.
