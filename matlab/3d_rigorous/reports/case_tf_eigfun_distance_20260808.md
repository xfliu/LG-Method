# Case report — why the degenerate pair $\lambda_2=\lambda_3$ on $T_F$ has unequal upper-bound errors

**Case ID:** `tf-eigfun-distance-20260808`
**Date:** 2026-08-08
**Status:** complete; question answered, no defect found in the LG pipeline
**Rigour:** floating-point diagnostic only (NOT interval). Does not affect any certified bound.

## Question

On $T_F$ the exact eigenvalue $\mu = \tfrac{\pi^2}{4}\cdot 140 = 345.4361540381$ is **double**,
yet the CG$^2$ Ritz upper bounds for $\lambda_2,\lambda_3$ at refinement level 2 have very
different errors:

| $k$ | Ritz $\lambda_{h,k}$ | exact | error |
|---|---|---|---|
| 2 | 372.4968338014 | 345.4361540381 | 27.0607 |
| 3 | 433.8690500095 | 345.4361540381 | 88.4329 |

User hypothesis: *the two eigenfunctions have different structure, so their approximation
errors differ.* Requested test: the distance from the exact 2-dimensional eigenspace $E$
to each computed eigenfunction.

**Verdict: hypothesis confirmed, and the mechanism identified exactly.**

## Problem setting

$-\Delta u = \lambda u$ on $T_F = \mathrm{conv}\{(0,0,0),(0,0,1),(\tfrac12,\tfrac12,\tfrac12),
(-\tfrac12,\tfrac12,\tfrac12)\}$, homogeneous Dirichlet BC.

Exact eigenpairs (Jia–Li–Zhang 2021, arXiv:2105.07547, Appendix *Exact eigenvalues of
homogeneous Dirichlet Laplacian on $T_F$*, eq. (g_sine)):

$$TS_{\mathbf k}(\mathbf s)=\frac1{24}\sum_{\sigma\in\mathcal G}(-1)^{|\sigma|}
e^{\frac{\pi i}{2}(\mathbf k\sigma)\cdot\mathbf s},\qquad
\mu_{\mathbf k}=\frac{\pi^2}{4}|\mathbf k|^2,$$

$\Lambda_0=\{\mathbf k\in\mathbb Z^4:\sum k_j=0,\ k_0\equiv k_1\equiv k_2\equiv k_3\ (\mathrm{mod}\ 4),\
k_0<k_1<k_2<k_3\}$, homogeneous coordinates
$s_1=\tfrac{-x_1+x_2+x_3}2,\ s_2=\tfrac{x_1-x_2+x_3}2,\ s_3=\tfrac{x_1+x_2-x_3}2,\ s_0=-(s_1+s_2+s_3)$.

### Multiplicity structure (enumeration of $\Lambda_0$)

| $\vert\mathbf k\vert^2$ | $\mu$ | \#indices | indices |
|---|---|---|---|
| 80 | 197.3920880218 | 1 | $(-6,-2,2,6)$ |
| **140** | **345.4361540381** | **2** | $(-9,-1,3,7),\ (-7,-3,1,9)$ |
| 160 | 394.7841760436 | 1 | $(-8,-4,4,8)$ |
| 208 | 513.2194288566 | 1 | $(-10,-2,2,10)$ |

The two indices at $|\mathbf k|^2=140$ are conjugate ($\mathbf k^B=-\mathrm{rev}(\mathbf k^A)$, so
$TS_{\mathbf k^B}=\overline{TS_{\mathbf k^A}}$), giving the **real** eigenspace

$$E=\mathrm{span}_{\mathbb R}\{\,\mathrm{Re}\,TS_{\mathbf k^A},\ \mathrm{Im}\,TS_{\mathbf k^A}\,\},
\qquad \mathbf k^A=(-9,-1,3,7),\quad \dim E = 2 .$$

$\lambda_1$ and $\lambda_4$ come from self-conjugate indices and are simple. This
independently confirms the multiplicities hard-coded at
`run_special_tet_lg_rigorous.m:27` (`lambda_exact = (pi^2/4)*[80;140;140;160]`).

## Discretisation

Uniform red refinement of $T_F$ (levels 2 and 3 = 64 and 512 tets), CG$^2$.
Diagnostic code is a **standalone P2-Lagrange re-implementation**
(`Code/diag_eigfun_distance.py`) of the CG space in `run_conforming_lg_rigorous.m`
(which uses the VFEM3D Bernstein basis). Quadrature: collapsed-coordinate Gauss,
$10^3$ points per tet.

**Validation gate (asserted in code):** the independently computed Ritz values must
match the MATLAB log. They do, to $3.3\times10^{-11}$ (level 2) and $2.7\times10^{-11}$
(level 3) — same discrete space, so the eigenfunctions being analysed are the ones the
LG runs actually used. Rayleigh quotients of the exact basis reproduce $\mu$ to $10^{-10}$;
$TS_{\mathbf k}$ vanishes on $\partial T_F$ to $10^{-12}$.

## Results

### Distance to the exact eigenspace, $\sin\theta_k = \|u_h^{(k)}-P_E u_h^{(k)}\|/\|u_h^{(k)}\|$

| level | mode | $\sin\theta$ ($L^2$) | $\sin\theta$ ($H^1$-semi) |
|---|---|---|---|
| 2 (64 tets) | $\lambda_1$ | 9.211879e-02 | 2.553331e-01 |
| 2 | $\lambda_2$ | 1.155619e-01 | 2.916010e-01 |
| 2 | $\lambda_3$ | **3.339311e-01** | **5.409302e-01** |
| 3 (512 tets) | $\lambda_1$ | 1.515718e-02 | 8.424525e-02 |
| 3 | $\lambda_2$ | 2.895353e-02 | 1.192213e-01 |
| 3 | $\lambda_3$ | 4.969776e-02 | 1.773688e-01 |

The third mode is genuinely **1.86× further** from $E$ than the second in the energy
seminorm at level 2 (1.49× at level 3).

### Exact eigenvalue-error identity

For any $u\in E$ with $\|u\|_0=1$: $\ \lambda_h-\lambda=|u_h-u|_1^2-\lambda\|u_h-u\|_0^2$.
Taking $u$ = normalised best approximation of $u_h$ in $E$:

| level | mode | $\vert u_h-u\vert_1^2$ | $\lambda\Vert u_h-u\Vert_0^2$ | difference | $\lambda_h-\lambda$ (from log) |
|---|---|---|---|---|---|
| 2 | $\lambda_1$ | 13.65323 | 1.67861 | 11.97462 | 11.97462 |
| 2 | $\lambda_2$ | 31.68933 | 4.62865 | 27.06068 | 27.06068 |
| 2 | $\lambda_3$ | 128.09071 | 39.65781 | 88.43290 | 88.43290 |
| 3 | $\lambda_1$ | 1.41064 | 0.04535 | 1.36528 | 1.36528 |
| 3 | $\lambda_2$ | 4.97661 | 0.28964 | 4.68697 | 4.68697 |
| 3 | $\lambda_3$ | 11.19312 | 0.85371 | 10.33941 | 10.33941 |

The identity reproduces `R.ub_err` exactly. Consequently the eigenvalue-error ratio *is*
the squared eigenfunction-error ratio:

| level | $(\sin\theta_3^{H^1}/\sin\theta_2^{H^1})^2$ | $(\lambda_{h,3}-\lambda)/(\lambda_{h,2}-\lambda)$ |
|---|---|---|
| 2 | 3.4412 | 3.2679 |
| 3 | 2.2133 | 2.2060 |

### Mechanism: residual mesh symmetry selects the two directions in $E$

$\mathrm{Isom}(T_F)$ has **order 8** ($T_F$ is a digonal disphenoid: opposite edges
$P_0P_1$ and $P_2P_3$ of length 1, the other four of length $\sqrt3/2$).

Red refinement splits the central octahedron along a fixed diagonal
(`special_tet_make_mesh.m` / `refine_tet_mesh.m`, the `m13`–`m24` choice), destroying
most of that group:

| mesh | order of the mesh-preserving subgroup |
|---|---|
| level 2 (64 tets) | **2** — identity and $\tau$ = vertex permutation $(3,2,1,0)$ |
| level 3 (512 tets) | **1** — trivial |

$\tau$ is the $180^\circ$ rotation with matrix $\begin{pmatrix}0&0&1\\0&-1&0\\1&0&0\end{pmatrix}$
($P_0\!\leftrightarrow\!P_3$, $P_1\!\leftrightarrow\!P_2$). It acts on the basis
$(\mathrm{Re}\,TS,\mathrm{Im}\,TS)$ of $E$ as $\begin{pmatrix}0&1\\1&0\end{pmatrix}$
(residual $1.5\times10^{-15}$), whose eigenvectors are $(1,\pm1)/\sqrt2$ with eigenvalues $\pm1$.
The basis is $L^2$-orthogonal with equal norms (normalised off-diagonal $2.6\times10^{-8}$),
so those eigenvectors sit at exactly $45^\circ$ and $135^\circ$.

Measured direction of $P_E u_h^{(k)}$ inside $E$, and measured parity
$\langle u_h\circ\tau,u_h\rangle/\|u_h\|^2$ at level 2:

| mode | direction in $E$ | parity under $\tau$ |
|---|---|---|
| $\lambda_1$ | — (simple) | $+1.000000000000$ |
| $\lambda_2$ | $45.000^\circ$ | $+1.000000000000$ (even) |
| $\lambda_3$ | $135.000^\circ$ | $-1.000000000000$ (odd) |

So the level-2 discrete modes are **forced** to be the symmetry-adapted pair
$(\mathrm{Re}\pm\mathrm{Im})TS_{\mathbf k^A}/\sqrt2$ — two structurally different functions.
At level 3, with no symmetry left, the directions are no longer pinned and drift to
$45.025^\circ/135.025^\circ$ (continuity), and the error asymmetry shrinks from 3.27× to 2.21×.

## Interpretation

1. The unequal errors are **not** a defect. Rayleigh–Ritz applied to a degenerate cluster
   returns $m$ discrete values that are the stationary values of the Rayleigh quotient over
   the $m$ directions of $E$, **automatically sorted by approximability**. Equal errors
   would require the FEM space to approximate $E$ isotropically, which needs the mesh to
   carry the symmetry that produces the degeneracy. Here it does not.
2. Which direction is "hard" is decided by the residual symmetry: at level 2 the
   $\tau$-even combination is well approximated ($\sin\theta_{H^1}=0.29$) and the $\tau$-odd
   one is not ($0.54$), and min–max necessarily assigns the former to $\lambda_2$ and the
   latter to $\lambda_3$.
3. The eigenvalue-error split is exactly the squared eigenfunction-error split — verified
   through an identity that reproduces `ub_err` to all printed digits.

## How to reproduce

```sh
cd /home/xfliu/Workspace/Paper_Plum_LG/Code
python3 diag_eigfun_distance.py 2>&1 | tee logs/diag_eigfun_distance_20260808.log
```

Log: `Code/logs/diag_eigfun_distance_20260808.log`.

Interactive 3D view of the four eigenfunctions (nodal surfaces / lobes):

```sh
python3 diag_eigfun_surfaces.py    # marching tetrahedra -> eigfun_surfaces.json
python3 make_discussion_html.py    # -> ../discussion_eigenfunction_lambda_2_3.html
```

A connected-component census on the sampling mesh confirms all four functions have
exactly **two** nodal domains (exact 2256/2237, `fem2` 2178/2093, `fem3` 2136/2135),
consistent with Courant's theorem for the second/third eigenvalue: the single dividing
sheet is displaced and more strongly folded for $u_h^{(3)}$, not replaced by more sheets.
Reference Ritz values are taken from
`Code/logs/special_tet_p2_ref{2,3}_lg_rigorous_20260609.log` and asserted in the script.

## Acceptance criteria

- [x] Independent code reproduces the MATLAB Ritz values (< 1e-7) — same discrete space.
- [x] Exact eigenfunctions verified: Dirichlet BC to 1e-12, Rayleigh quotient = $\mu$ to 1e-10.
- [x] Multiplicity 2 at $\mu=345.436$ derived from $\Lambda_0$, not assumed.
- [x] Distance from $E$ computed for each mode in $L^2$ and $H^1$.
- [x] Eigenvalue-error identity reproduces `R.ub_err` exactly.
- [x] Mechanism identified and verified (mesh symmetry group, parity of the modes).

## Notes / caveats / open items

- Diagnostic is plain floating point. It explains the observed Ritz errors; it makes no
  rigorous claim and does not touch any certified bound.
- The convergence rates of $\sin\theta$ (level 2 → 3) are preasymptotic and non-uniform
  ($L^2$: 2.60 / 2.00 / 2.75; $H^1$: 1.60 / 1.29 / 1.61) — expected, since the level-2
  errors are 9–33% and the modes are highly oscillatory ($|\mathbf k|$ up to 9).
- **Open (paper-relevant, user decision):** if the manuscript wants the degenerate pair to
  be approximated symmetrically, the refinement would need to preserve
  $\mathrm{Isom}(T_F)$ (e.g. a symmetry-adapted diagonal choice, or averaging over the
  three diagonals). Not attempted here — no paper edits made.
