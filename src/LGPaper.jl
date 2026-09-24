"""
    LGPaper

Julia (VFEM.jl) implementation of the two-stage algorithm of

    X. Liu and M. Plum, "A Two-Stage Finite Element Approach for
    High-precision Guaranteed Lower Eigenvalue Bounds",

for the three-dimensional Dirichlet Laplacian on tetrahedral meshes.

* Stage 1 (`cr_projection_bounds`): rough lower bounds from the
  Crouzeix–Raviart projection estimate (Theorem 3.1 of the paper)
      λ_k ≥ λ_{k,h}^{CR} / (1 + C_h² λ_{k,h}^{CR}),   C_h = 0.3804 h_max.
* Stage 2 (`lg_bounds`): Lehmann–Goerisch bounds (Theorem 4.4 / Algorithm 1)
  with conforming CG^p trial functions, scalar RT^p auxiliary fields, the shift
  technique of Remark 4.7 (shift λ̂) and eigenvalue groups that share one
  separator ρ.

Two modes are provided:

* `interval = false` — plain floating point (fast; for exploration).
* `interval = true`  — the small Lehmann–Goerisch matrices A0, A1 (CG) and the
  DG mass matrix are assembled in interval arithmetic by VFEM.jl, the auxiliary
  field component w² is enclosed by a verified element-wise solve of the
  constraint equation, and the matrix pencils are certified with Veigs.jl.
  The RT mass and divergence matrices of VFEM.jl are available in floating
  point only and are wrapped into thin intervals; see the README for the
  difference to the MATLAB/INTLAB computation reported in the paper.
"""
module LGPaper

using LinearAlgebra, SparseArrays, Printf
using IntervalArithmetic
using IntervalArithmetic: interval, inf, sup, mid, isguaranteed
using VFEM
import Veigs

export mesh_from_nodes_elements, regular_tetrahedron_mesh, fundamental_tetrahedron_mesh,
       cr_projection_bounds, lg_bounds, LGResult, print_result, LIU_CR_CONST_3D,
       FUNDAMENTAL_TET_EXACT

"C_h = C · h_max for the Crouzeix–Raviart interpolation on tetrahedra, as used in the paper."
const LIU_CR_CONST_3D = 0.3804

"Exact Dirichlet eigenvalues λ_1..λ_5 of the fundamental tetrahedron T_F: (π²/4)·|k|²."
const FUNDAMENTAL_TET_EXACT = (pi^2 / 4) .* [80.0, 140.0, 140.0, 160.0, 208.0]

# ---------------------------------------------------------------------------
# Meshes
# ---------------------------------------------------------------------------

"""
    mesh_from_nodes_elements(nodes, elements) -> Mesh3D

Build a VFEM.jl `Mesh3D` from an `n×3` node table and an `m×4` connectivity table.
"""
function mesh_from_nodes_elements(nodes::AbstractMatrix{<:Real}, elements::AbstractMatrix{<:Integer})
    nodes = Matrix{Float64}(nodes)
    elements = sort(Matrix{Int}(elements), dims = 2)
    facets = get_facet_list(elements)
    edges = get_edge_list(elements)
    f2e, e2f = facet_element_connectivity(elements, facets)
    return Mesh3D(nodes, elements, facets, edges, f2e, e2f,
                  size(nodes, 1), size(elements, 1), size(facets, 1), size(edges, 1))
end

"""
    regular_tetrahedron_mesh(level) -> Mesh3D

Unit regular tetrahedron (all edges of length one), red-refined `level` times
(8^level tetrahedra). Section 5.4 of the paper uses `level = 1`.
"""
function regular_tetrahedron_mesh(level::Integer)
    nodes = [0.0 0.0 0.0;
             1.0 0.0 0.0;
             0.5 sqrt(3) / 2 0.0;
             0.5 sqrt(3) / 6 sqrt(6) / 3]
    m = mesh_from_nodes_elements(nodes, reshape([1, 2, 3, 4], 1, 4))
    for _ in 1:level
        m = red_refine_mesh_3d(m)
    end
    return m
end

"""
    fundamental_tetrahedron_mesh(level) -> Mesh3D

The special tetrahedron T_F = conv{(0,0,0), (0,0,1), (1/2,1/2,1/2), (-1/2,1/2,1/2)}
of Section 5.3, red-refined `level` times.
"""
fundamental_tetrahedron_mesh(level::Integer) = special_tetrahedron_red_mesh(level)

# ---------------------------------------------------------------------------
# Stage 1: Crouzeix–Raviart projection bounds
# ---------------------------------------------------------------------------

"""
    cr_projection_bounds(m, neig; C = LIU_CR_CONST_3D)

Lower bounds λ_k ≥ λ_{k,h}/(1 + C_h² λ_{k,h}), C_h = C·h_max, from the
Crouzeix–Raviart eigenvalues on `m` (floating point). Returns a named tuple
`(lower, cr_eig, Ch, hmax, ndof)`.
"""
function cr_projection_bounds(m::Mesh3D, neig::Integer; C::Real = LIU_CR_CONST_3D)
    Mcr, Acr, info = create_matrix_crouzeix_raviart_3d(m)
    int = info.interior_dofs
    length(int) ≥ neig || error("not enough interior CR facets ($(length(int))) for $neig eigenvalues")
    vals = eigen(Symmetric(Matrix(Acr[int, int])), Symmetric(Matrix(Mcr[int, int]))).values
    vals = sort(real.(vals))[1:neig]
    hmax = find_mesh_hmax_3d(m)
    Ch = C * hmax
    lows = vals ./ (1 .+ Ch^2 .* vals)
    return (lower = lows, cr_eig = vals, Ch = Ch, hmax = hmax, ndof = length(int))
end

# ---------------------------------------------------------------------------
# Helpers for the interval mode
# ---------------------------------------------------------------------------

_thin(S::SparseMatrixCSC{Float64}) =
    SparseMatrixCSC(S.m, S.n, copy(S.colptr), copy(S.rowval), interval.(S.nzval))
_thin(A::AbstractMatrix{Float64}) = interval.(A)

_inside(a, b) = inf(b) < inf(a) && sup(a) < sup(b)
_inflate(x) = x * interval(1 - 1e-10, 1 + 1e-10) + interval(-1e-300, 1e-300)

"""
    verified_solve(A, B)

Enclosure of the solution `X` of `A X = B` for a small dense interval matrix `A`
(Krawczyk iteration with ε-inflation). Throws if no inclusion is obtained.
"""
function verified_solve(A::AbstractMatrix{<:Interval}, B::AbstractMatrix{<:Interval})
    n = size(A, 1)
    Am = mid.(A)
    R = interval.(inv(Am))
    Xa = interval.(Am \ mid.(B))
    Z = R * (B - A * Xa)
    C = interval.(Matrix{Float64}(I, n, n)) - R * A
    X = Z
    for _ in 1:30
        Xi = _inflate.(X)
        Y = Z + C * Xi
        if all(_inside.(Y, Xi))
            return Xa + Y
        end
        X = Y
    end
    error("verified_solve: Krawczyk iteration did not contract (matrix too ill-conditioned)")
end

# ---------------------------------------------------------------------------
# Stage 2: Lehmann–Goerisch bounds
# ---------------------------------------------------------------------------

"""
    LGResult

Result of `lg_bounds`. `lb[k]`/`ub[k]` are the lower/upper bounds of λ_k that
are certified in interval mode (inf/sup of the enclosures) or computed in
floating point otherwise. `groups` records, per eigenvalue group, the separator
ρ, the Lehmann eigenvalues ν, whether B was verified positive definite and the
margin `inf(ρ+λ̂) − sup(Λ_max)` of Lemma 4.5 (positive ⇒ B positive definite).
"""
struct LGResult
    p::Int
    RT_order::Int
    lhat::Float64
    interval::Bool
    NumElt::Int
    DimCG::Int
    DimInt::Int
    DimRT::Int
    DimDG::Int
    ritz::Vector{Float64}
    ub::Vector{Float64}
    lb::Vector{Float64}
    groups::Vector{NamedTuple}
    a2_check::Vector{Float64}
    constraint_residual::Float64
    guaranteed::Bool
end

# A0u is the Goerisch matrix of the *unshifted* operator; the pencil is formed for
# the shifted operator −Δ + λ̂ with A0 = A0u + λ̂ A1 and ρ_sh = ρ + λ̂ (Remark 4.7).
function _certify_float(A0u, A1, A2, idx, rho, lhat)
    rho_sh = rho + lhat
    A0us = Symmetric((A0u[idx, idx] + A0u[idx, idx]') / 2)
    A1s = Symmetric((A1[idx, idx] + A1[idx, idx]') / 2)
    A2s = Symmetric((A2[idx, idx] + A2[idx, idx]') / 2)
    A0s = Symmetric(A0us + lhat * A1s)
    AL = Symmetric(A0s - rho_sh * A1s)
    BL = Symmetric(A0s - 2rho_sh * A1s + rho_sh^2 * A2s)
    q = length(idx)
    Λ = eigvals(A0us, A1s)
    margin = rho - maximum(Λ)          # Lemma 4.5: ρ > Λ_n  ⇒  B positive definite
    B_spd = isposdef(BL)
    nu = B_spd ? sort(eigvals(AL, BL)) : fill(NaN, q)
    lb = [B_spd ? rho_sh - rho_sh / (1 - nu[q + 1 - j]) - lhat : NaN for j in 1:q]
    return (idx = idx, rho = rho, nu = nu, B_spd = B_spd, margin = margin, lb = lb)
end

function _certify_interval(A0u, A1, A2, idx, rho, lhat)
    rho_i = interval(Float64(rho)); lhat_i = interval(Float64(lhat))
    rho_sh = rho_i + lhat_i
    A0us = Veigs.sym_hull(A0u[idx, idx]); A1s = Veigs.sym_hull(A1[idx, idx]); A2s = Veigs.sym_hull(A2[idx, idx])
    A0s = Veigs.sym_hull(A0us + lhat_i * A1s)
    AL = Veigs.sym_hull(A0s - rho_sh * A1s)
    BL = Veigs.sym_hull(A0s - interval(2.0) * rho_sh * A1s + rho_sh^2 * A2s)
    q = length(idx)
    B_spd = try
        Veigs.verified_isspd(BL)
    catch
        false
    end
    Λ, _ = Veigs.veig(A0us, A1s)
    margin = inf(rho_i) - sup(Λ[end])   # Lemma 4.5: ρ > Λ_n  ⇒  B positive definite
    if !B_spd
        return (idx = idx, rho = rho, nu = Interval{Float64}[], B_spd = false, margin = margin,
                lb = fill(NaN, q))
    end
    nu = q == 1 ? [AL[1, 1] / BL[1, 1]] : first(Veigs.veig(AL, BL))
    nu = sort(nu; by = inf)
    lb = Vector{Float64}(undef, q)
    for j in 1:q
        ν = nu[q + 1 - j]
        lb[j] = inf(rho_sh - rho_sh / (interval(1.0) - ν) - lhat_i)
    end
    return (idx = idx, rho = rho, nu = nu, B_spd = true, margin = margin, lb = lb)
end

"""
    lg_bounds(m, p, neig, groups; lhat = 200.0, RT_order = p, interval = true, verbose = true)

Two-sided bounds for the first `neig` Dirichlet eigenvalues on `m` from the
Lehmann–Goerisch method with CG^p trial functions and RT^`RT_order` auxiliary
fields (Algorithm 1 of the paper, with shift `lhat`).

`groups` is a vector of pairs `idx => rho`: the eigenvalues with indices `idx`
(a subset of `1:neig`, e.g. `2:4`) are bounded from below by one Lehmann–Goerisch
pencil using the a-priori separator `rho ≤ λ_{max(idx)+1}`.
"""
function lg_bounds(m::Mesh3D, p::Integer, neig::Integer, groups;
                   lhat::Real = 200.0, RT_order::Integer = p,
                   interval::Bool = true, verbose::Bool = true)
    lhat = Float64(lhat)
    t0 = time()
    # ---- conforming CG^p Rayleigh–Ritz problem (floating point) ----
    A, M, info, L2G = create_matrix_lagrange_3d(m, p)
    int = info.interior_dofs
    DimInt = length(int)
    DimInt ≥ neig || error("only $DimInt interior CG DOFs for $neig eigenvalues")
    Kint = Symmetric(Matrix(A[int, int])); Mint = Symmetric(Matrix(M[int, int]))
    F = eigen(Kint, Mint)
    ritz = F.values[1:neig]
    V = zeros(info.DimCG, neig)
    for k in 1:neig
        v = F.vectors[:, k]
        v ./= sqrt(dot(v, Mint * v))
        V[int, k] .= v
    end
    verbose && @printf("[CG^%d] %d tets, DimCG=%d, interior=%d  (%.1fs)\n", p, m.NumElt, info.DimCG, DimInt, time() - t0)
    verbose && (print("[CG^$p] Ritz upper bounds: "); foreach(x -> @printf("%.10f ", x), ritz); println())

    # ---- approximate auxiliary fields w = (w1, w2): shifted saddle point ----
    rt = create_matrix_rt_3d(m, RT_order)
    rt.DegK == size(L2G, 2) || error("RT_order must equal p (same local Bernstein basis for CG and DG)")
    DegK = rt.DegK
    Vdg = zeros(rt.DimDG, neig)
    for e in 1:m.NumElt
        Vdg[(e - 1) * DegK .+ (1:DegK), :] .= V[L2G[e, :], :]
    end
    SP = [rt.A_rt rt.B_rt'; rt.B_rt -lhat * rt.M_dg]
    rhs = [zeros(rt.DimRT, neig); -(rt.M_dg * Vdg)]
    Fsp = lu(SP)
    sol = Fsp \ rhs
    sol .+= Fsp \ (rhs - SP * sol)
    W1 = sol[1:rt.DimRT, :]
    W2 = sol[rt.DimRT + 1:end, :]
    res = norm(rt.B_rt * W1 - lhat * (rt.M_dg * W2) + rt.M_dg * Vdg) / norm(rt.M_dg * Vdg)
    verbose && @printf("[RT^%d] DimRT=%d, DimDG=%d, constraint residual %.2e  (%.1fs)\n",
                       RT_order, rt.DimRT, rt.DimDG, res, time() - t0)

    if !interval
        A0 = V' * (A * V); A1 = V' * (M * V)
        A2 = W1' * (rt.A_rt * W1) + lhat * (W2' * (rt.M_dg * W2))
        a2_check = [A2[k, k] * (ritz[k] + lhat) for k in 1:neig]
        grp = [_certify_float(A0, A1, A2, collect(idx), Float64(rho), lhat) for (idx, rho) in groups]
        ub = copy(ritz)
        guaranteed = false
    else
        # ---- interval Lehmann–Goerisch matrices ----
        Ai, Mi, _, _ = create_matrix_lagrange_3d(m, p; T = Interval{Float64})
        Vi = _thin(V)
        A0 = Veigs.sym_hull(Vi' * (Ai * Vi))
        A1 = Veigs.sym_hull(Vi' * (Mi * Vi))
        Mdg_i, _ = create_matrix_dg_3d(m, RT_order; T = Interval{Float64})
        Art_i = _thin(rt.A_rt)          # floating-point RT mass, thin intervals (see README)
        Brt_i = _thin(rt.B_rt)          # floating-point divergence matrix, thin intervals
        W1i = _thin(W1)
        Vdg_i = _thin(Vdg)
        lhat_i = IntervalArithmetic.interval(lhat)
        # w2 from the constraint  B_rt w1 − λ̂ M_dg w2 + M_dg v = 0  (verified, element by element)
        R = (Brt_i * W1i + Mdg_i * Vdg_i) ./ lhat_i
        W2i = Matrix{Interval{Float64}}(undef, rt.DimDG, neig)
        for e in 1:m.NumElt
            blk = (e - 1) * DegK .+ (1:DegK)
            W2i[blk, :] .= verified_solve(Matrix(Mdg_i[blk, blk]), R[blk, :])
        end
        A2 = Veigs.sym_hull(W1i' * (Art_i * W1i) + lhat_i * (W2i' * (Mdg_i * W2i)))
        a2_check = [mid(A2[k, k]) * (ritz[k] + lhat) for k in 1:neig]
        verbose && @printf("[LG] interval matrices assembled, max rad(A2) = %.2e  (%.1fs)\n",
                           maximum(sup.(A2) .- inf.(A2)) / 2, time() - t0)
        Λ, _ = Veigs.veig(A0, A1)
        ub = sup.(Λ)[1:neig]
        grp = [_certify_interval(A0, A1, A2, collect(idx), Float64(rho), lhat) for (idx, rho) in groups]
        guaranteed = all(isguaranteed, A0) && all(isguaranteed, A1) && all(isguaranteed, A2) &&
                     all(isguaranteed, Λ) && all(g -> all(isguaranteed, g.nu), grp)
    end
    lb = fill(NaN, neig)
    for g in grp
        lb[g.idx] .= g.lb
        verbose && @printf("[LG] group %s rho=%.6f: B_spd=%s, margin=%.3e\n",
                           string(g.idx), g.rho, string(g.B_spd), g.margin)
    end
    verbose && @printf("[LG] done in %.1fs\n", time() - t0)
    return LGResult(p, RT_order, lhat, interval, m.NumElt, info.DimCG, DimInt, rt.DimRT, rt.DimDG,
                    ritz, ub, lb, grp, a2_check, res, guaranteed)
end

"""
    print_result(io, r::LGResult; exact = nothing, reference = nothing)

Tabulate the bounds; optionally against exact eigenvalues and/or the reference
values reported in the paper (`reference = (lb = [...], ub = [...])`).
"""
function print_result(io::IO, r::LGResult; exact = nothing, reference = nothing, label = "")
    println(io, "Lehmann–Goerisch bounds ", label, " — CG^", r.p, "+RT^", r.RT_order,
            ", ", r.NumElt, " tets, λ̂=", r.lhat, ", ", r.interval ? "interval mode" : "floating-point mode")
    @printf(io, "  DimCG=%d interior=%d DimRT=%d DimDG=%d\n", r.DimCG, r.DimInt, r.DimRT, r.DimDG)
    print(io, "  k       lower bound            upper bound")
    exact === nothing || print(io, "            exact              lb<exact<ub")
    reference === nothing || print(io, "     lb-paper      ub-paper")
    println(io)
    for k in eachindex(r.lb)
        @printf(io, "  %d  %20.12f  %20.12f", k, r.lb[k], r.ub[k])
        if exact !== nothing
            ok = r.lb[k] < exact[k] < r.ub[k]
            @printf(io, "  %20.12f  %s", exact[k], ok ? "yes" : "NO ")
        end
        if reference !== nothing
            @printf(io, "  %+.3e  %+.3e", r.lb[k] - reference.lb[k], r.ub[k] - reference.ub[k])
        end
        println(io)
    end
    for g in r.groups
        @printf(io, "  group %-9s rho=%.10f  B positive definite: %s  margin(Lemma 4.5)=%.3e\n",
                string(g.idx), g.rho, string(g.B_spd), g.margin)
    end
    r.interval && println(io, "  all interval flags guaranteed: ", r.guaranteed)
end
print_result(r::LGResult; kw...) = print_result(stdout, r; kw...)

end # module
