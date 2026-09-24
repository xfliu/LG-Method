#!/usr/bin/env julia

using LinearAlgebra
using SparseArrays
using Printf
using Dates
using VFEM

const EXACT = (pi^2 / 4.0) .* [80.0, 140.0, 140.0, 160.0, 208.0]

function ijkl(M::Int)
    out = NTuple{4,Int}[]
    for i in M:-1:0
        for j in (M - i):-1:0
            for k in (M - i - j):-1:0
                push!(out, (i, j, k, M - i - j - k))
            end
        end
    end
    out
end

function special_tet_mesh(level::Int)
    nodes = [
        0.0 0.0 0.0
        0.0 0.0 1.0
        0.5 0.5 0.5
       -0.5 0.5 0.5
    ]
    tets = reshape([1, 2, 3, 4], 1, 4)
    for _ in 1:level
        nodes, tets = refine_mesh(nodes, tets)
    end
    complete_mesh(nodes, tets)
end

function refine_mesh(nodes::Matrix{Float64}, tets::Matrix{Int})
    edge_pairs = [(1,2), (1,3), (1,4), (2,3), (2,4), (3,4)]
    all_edges = Tuple{Int,Int}[]
    for e in 1:size(tets, 1)
        v = tets[e, :]
        for (q, (a, b)) in enumerate(edge_pairs)
            push!(all_edges, minmax(v[a], v[b]))
        end
    end
    unique_edges = sort!(collect(Set(all_edges)))
    edge_to_mid = Dict(edge => size(nodes, 1) + i for (i, edge) in enumerate(unique_edges))
    new_nodes = [nodes[i, :] for i in 1:size(nodes, 1)]
    for edge in unique_edges
        push!(new_nodes, 0.5 .* (nodes[edge[1], :] .+ nodes[edge[2], :]))
    end

    mid = Matrix{Int}(undef, size(tets, 1), 6)
    cursor = 0
    for e in 1:size(tets, 1)
        v = tets[e, :]
        for (q, (a, b)) in enumerate(edge_pairs)
            cursor += 1
            mid[e, q] = edge_to_mid[minmax(v[a], v[b])]
        end
    end

    refined = Matrix{Int}(undef, 8 * size(tets, 1), 4)
    r = 0
    for e in 1:size(tets, 1)
        v1, v2, v3, v4 = tets[e, :]
        m12, m13, m14, m23, m24, m34 = mid[e, :]
        children = [
            (v1, m12, m13, m14),
            (v2, m12, m23, m24),
            (v3, m13, m23, m34),
            (v4, m14, m24, m34),
            (m12, m13, m14, m24),
            (m12, m13, m23, m24),
            (m13, m14, m24, m34),
            (m13, m23, m24, m34),
        ]
        for c in children
            r += 1
            refined[r, :] = sort(collect(c))
        end
    end

    node_mat = reduce(vcat, reshape.(new_nodes, 1, 3))
    node_mat, refined
end

function complete_mesh(nodes::Matrix{Float64}, tets::Matrix{Int})
    facets = get_facet_list(tets)
    edges = get_edge_list(tets)
    f2e, e2f = facet_element_connectivity(tets, facets)
    Mesh3D(nodes, tets, facets, edges, f2e, e2f,
           size(nodes, 1), size(tets, 1), size(facets, 1), size(edges, 1))
end

function tet_volume(P)
    abs(det([P[2, :] .- P[1, :] P[3, :] .- P[1, :] P[4, :] .- P[1, :]])) / 6.0
end

function grad_lambda(P)
    D = ones(4, 4)
    D[:, 2:4] .= P
    invD = inv(D)
    transpose(invD[2:4, :])
end

function grad_monomial_matrix(p::Int, P)
    alpha = ijkl(p)
    beta = ijkl(p - 1)
    beta_idx = Dict{NTuple{4,Int},Int}(b => i for (i, b) in enumerate(beta))
    grads = grad_lambda(P)
    G = zeros(length(beta), length(alpha), 3)
    for (aidx, a) in enumerate(alpha)
        aa = collect(a)
        for ell in 1:4
            if aa[ell] > 0
                bb = copy(aa)
                bb[ell] -= 1
                bidx = beta_idx[Tuple(bb)]
                for d in 1:3
                    G[bidx, aidx, d] += aa[ell] * grads[ell, d]
                end
            end
        end
    end
    G
end

function cg_l2g(m::Mesh3D, p::Int)
    p in (2, 3) || error("Only p=2,3 implemented for this cross-check.")
    alpha = ijkl(p)
    deg = length(alpha)
    l2g = zeros(Int, m.NumElt, deg)
    bd = falses(10^7)

    if p == 2
        dim = m.NumNode + m.NumEdge
        bd = falses(dim)
        edge_lookup = Dict{Tuple{Int,Int},Int}()
        for i in 1:m.NumEdge
            edge_lookup[minmax(m.EdgeList[i, 1], m.EdgeList[i, 2])] = i
        end
        for e in 1:m.NumElt
            nodes = m.ElementList[e, :]
            for (d, a) in enumerate(alpha)
                nz = findall(>(0), collect(a))
                if length(nz) == 1
                    l2g[e, d] = nodes[nz[1]]
                else
                    key = minmax(nodes[nz[1]], nodes[nz[2]])
                    l2g[e, d] = m.NumNode + edge_lookup[key]
                end
            end
        end
    else
        dim = m.NumNode + 2 * m.NumEdge + m.NumF
        bd = falses(dim)
        edge_lookup = Dict{Tuple{Int,Int},Int}()
        for i in 1:m.NumEdge
            edge_lookup[minmax(m.EdgeList[i, 1], m.EdgeList[i, 2])] = i
        end
        for e in 1:m.NumElt
            nodes = m.ElementList[e, :]
            for (d, a) in enumerate(alpha)
                aa = collect(a)
                nz = findall(>(0), aa)
                if length(nz) == 1
                    l2g[e, d] = nodes[nz[1]]
                elseif length(nz) == 2
                    A, B = nz
                    key = minmax(nodes[A], nodes[B])
                    eidx = edge_lookup[key]
                    frac_num = nodes[A] < nodes[B] ? aa[B] : aa[A]
                    l2g[e, d] = m.NumNode + 2 * (eidx - 1) + frac_num
                else
                    zero_pos = findfirst(==(0), aa)
                    l2g[e, d] = m.NumNode + 2 * m.NumEdge + m.Element2Facet[e, zero_pos]
                end
            end
        end
    end

    for f in 1:m.NumF
        if m.Facet2Element[f, 2] == 0
            e = m.Facet2Element[f, 1]
            fset = Set(m.FacetList[f, :])
            loc = findfirst(i -> !(m.ElementList[e, i] in fset), 1:4)
            for (d, a) in enumerate(alpha)
                if a[loc] == 0
                    bd[l2g[e, d]] = true
                end
            end
        end
    end
    l2g, findall(bd), dim
end

function assemble_cg(m::Mesh3D, p::Int)
    alpha = ijkl(p)
    deg = length(alpha)
    l2g, bd, dim = cg_l2g(m, p)
    Mref = inner_prod_matrix_reference(p, p)
    Mref_grad = inner_prod_matrix_reference(p - 1, p - 1)
    rows = Int[]
    cols = Int[]
    kvals = Float64[]
    mvals = Float64[]
    sizehint!(rows, m.NumElt * deg^2)
    sizehint!(cols, m.NumElt * deg^2)
    sizehint!(kvals, m.NumElt * deg^2)
    sizehint!(mvals, m.NumElt * deg^2)

    for e in 1:m.NumElt
        P = m.NodeList[m.ElementList[e, :], :]
        vol = tet_volume(P)
        G = grad_monomial_matrix(p, P)
        Kloc = zeros(deg, deg)
        for d in 1:3
            Kloc .+= transpose(G[:, :, d]) * Mref_grad * G[:, :, d]
        end
        Kloc .*= vol
        Mloc = Mref .* vol
        dofs = l2g[e, :]
        for j in 1:deg, i in 1:deg
            push!(rows, dofs[i])
            push!(cols, dofs[j])
            push!(kvals, Kloc[i, j])
            push!(mvals, Mloc[i, j])
        end
    end

    K = sparse(rows, cols, kvals, dim, dim)
    M = sparse(rows, cols, mvals, dim, dim)
    interior = setdiff(collect(1:dim), bd)
    K[interior, interior], M[interior, interior], length(interior), dim
end

function cg_eigs(m::Mesh3D, p::Int; neig::Int=4)
    K, M, ndof_int, dim = assemble_cg(m, p)
    vals = eigen(Symmetric(Matrix(K)), Symmetric(Matrix(M))).values
    vals = sort(real(vals))
    vals[1:neig], ndof_int, dim
end

function assemble_cr(m::Mesh3D)
    interior_faces = findall(f -> m.Facet2Element[f, 2] != 0, 1:m.NumF)
    fmap = zeros(Int, m.NumF)
    for (i, f) in enumerate(interior_faces)
        fmap[f] = i
    end
    ndof = length(interior_faces)
    K = spzeros(Float64, ndof, ndof)
    M = spzeros(Float64, ndof, ndof)
    Mloc_ref = fill(-1.0 / 20.0, 4, 4)
    for i in 1:4
        Mloc_ref[i, i] = 2.0 / 5.0
    end
    for e in 1:m.NumElt
        P = m.NodeList[m.ElementList[e, :], :]
        vol = tet_volume(P)
        grads = grad_lambda(P)
        Kloc = 9.0 * vol .* (grads * transpose(grads))
        Mloc = vol .* Mloc_ref
        dofs = [fmap[m.Element2Facet[e, i]] for i in 1:4]
        for j in 1:4, i in 1:4
            if dofs[i] > 0 && dofs[j] > 0
                K[dofs[i], dofs[j]] += Kloc[i, j]
                M[dofs[i], dofs[j]] += Mloc[i, j]
            end
        end
    end
    K, M, ndof
end

function hmax(m::Mesh3D)
    mx = 0.0
    for e in 1:m.NumElt
        P = m.NodeList[m.ElementList[e, :], :]
        for i in 1:4, j in (i + 1):4
            mx = max(mx, norm(P[i, :] .- P[j, :]))
        end
    end
    mx
end

function cr_eigs(m::Mesh3D; neig::Int=5)
    K, M, ndof = assemble_cr(m)
    vals = eigen(Symmetric(Matrix(K)), Symmetric(Matrix(M))).values
    vals = sort(real(vals))[1:neig]
    ch = hmax(m) / sqrt(10.0)
    lbs = vals ./ (1.0 .+ ch^2 .* vals)
    vals, lbs, ndof, ch
end

function compare_vec(label, got, ref)
    println(label)
    for i in eachindex(got)
        @printf("  %d  julia %.12f  reference %.12f  diff %.3e\n",
                i, got[i], ref[i], got[i] - ref[i])
    end
end

function main()
    expected_cg = Dict(
        (2, 2) => [209.366708008859, 372.496833801567, 433.869050009704, 442.7695282258],
        (2, 3) => [198.757372321342, 350.123121926045, 355.775564719856, 401.9074231919],
        (3, 2) => [198.558427778767, 351.378479802418, 353.835740671169, 403.9673678269],
        (3, 3) => [197.417591322189, 345.586675502172, 345.738849286263, 395.0636213965],
    )
    expected_cr_nu = Dict(
        2 => [154.809042717544, 202.486779233889, 202.675696872475, 224.906408971317, 278.520983378474],
        3 => [185.289757737926, 307.383881911103, 308.146613030419, 353.103771429679, 428.250610605057],
    )
    expected_cr_lb = Dict(
        2 => [52.743809497581, 57.344072464712, 57.359213859523, 59.009952524146, 62.148883059255],
        3 => [103.158428580423, 132.447945026037, 132.589357221885, 140.274023601466, 150.785070367920],
    )

    open(joinpath(@__DIR__, "..", "results", "special_tet_crosscheck_results.txt"), "w") do io
        redirect_stdout(io) do
            println("Special tetrahedron Julia/VFEM cross-check")
            println("timestamp: ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
            println("Julia: ", VERSION)
            println()

            for p in (2, 3), ref in (2, 3)
                m = special_tet_mesh(ref)
                vals, ndof_int, dim = cg_eigs(m, p)
                @printf("CG%d ref%d: tets=%d dim=%d interior=%d\n", p, ref, m.NumElt, dim, ndof_int)
                compare_vec("  upper bounds", vals, expected_cg[(p, ref)])
                println()
            end

            for ref in (2, 3)
                m = special_tet_mesh(ref)
                vals, lbs, ndof, ch = cr_eigs(m)
                @printf("CR ref%d: tets=%d interior_faces=%d C_h=%.12f\n", ref, m.NumElt, ndof, ch)
                compare_vec("  raw CR eigenvalues", vals, expected_cr_nu[ref])
                compare_vec("  Liu lower bounds", lbs, expected_cr_lb[ref])
                println()
            end
        end
    end

    print(read(joinpath(@__DIR__, "..", "results", "special_tet_crosscheck_results.txt"), String))
end

main()
