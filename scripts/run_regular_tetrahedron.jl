#!/usr/bin/env julia
#
# Section 5.4 of the paper: the unit regular tetrahedron.  Reproduces Table 7
# (CG^8 + RT^8 on the 8-tetrahedron mesh, lambda_hat = 200).
#
# Stage 1.  The separators used in the paper,
#     rho_2 = 264.5223894731 < lambda_2,   rho_5 = 386.4419682205 < lambda_5,
# come from the Crouzeix–Raviart projection bound on a level-6 mesh with 127661
# tetrahedra (MATLAB computation, see Section 5.4).  This script recomputes the
# CR projection bounds on the uniform red refinements up to `--cr-level L`
# (default 3; level 4 has 4096 tetrahedra and takes a few minutes with the dense
# eigensolver) to show the first stage in Julia, and uses the paper's values for
# the second stage unless `--rho-from-cr` is given.
#
# Usage (from the repository root):
#   julia --project=. scripts/run_regular_tetrahedron.jl [--float] [--cr-level 4] [--rho-from-cr] [p]

using Printf, Dates
using LGPaper

const INTERVAL = !("--float" in ARGS)
function optval(flag, default)
    i = findfirst(==(flag), ARGS)
    i === nothing ? default : parse(Int, ARGS[i + 1])
end
cr_level = optval("--cr-level", 3)
rho_from_cr = "--rho-from-cr" in ARGS
pos = filter(a -> !startswith(a, "--") && all(isdigit, a), ARGS)
p = isempty(pos) ? 8 : parse(Int, pos[end])
p == optval("--cr-level", -1) && (p = 8)   # guard against reading the option value as p

const LHAT = 200.0
const RHO2_PAPER = 264.5223894731
const RHO5_PAPER = 386.4419682205
const PAPER_TABLE7 = (lb = [150.970084173612, 278.213131338008, 278.304780884805, 278.304780884805],
                      ub = [150.972775230721, 278.699425042544, 278.699425042544, 278.725307242761])

mkpath(joinpath(@__DIR__, "..", "results"))
out = joinpath(@__DIR__, "..", "results",
               "regular_tetrahedron_p$(p)_$(INTERVAL ? "interval" : "float")_$(Dates.format(now(), "yyyymmdd")).txt")
open(out, "w") do io
    println(io, "Unit regular tetrahedron — two-stage Lehmann–Goerisch bounds (LGPaper.jl on VFEM.jl)")
    println(io, "date: ", now(), "   Julia ", VERSION, "   mode: ", INTERVAL ? "interval" : "float")
    println(io, "\nStage 1: Crouzeix–Raviart projection bounds, C_h = 0.3804 h_max")
    println(io, " level   tets   interior facets   h_max      lambda_h^CR(1..5)                               lower bounds(1..5)")
    rho = (RHO2_PAPER, RHO5_PAPER)
    for level in 1:cr_level
        m = regular_tetrahedron_mesh(level)
        c = cr_projection_bounds(m, 5)
        line = @sprintf("  %d   %6d   %8d   %.6f   %s   %s", level, m.NumElt, c.ndof, c.hmax,
                        join([@sprintf("%.4f", x) for x in c.cr_eig], " "),
                        join([@sprintf("%.4f", x) for x in c.lower], " "))
        println(line); println(io, line)
        rho = (c.lower[2], c.lower[5])
    end
    if !rho_from_cr
        rho = (RHO2_PAPER, RHO5_PAPER)
        println(io, "\nSeparators used below (paper, CR level-6 mesh with 127661 tets): rho_2 = $(rho[1]), rho_5 = $(rho[2])")
    else
        println(io, "\nSeparators used below (CR level $cr_level above): rho_2 = $(rho[1]), rho_5 = $(rho[2])")
    end
    println("\n===== regular tetrahedron, 8 tets, CG^$p + RT^$p =====")
    m = regular_tetrahedron_mesh(1)
    r = lg_bounds(m, p, 4, [1:1 => rho[1], 2:4 => rho[2]]; lhat = LHAT, interval = INTERVAL)
    ref = (p == 8 && !rho_from_cr) ? PAPER_TABLE7 : nothing
    for stream in (stdout, io)
        println(stream)
        print_result(stream, r; reference = ref, label = "regular tetrahedron")
        println(stream, "  relative gaps (ub-lb)/ub: ",
                join([@sprintf("%.6f%%", 100 * (r.ub[k] - r.lb[k]) / r.ub[k]) for k in 1:4], "  "))
    end
end
println("results written to ", out)
