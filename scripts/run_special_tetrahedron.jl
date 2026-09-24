#!/usr/bin/env julia
#
# Section 5.3 of the paper: the special tetrahedron T_F with analytic
# eigenvalues.  Reproduces
#   Table 5 (CG^2/CG^3 + RT^2/RT^3 on the level-2 and level-3 meshes) and
#   Table 6 (CG^6 + RT^6 on levels 0, 1, 2),
# using the exact separators rho_2 = lambda_2 = 35 pi^2, rho_4 = lambda_4 = 40 pi^2
# and the shift lambda_hat = 200.
#
# Usage (from the repository root):
#   julia --project=. scripts/run_special_tetrahedron.jl            # interval mode (default)
#   julia --project=. scripts/run_special_tetrahedron.jl --float    # floating point only
#   julia --project=. scripts/run_special_tetrahedron.jl 6 2        # a single (p, level) pair
#
# Output: results/special_tetrahedron_<mode>_<date>.txt

using Printf, Dates
using LGPaper

const INTERVAL = !("--float" in ARGS)
args = filter(a -> !startswith(a, "--"), ARGS)
configs = isempty(args) ? [(2, 2), (2, 3), (3, 2), (3, 3), (6, 0), (6, 1), (6, 2)] :
                          [(parse(Int, args[1]), parse(Int, args[2]))]

const EXACT = FUNDAMENTAL_TET_EXACT
const RHO2 = EXACT[2]           # 35 pi^2
const RHO4 = EXACT[4]           # 40 pi^2
const LHAT = 200.0

# Certified values reported in the paper (MATLAB + INTLAB), for comparison.
const PAPER = Dict(
    (2, 2) => (lb = [171.349851167909, -881.761171388929, 65.218698493971],
               ub = [209.366708008796, 372.496833801420, 433.869050009548]),
    (2, 3) => (lb = [193.964528909108, 238.851020500837, 297.073951225893],
               ub = [198.757372321029, 350.123121925677, 355.775564719465]),
    (3, 2) => (lb = [194.723239793176, 264.794131299928, 289.115772383997],
               ub = [198.558427777527, 351.378479800498, 353.835740669255]),
    (3, 3) => (lb = [197.326517688624, 342.227090878401, 343.846088479176],
               ub = [197.417591315668, 345.586675494872, 345.738849278548]),
    (6, 0) => (lb = [195.255082156927, -762.741339957103, -762.741339957103],
               ub = [198.585756217428, 413.560898878748, 413.560898878748]),
    (6, 1) => (lb = [197.309973466981, 310.499335665668, 340.704419467296],
               ub = [197.431561337320, 345.975274300596, 349.824738185273]),
    (6, 2) => (lb = [197.391970811763, 345.382905771827, 345.423551946928],
               ub = [197.392118304215, 345.437290497899, 345.441255265817]),
)

mkpath(joinpath(@__DIR__, "..", "results"))
out = joinpath(@__DIR__, "..", "results",
               "special_tetrahedron_$(INTERVAL ? "interval" : "float")_$(Dates.format(now(), "yyyymmdd")).txt")
open(out, "w") do io
    println(io, "Special tetrahedron T_F — Lehmann–Goerisch bounds (LGPaper.jl on VFEM.jl)")
    println(io, "date: ", now(), "   Julia ", VERSION, "   mode: ", INTERVAL ? "interval" : "float")
    @printf(io, "rho_2 = lambda_2 = %.12f, rho_4 = lambda_4 = %.12f, lambda_hat = %.1f\n\n", RHO2, RHO4, LHAT)
    for (p, level) in configs
        println("\n===== T_F  p=$p  level=$level =====")
        m = fundamental_tetrahedron_mesh(level)
        r = lg_bounds(m, p, 3, [1:1 => RHO2, 2:3 => RHO4]; lhat = LHAT, interval = INTERVAL)
        for stream in (stdout, io)
            print_result(stream, r; exact = EXACT[1:3], reference = get(PAPER, (p, level), nothing),
                         label = "T_F level $level")
            println(stream)
        end
        flush(io)
    end
end
println("results written to ", out)
