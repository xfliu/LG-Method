#!/usr/bin/env julia
#
# Section 5.4, "Efficiency of the LG lower-bound stage": reproduces Table 8.
#
# The Crouzeix–Raviart data below are the CR eigenvalues on the two finest
# regular-tetrahedron meshes (levels 5 and 6, 33496 and 247130 interior facets)
# computed with the MATLAB implementation; they are recorded here as input data.
# The certified LG interval widths are those of Table 7 (results of
# scripts/run_regular_tetrahedron.jl in interval mode agree to ~1e-9).
#
# Usage:  julia --project=. scripts/cr_efficiency_table.jl

using Printf

nu_cr5 = [150.1439321381, 275.8058811457, 275.8994274810, 275.9456931688]
nu_cr6 = [150.7548672276, 277.9517192483, 277.9557714481, 277.9628284943]
hmax_cr5 = 0.073618;  hmax_cr6 = 0.035528
ndof_cr5 = 33496;     ndof_cr6 = 247130
C = 0.3804
lb_cr5 = nu_cr5 ./ (1 .+ (C * hmax_cr5)^2 .* nu_cr5)
lb_cr6 = nu_cr6 ./ (1 .+ (C * hmax_cr6)^2 .* nu_cr6)

# Table 7 (CG^8 + RT^8 on 8 tetrahedra)
lb_lg = [150.9700841736, 278.2131313380, 278.3047808848, 278.3047808848]
ub_lg = [150.9727752307, 278.6994250425, 278.6994250425, 278.7253072428]
ndof_lg = 3960 + 1320          # dim RT^8 + dim DG^8

width = ub_lg .- lb_lg
err5 = ub_lg .- lb_cr5
err6 = ub_lg .- lb_cr6
alpha = log.(err5 ./ err6) ./ log(ndof_cr6 / ndof_cr5)   # error ~ N^(-alpha)
q = 3 .* alpha                                            # equivalent mesh-size order
N_pred = ndof_cr6 .* (err6 ./ width) .^ (1 ./ alpha)
eta = N_pred ./ ndof_lg

println("Predicted direct-CR cost for matching the CG^8+RT^8 LG interval widths (Table 8)")
println(" k   CR order q_k   eps_k^LG      N_pred_CR,k    eta_k^eff")
for k in 1:4
    @printf(" %d   %.3f          %.3e     %.2e      %.2e\n", k, q[k], width[k], N_pred[k], eta[k])
end
@printf("max eta = %.2e\n", maximum(eta))
