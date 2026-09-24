S=load('regular_tet_ref1_p8_lg_rigorous_results.mat'); R=S.R;
fprintf('label=%s p=%d ntets=%d\n', R.label, R.p, R.n_tets);
fprintf('dim_cg=%d dim_int=%d dim_rt=%d dim_dg=%d\n', R.dim_cg, R.dim_int, R.dim_rt, R.dim_dg);
fprintf('A2_max_rad=%.3e\n', R.A2_max_rad);
for k=1:4
  fprintf('k=%d lb=%.10f ub=%.10f gap=%.4f%% width=%.4e\n', k, R.lb(k), R.ub(k), (R.ub(k)-R.lb(k))/R.ub(k)*100, R.ub(k)-R.lb(k));
end
