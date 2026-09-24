% Postprocess the regular-tetrahedron CR-vs-LG efficiency estimate.
%
% The direct CR prediction is fitted from the observed CR lower-bound errors
% on the level-5 and level-6 regular-tetrahedron meshes.  The target accuracy
% is the certified width of the one-refinement CG12+RT12 LG interval.

code_dir = fileparts(mfilename('fullpath'));
ref = load(fullfile(code_dir, 'regular_tet_ref1_p12_lg_results.mat'));

R = ref.R;
lambda_ref = R.ub(:);
target_width = R.ub(:) - R.lb(:);

nu_cr5 = [
    150.1439321381;
    275.8058811457;
    275.8994274810;
    275.9456931688
];
nu_cr6 = [
    150.7548672276;
    277.9517192483;
    277.9557714481;
    277.9628284943
];
hmax_cr5 = 0.073618;
hmax_cr6 = 0.035528;
Ch_cr5 = 0.3804 * hmax_cr5;
Ch_cr6 = 0.3804 * hmax_cr6;
lb_cr5 = nu_cr5 ./ (1 + Ch_cr5^2 .* nu_cr5);
lb_cr6 = nu_cr6 ./ (1 + Ch_cr6^2 .* nu_cr6);

ndof_cr5 = 33496;
ndof_cr6 = 247130;
ndof_lg = R.dim_rt + R.dim_dg;
ndof_lg_total = ndof_lg + R.dim_int;

err_cr5 = lambda_ref - lb_cr5;
err_cr6 = lambda_ref - lb_cr6;

alpha_dof = log(err_cr5 ./ err_cr6) ./ log(ndof_cr6 / ndof_cr5);
order_h_equiv = 3 * alpha_dof;
pred_cr_dof = ndof_cr6 * (err_cr6 ./ target_width) .^ (1 ./ alpha_dof);
eff_lg = pred_cr_dof / ndof_lg;
eff_lg_total = pred_cr_dof / ndof_lg_total;

fprintf('Regular tetrahedron CR-vs-LG efficiency estimate\n');
fprintf('CR level 5 DOFs = %d, CR level 6 DOFs = %d\n', ndof_cr5, ndof_cr6);
fprintf('LG lower-bound mixed dimension = %d; including CG interior space = %d\n', ...
    ndof_lg, ndof_lg_total);
fprintf('%2s %12s %14s %14s %14s %14s\n', ...
    'k', 'h-order', 'target width', 'CR pred DOF', 'E_LG', 'E_total');
for k = 1:numel(lambda_ref)
    fprintf('%2d %12.6f %14.6e %14.6e %14.6e %14.6e\n', ...
        k, order_h_equiv(k), target_width(k), pred_cr_dof(k), ...
        eff_lg(k), eff_lg_total(k));
end

save(fullfile(code_dir, 'regular_tet_cr_efficiency_results.mat'), ...
    'lambda_ref', 'target_width', 'lb_cr5', 'lb_cr6', ...
    'ndof_cr5', 'ndof_cr6', 'ndof_lg', 'ndof_lg_total', ...
    'err_cr5', 'err_cr6', 'alpha_dof', 'order_h_equiv', ...
    'pred_cr_dof', 'eff_lg', 'eff_lg_total');
