function R = run_special_tet_lg_rigorous(p, n_refine)
%RUN_SPECIAL_TET_LG_RIGOROUS
% Fully rigorous CG^p + RT^p Lehmann-Goerisch bounds on the fundamental
% tetrahedron T_F, whose Dirichlet-Laplacian eigenvalues are known EXACTLY:
%   mu_k = (pi^2/4) * |k|^2,  giving [80,140,140,160,208]*(pi^2/4).
%
% T_F is the falsifiable end-to-end check: every certified lower bound MUST lie
% strictly below the exact mu_k.  Exact separators are used for rho (the exact
% mu of the next eigenvalue), isolating the LG interval step from the CR stage.
%
% Rigorous replacement for run_special_tet_p12_exactrho_lg.m (which assembled
% A0,A1,A2 in double, added an ad-hoc double "A2 safety" diagonal, then wrapped
% with intval).  Here everything is interval and no safety bump is used.

    if nargin < 1, p = 12; end
    if nargin < 2, n_refine = 0; end
    log_file = sprintf('logs/special_tet_p%d_ref%d_lg_rigorous_%s.log', p, n_refine, datestr(now,'yyyymmdd'));
    diary(log_file); diary on;
    c = onCleanup(@() diary('off'));

    vfem_root = '/home/xfliu/Workspace/New_Project_SchrodingerEVP/VFEM_LIB/VFEM3D';
    old_pwd = pwd; cd(vfem_root); my_intlab_mode_config; cd(old_pwd);
    addpath(vfem_root); addpath(genpath(fullfile(vfem_root, 'VFEM_Lib')));
    addpath(fullfile(fileparts(vfem_root), 'vfem2d', 'veigs'));

    n_eig = 3; lhat = 200.0;
    lambda_exact = (pi^2/4) * [80; 140; 140; 160];
    rho2 = lambda_exact(2);   % exact separator for lambda_1
    rho4 = lambda_exact(4);   % exact separator for the lambda_2,3 cluster

    mesh = special_tet_make_mesh(n_refine);

    groups = {struct('idx',1,'rho',rho2), struct('idx',2:3,'rho',rho4)};
    R = run_conforming_lg_rigorous(mesh, p, n_eig, lhat, groups, ...
        sprintf('TF-ref%d-P%d', n_refine, p));
    R.domain = 'fundamental tetrahedron T_F';
    R.n_refine = n_refine;
    R.lambda_exact = lambda_exact(1:n_eig);
    R.rho2 = rho2; R.rho4 = rho4;
    R.ub = R.ub(:); R.lb = R.lb(:);
    R.ub_err = R.ub - lambda_exact(1:n_eig);
    R.lb_err = lambda_exact(1:n_eig) - R.lb;
    R.lb_below_exact = R.lb < lambda_exact(1:n_eig);
    R.ub_above_exact = R.ub > lambda_exact(1:n_eig);

    fprintf('\nCertified T_F P%d ref=%d rigorous LG bounds vs EXACT:\n', p, n_refine);
    fprintf('  k     lower bound        exact              upper bound        lb<exact ub>exact\n');
    for k = 1:n_eig
        fprintf('  %d   %.12f   %.12f   %.12f      %d        %d\n', ...
            k, R.lb(k), lambda_exact(k), R.ub(k), R.lb_below_exact(k), R.ub_above_exact(k));
    end
    if all(R.lb_below_exact) && all(R.ub_above_exact)
        fprintf('  PASS: all certified bounds bracket the exact eigenvalues.\n');
    else
        fprintf('  *** FAIL: a certified bound does not bracket the exact eigenvalue. ***\n');
    end
    save(sprintf('special_tet_p%d_ref%d_lg_rigorous_results.mat', p, n_refine), 'R', '-v7.3');
end
