function run_tf_tables()
% Compute all T_F configurations for paper section 5.3:
%   Table special-tet-convergence : CG2/CG3 on level-2 (ref2=64 tets) and level-3 (ref3=512)
%   Table special-tet-p6-exactrho : CG6 on ref0 (1), ref1 (8), ref2 (64)
cfgs = [2 2; 2 3; 3 2; 3 3; 6 0; 6 1; 6 2];
for i = 1:size(cfgs,1)
    p = cfgs(i,1); r = cfgs(i,2);
    fprintf('\n===== TF p=%d ref=%d =====\n', p, r);
    try
        run_special_tet_lg_rigorous(p, r);
    catch ME
        fprintf('TF p=%d ref=%d FAILED: %s\n', p, r, ME.message);
    end
    fprintf('===== TF p=%d ref=%d DONE =====\n', p, r);
end
fprintf('TF TABLES BATCH DONE\n');
end
