function extract_tf()
ex = (pi^2/4)*[80;140;140;160];
load_lbub = @(p,r) getfield(load(sprintf('special_tet_p%d_ref%d_lg_rigorous_results.mat',p,r)),'R');
ord = @(ec,ef) log(ec./ef)/log(2);
fprintf('=== Table 3 (CG2/CG3, level2=ref2, level3=ref3) ===\n');
for p = [2 3]
    R2 = load_lbub(p,2); R3 = load_lbub(p,3);
    for r = [2 3]
        R = load_lbub(p,r);
        fprintf('p=%d lvl=%d: ', p, r);
        for k=1:3, fprintf('[%.2f, %.2f]  ', R.lb(k), R.ub(k)); end
        fprintf('\n');
    end
    ubo = ord(R2.ub(1:3)-ex(1:3), R3.ub(1:3)-ex(1:3));
    lbo = ord(ex(1:3)-R2.lb(1:3), ex(1:3)-R3.lb(1:3));
    fprintf('  UB order (l2->l3): %.1f %.1f %.1f | LB order: %.1f %.1f %.1f\n', ubo, lbo);
end
fprintf('=== Table 4 (CG6, ref0/1/2) ===\n');
R0=load_lbub(6,0); R1=load_lbub(6,1); R2=load_lbub(6,2);
for r=[0 1 2]
    R=load_lbub(6,r); fprintf('ref=%d (#elt=%d): ', r, R.n_tets);
    for k=1:3, fprintf('[%.3f, %.3f]  ', R.lb(k), R.ub(k)); end
    fprintf('\n');
end
ubo01=ord(R0.ub(1:3)-ex(1:3),R1.ub(1:3)-ex(1:3)); lbo01=ord(ex(1:3)-R0.lb(1:3),ex(1:3)-R1.lb(1:3));
ubo12=ord(R1.ub(1:3)-ex(1:3),R2.ub(1:3)-ex(1:3)); lbo12=ord(ex(1:3)-R1.lb(1:3),ex(1:3)-R2.lb(1:3));
fprintf('  order ref0->1: UB %.1f %.1f %.1f LB %.1f %.1f %.1f\n', ubo01,lbo01);
fprintf('  order ref1->2: UB %.1f %.1f %.1f LB %.1f %.1f %.1f\n', ubo12,lbo12);
fprintf('EXTRACT DONE\n');
end
