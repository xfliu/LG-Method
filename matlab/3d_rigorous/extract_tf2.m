function extract_tf2()
% Show per-group certification (nu sign, B_spd) so spurious entries are flagged.
cfgs = [2 2; 2 3; 3 2; 3 3; 6 0; 6 1; 6 2];
for i = 1:size(cfgs,1)
    p=cfgs(i,1); r=cfgs(i,2);
    f = sprintf('special_tet_p%d_ref%d_lg_rigorous_results.mat',p,r);
    if ~exist(f,'file'), continue; end
    S=load(f); R=S.R;
    fprintf('p=%d ref=%d:\n', p, r);
    for g=1:numel(R.cert)
        cg=R.cert{g};
        nus=''; if ~isempty(cg.nu), nus=sprintf('%.3g ', inf(cg.nu)); end
        fprintf('  group spd=%d nu(inf)=[%s] err="%s"\n', cg.B_spd, nus, cg.error);
    end
    fprintf('  lb: '); fprintf('%.4f ', R.lb); fprintf('\n');
end
fprintf('EXTRACT2 DONE\n');
end
