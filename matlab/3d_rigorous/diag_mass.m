vfem_root = '/home/xfliu/Workspace/New_Project_SchrodingerEVP/VFEM_LIB/VFEM3D';
old_pwd = pwd; cd(vfem_root); my_intlab_mode_config; cd(old_pwd);
addpath(vfem_root); addpath(genpath(fullfile(vfem_root, 'VFEM_Lib')));
global INTERVAL_MODE
LN = [0,0,0; 1,0,0; 0.5,sqrt(3)/2,0; 0.5,sqrt(3)/6,sqrt(6)/3];
for p = [6 12]
    INTERVAL_MODE = 1;
    LN_iv = intval(LN); Kvol = get_volume(LN_iv);
    Aref = getInnerProdMatrix_Reference(p,p);
    M1 = Aref * Kvol;                       % what interval_lg_rt_terms used
    M2 = getInnerProdMatrix(p, p, Kvol);    % what the note mandates
    fprintf('p=%d: ref*Kvol maxrad=%.2e   getInnerProdMatrix maxrad=%.2e   max|mid diff|=%.2e\n',...
        p, full(max(rad(M1(:)))), full(max(rad(M2(:)))), full(max(abs(mid(M1(:))-mid(M2(:))))));
    % verifylss test on lhat*M
    b = intval(ones(size(M1,1),1));
    okfun = @(M) safe_solve(200*M, b);
    fprintf('   verifylss radius: ref*Kvol -> %.2e ;  getInnerProdMatrix -> %.2e\n', okfun(M1), okfun(M2));
    INTERVAL_MODE = 0;
end
fprintf('DIAG DONE\n');

function r = safe_solve(A, b)
    try
        x = verifylss(full(A), b);
        r = max(rad(x(:)));
    catch ME
        r = Inf; fprintf('   (verifylss threw: %s)\n', ME.message);
    end
end
