% Smoke test: do VFEM primitives return genuine intervals in INTERVAL_MODE?
vfem_root = '/home/xfliu/Workspace/New_Project_SchrodingerEVP/VFEM_LIB/VFEM3D';
old_pwd = pwd; cd(vfem_root); my_intlab_mode_config; cd(old_pwd);
addpath(vfem_root); addpath(genpath(fullfile(vfem_root, 'VFEM_Lib')));
global INTERVAL_MODE

LN = [0 0 0; 1 0 0; 0.5 sqrt(3)/2 0; 0.5 sqrt(3)/6 sqrt(6)/3];
p = 4;

% --- double reference ---
INTERVAL_MODE = 0;
Amm_d = getInnerProdMatrix_Reference(p,p);
vol_d = get_volume(LN);
MG_d  = get_GradMat(p, LN);
Tr_d  = GetElementTransMat(LN, vol_d, p);
Dv_d  = get_DivMat(p, LN);
Up_d  = get_MatDegreeUpByX(p, LN);

% --- interval ---
INTERVAL_MODE = 1;
Amm_i = getInnerProdMatrix_Reference(p,p);
vol_i = get_volume(intval(LN));
MG_i  = get_GradMat(p, LN);
Tr_i  = GetElementTransMat(LN, vol_i, p);
Dv_i  = get_DivMat(p, LN);
Up_i  = get_MatDegreeUpByX(p, LN);

chk = @(name,iv,dv) fprintf('%-26s isintval=%d  maxrad=%.2e  max|mid-dbl|=%.2e\n', ...
    name, isa(iv,'intval'), full(max(rad(iv(:)))), full(max(abs(mid(iv(:))-double(dv(:))))) );
chk('getInnerProdMatrix_Reference', Amm_i, Amm_d);
chk('get_volume(intval LN)',        vol_i, vol_d);
chk('get_GradMat',                  MG_i,  MG_d);
chk('GetElementTransMat',           Tr_i,  Tr_d);
chk('get_DivMat',                   Dv_i,  Dv_d);
chk('get_MatDegreeUpByX',           Up_i,  Up_d);
INTERVAL_MODE = 0;
fprintf('SMOKE OK\n');
