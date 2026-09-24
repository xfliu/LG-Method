my_intlab_mode_config

global INTERVAL_MODE

INTERVAL_MODE = 0;
%clear;clc;INTERVAL_MODE

% % Mesh 1:
% mesh_path   = './mesh_data/UnitSquare8x8/';
% mesh = read_mesh_from_folder(mesh_path);
% neig = 4;
% rho = 100;

% % Mesh 2:
%mesh = read_dolfin_mesh("./mesh_data/dumbbell_graded_R2.xml");
%rho = 57;
%neig = 6;

%rho = 44.0;
%neig = 2;

% % Mesh 0
%mesh = read_dolfin_mesh("./mesh_data/dumbbell_graded_mesh_R2.xml");
rho = 74.0;
neig = 6;

% % Mesh 3 w=w=0.390625;
%mesh = read_dolfin_mesh("./mesh_data/dumbbell_graded_mesh3_R2.xml");
rho = 66.0;
neig = 6;

% % Mesh 3 w=w=0.390625;
% mesh = read_dolfin_mesh("./mesh_data/dumbbell_graded_0.40625.xml");
rho = 49.2;
neig = 4;

% % Mesh 3 w=w=0.390625;
mesh = read_dolfin_mesh("./mesh_data/dumbbell_graded_0.40625.xml");
rho = 74.0;
neig = 6;
disp(size(mesh.elements));

%Note : lower bound for lambda 2: R1_2, deg=2, eig_2 > 19.16932
%Note : upper bound for lambda 1: R2, deg=3, eig_1 < 19.1689

%write_mesh_xml(vert, tri, [mesh_path, 'usd.xml'])

% vert = I_intval(mesh.nodes);
% tri  = mesh.elements;
% 
% trimesh(tri, vert(:,1), vert(:,2), 'color', 'k')
% axis equal
% axis off

% mesh_size = [size(vert, 1), size(tri, 1)]

%================================================================
% disp('Compute Laplaian eig using CR element');
% [CR_eig, Ch] = CR_steklov_eig(vert, edge, tri, bd, neig);
% Ch = Ch(1);
% disp('compute steklov eig lower bound by our thm')
% CR_eig_low = CR_eig ./ (1 + CR_eig .* Ch^2);  
%================================================================
disp('Compute Laplace Eig Using Lagrange Element (Continuous Galerkin =CG))');
Lagrange_order = 1;
[LA_eig, LA_eigf, LA_A, LA_M] = Lagrange_laplace_eig(Lagrange_order, mesh, neig);

disp("Size of matrix:")
disp(size(LA_A,1))

if INTERVAL_MODE 
    A2 = LA_eigf'*LA_A*LA_eigf;
    M2 = LA_eigf'*LA_M*LA_eigf;
    LA_eig = veig(hull(A2, A2'), hull(M2, M2'));
end

format long 
LA_eig

LA_eigf = LA_eigf(:,3:4);

%% Apply the Lehmann-Goerisch method to obtain sharp lower bounds.
%================================================================
disp('Compute Laplace Eig Lower Bound by the Lehmann-Goerisch method')
A0 = LA_eigf' * LA_A * LA_eigf;
A1 = LA_eigf' * LA_M * LA_eigf;

RT_order = Lagrange_order;
mat_b_w_w = RT_Hdiv_problem(mesh, RT_order, LA_eigf);
A_lg = mat_b_w_w;

%A2_diff = A_lg - A2

AL = A0 - rho * A1;
BL = A0 - 2*rho*A1 + rho*rho*A_lg;

if INTERVAL_MODE 
    sym_AL = hull(AL, AL');
    sym_BL = hull(BL, BL');
    LG_eig_low = veig(sym_AL, sym_BL);
    [var_,idx] = sort(I_inf(LG_eig_low));
    LG_eig_low = LG_eig_low(idx);
else
    LG_eig_low = sort(eig(AL, BL));
    LG_eig_low
end
LG_eig_low = rho - rho./(1-LG_eig_low(end:-1:1));

format long
LA_eig
LG_eig_low

eig_lower_and_upper = [LG_eig_low, LA_eig(3:4)]

save result_0.mat  eig_lower_and_upper 
save LG_bound_0.mat  LG_eig_low

