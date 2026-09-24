function R = run_regular_tet_ref1_lg_rigorous(p)
%RUN_REGULAR_TET_REF1_LG_RIGOROUS
% Fully rigorous CG^p + RT^p Lehmann-Goerisch bounds for lambda_1..4 of the
% Dirichlet Laplacian on the unit regular tetrahedron, one uniform refinement
% (8 tetrahedra).  Rigorous replacement for run_regular_tet_ref1_p12_lg.m:
% A0,A1,A2 are assembled in interval arithmetic, not intval(double).

    if nargin < 1, p = 12; end
    log_file = sprintf('logs/regular_tet_ref1_p%d_lg_rigorous_%s.log', p, datestr(now,'yyyymmdd'));
    diary(log_file); diary on;
    c = onCleanup(@() diary('off'));

    vfem_root = '/home/xfliu/Workspace/New_Project_SchrodingerEVP/VFEM_LIB/VFEM3D';
    old_pwd = pwd; cd(vfem_root); my_intlab_mode_config; cd(old_pwd);
    addpath(vfem_root); addpath(genpath(fullfile(vfem_root, 'VFEM_Lib')));
    addpath(fullfile(fileparts(vfem_root), 'vfem2d', 'veigs'));

    n_eig = 4; lhat = 200.0;
    rho2 = 264.5223894731;   % CR level-6 prebound, rho_2 < lambda_2
    rho5 = 386.4419682205;   % CR level-6 prebound, rho_5 < lambda_5

    NodeList = [0,0,0; 1,0,0; 0.5,sqrt(3)/2,0; 0.5,sqrt(3)/6,sqrt(6)/3];
    mesh0 = make_one_tet_mesh(NodeList, [1,2,3,4]);
    mesh = refine_tet_mesh(mesh0);   % 1 -> 8 tetrahedra

    groups = {struct('idx',1,'rho',rho2), struct('idx',2:4,'rho',rho5)};
    R = run_conforming_lg_rigorous(mesh, p, n_eig, lhat, groups, ...
        sprintf('reg-8elt-P%d', p));
    R.domain = 'unit regular tetrahedron, one uniform refinement (8 tets)';
    R.rho2 = rho2; R.rho5 = rho5; R.rel_gap = (R.ub - R.lb)./R.ub*100;

    fprintf('\nCertified 8-element P%d rigorous LG bounds:\n', p);
    fprintf('  k     lower bound        upper bound        rel gap\n');
    for k = 1:n_eig
        fprintf('  %d   %.12f   %.12f   %.6f%%\n', k, R.lb(k), R.ub(k), R.rel_gap(k));
    end
    save(sprintf('regular_tet_ref1_p%d_lg_rigorous_results.mat', p), 'R', '-v7.3');
end

function mesh = make_one_tet_mesh(NodeList, ElementList)
    ElementList = sort(ElementList, 2);
    FacetList = get_FacetList(ElementList);
    mesh.ElementList = ElementList; mesh.NodeList = NodeList;
    mesh.FacetList = FacetList; mesh.NumF = size(FacetList,1);
    mesh.NumNode = size(NodeList,1); mesh.NumElt = size(ElementList,1);
    edges = [ElementList(:,[1 2]); ElementList(:,[1 3]); ElementList(:,[1 4]); ...
             ElementList(:,[2 3]); ElementList(:,[2 4]); ElementList(:,[3 4])];
    mesh.EdgeList = unique(sort(edges,2),'rows'); mesh.NumEdge = size(mesh.EdgeList,1);
    [mesh.Facet2Element, mesh.Element2Facet] = local_f2e(ElementList, FacetList);
end

function [F2E, E2F] = local_f2e(ElementList, FacetList)
    NumElt = size(ElementList,1); NumF = size(FacetList,1);
    lf = [2 3 4; 1 3 4; 1 2 4; 1 2 3];
    E2F = zeros(NumElt,4); F2E = zeros(NumF,2);
    for e = 1:NumElt
        for k = 1:4
            f = sort(ElementList(e, lf(k,:)),2);
            idx = find(ismember(FacetList, f, 'rows'),1);
            E2F(e,k) = idx;
            if F2E(idx,1)==0, F2E(idx,1)=e; else, F2E(idx,2)=e; end
        end
    end
end
