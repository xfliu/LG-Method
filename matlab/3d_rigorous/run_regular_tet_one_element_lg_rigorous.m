function R = run_regular_tet_one_element_lg_rigorous(p)
%RUN_REGULAR_TET_ONE_ELEMENT_LG_RIGOROUS
% Fully rigorous one-element bubble + RT Lehmann-Goerisch bounds for the
% Dirichlet Laplacian on the unit regular tetrahedron.
%
% Differences from the archived run_regular_tet_one_element_p12_lg.m:
%   * A0, A1 are assembled from INTERVAL element geometry (interval Kp, Mp),
%     not intval(double(Kp/Mp)).
%   * A2 = b_G(w_i,w_j) is assembled element-wise in interval arithmetic from
%     interval RT mass, divergence and DG mass blocks (interval_lg_rt_terms),
%     not w1'*intval(double(A_rt))*w1.
%   * The shifted DG component w2 is solved rigorously from the constraint.
%   * B positive-definiteness and the upper Ritz bounds are certified, not
%     only checked at the matrix midpoint.
% The approximate RT field w1 is still produced by a double saddle-point solve
% (theory allows w1 to be any fixed field; only the constraint must hold).

    if nargin < 1, p = 12; end
    log_file = sprintf('logs/regular_tet_one_element_p%d_lg_rigorous_%s.log', p, datestr(now,'yyyymmdd'));
    diary(log_file); diary on;
    c = onCleanup(@() diary('off'));

    vfem_root = '/home/xfliu/Workspace/New_Project_SchrodingerEVP/VFEM_LIB/VFEM3D';
    old_pwd = pwd; cd(vfem_root); my_intlab_mode_config; cd(old_pwd);
    addpath(vfem_root); addpath(genpath(fullfile(vfem_root, 'VFEM_Lib')));
    addpath(fullfile(fileparts(vfem_root), 'vfem2d', 'veigs'));
    global INTERVAL_MODE; INTERVAL_MODE = 0;

    qdeg = p - 4;
    n_eig = 4;
    lhat = 200.0;
    rho2 = 264.5223894731;   % CR level-6 prebound, rho_2 < lambda_2  (for lambda_1)
    rho5 = 386.4419682205;   % CR level-6 prebound, rho_5 < lambda_5  (for cluster)

    NodeList = [0,0,0; 1,0,0; 0.5,sqrt(3)/2,0; 0.5,sqrt(3)/6,sqrt(6)/3];
    ElementList = [1,2,3,4];
    mesh = make_one_tet_mesh(NodeList, ElementList);
    LN = NodeList(ElementList, :);

    DegP = get_DOF(3, p);
    DegQ = get_DOF(3, qdeg);
    ijkl_p = get_IJKL(p);
    ijkl_q = get_IJKL(qdeg);

    % Bubble embedding P_{p-4} -> P_p : multiply each monomial by l1 l2 l3 l4.
    C = zeros(DegP, DegQ);
    for j = 1:DegQ
        alpha = ijkl_q(j,:) + [1,1,1,1];
        row = find(ismember(ijkl_p, alpha, 'rows'), 1);
        C(row, j) = 1.0;
    end

    % ---- interval element stiffness/mass of the P_p space on this tetrahedron ----
    INTERVAL_MODE = 1;
    LN_iv  = intval(LN);
    vol_iv = get_volume(LN_iv);
    Amm_p_iv  = getInnerProdMatrix_Reference(p, p);
    Amm_pm_iv = getInnerProdMatrix_Reference(p-1, p-1);
    MG_iv = get_GradMat(p, LN_iv);
    Kp_iv = intval(zeros(DegP));
    for d = 1:3
        Kp_iv = Kp_iv + MG_iv(:,:,d)' * Amm_pm_iv * MG_iv(:,:,d);
    end
    Kp_iv = vol_iv * Kp_iv;
    Mp_iv = vol_iv * Amm_p_iv;
    INTERVAL_MODE = 0;

    % ---- double Ritz problem in the bubble space to pick trial functions ----
    Kp = mid(Kp_iv); Mp = mid(Mp_iv);
    Kb = C'*Kp*C; Kb = (Kb+Kb')/2;
    Mb = C'*Mp*C; Mb = (Mb+Mb')/2;
    scale = 1 ./ sqrt(diag(Mb));  S = diag(scale);
    Ks = (S*Kb*S); Ks = (Ks+Ks')/2;
    Ms = (S*Mb*S); Ms = (Ms+Ms')/2;
    [Y, D] = eig(Ks, Ms);
    [lam_ritz, si] = sort(real(diag(D)));
    Y = real(Y(:, si(1:n_eig)));
    lam_ritz = lam_ritz(1:n_eig);
    Z = S*Y;
    for k = 1:n_eig
        Z(:,k) = Z(:,k) / sqrt(Z(:,k)'*Mb*Z(:,k));
    end
    Vdg = C*Z;    % P_p DG coefficient vectors on the single element

    fprintf('One-element P%d bubble rigorous LG, unit regular tetrahedron\n', p);
    fprintf('  Dim bubble=%d, Dim P_%d=%d\n', DegQ, p, DegP);
    fprintf('  Ritz upper bounds (double): '); fprintf('%.10f ', lam_ritz); fprintf('\n');
    fprintf('  rho2=%.10f rho5=%.10f lhat=%.1f\n', rho2, rho5, lhat);

    % ---- approximate RT auxiliary field w1 (double saddle-point) ----
    rt = build_scalar_rt_matrices(mesh, p);
    A_rt = rt.A_rt; B_rt = rt.B_rt; M_dg = rt.M_dg;
    DimRT = rt.DimRT; DimDG = rt.DimDG;
    fprintf('  DimRT=%d DimDG=%d\n', DimRT, DimDG);
    D_V = lhat*M_dg; SP = [A_rt, B_rt'; B_rt, -D_V];
    [Lsp,Usp,Psp,Qsp] = lu(SP);
    w1_d = zeros(DimRT, n_eig);
    w2_d = zeros(DimDG, n_eig);
    for i = 1:n_eig
        rhs = [zeros(DimRT,1); -M_dg*Vdg(:,i)];
        sol = Qsp*(Usp\(Lsp\(Psp*rhs)));
        sol = sol + Qsp*(Usp\(Lsp\(Psp*(rhs - SP*sol))));
        w1_d(:,i) = sol(1:DimRT);
        w2_d(:,i) = sol(DimRT+1:end);
    end
    fprintf('  ||w1|| per col: '); fprintf('%.2e ', sqrt(sum(w1_d.^2,1))); fprintf('\n');

    % ---- rigorous A0, A1 (interval) ----
    INTERVAL_MODE = 1;
    Vdg_iv = intval(Vdg); lhat_iv = intval(lhat);
    A0_iv = Vdg_iv' * ((Kp_iv + lhat_iv*Mp_iv) * Vdg_iv);
    A1_iv = Vdg_iv' * (Mp_iv * Vdg_iv);
    A0_iv = hull(A0_iv, A0_iv'); A1_iv = hull(A1_iv, A1_iv');

    % ---- rigorous A2 (full form, physical-element Bernstein assembly) ----
    fprintf('  Assembling rigorous Goerisch A2 (full form, tight via b=MFB\\w1)...\n');
    rtterms = interval_lg_rt_terms(mesh, p, w1_d, Vdg, lhat);
    A2_iv = rtterms.A2_iv;
    fprintf('  A2(k,k)*(lambda_k+lhat) [should be ~1]: ');
    for k = 1:n_eig, fprintf('%.6f ', mid(A2_iv(k,k))*(lam_ritz(k)+lhat)); end
    fprintf('\n  A2 max interval radius: %.3e\n', max(rad(A2_iv(:))));

    % ---- certified upper bounds (Rayleigh-Ritz values enclosed) ----
    [ub_enc, ~] = veig(A0_iv - lhat_iv*A1_iv, A1_iv, 1:n_eig);
    ub = sup(ub_enc(1:n_eig));

    % ---- certified lower bounds ----
    r1  = lg_certify(A0_iv, A1_iv, A2_iv, 1,   rho2, lhat);
    rcl = lg_certify(A0_iv, A1_iv, A2_iv, 2:4, rho5, lhat);
    INTERVAL_MODE = 0;
    lb = [r1.lb; rcl.lb];

    fprintf('\n  B PD: lambda1 group spd=%d (margin %.3e); cluster spd=%d (margin %.3e)\n', ...
        r1.B_spd, r1.B_pd_margin, rcl.B_spd, rcl.B_pd_margin);
    fprintf('\nCertified one-element P%d rigorous LG bounds:\n', p);
    fprintf('  k     lower bound        upper bound        rel gap\n');
    for k = 1:n_eig
        fprintf('  %d   %.12f   %.12f   %.6f%%\n', k, lb(k), ub(k), (ub(k)-lb(k))/ub(k)*100);
    end

    R = struct('domain','unit regular tetrahedron', ...
        'method','rigorous one-element P bubble + RT Lehmann-Goerisch', ...
        'p',p,'qdeg',qdeg,'n_eig',n_eig,'lhat',lhat,'rho2',rho2,'rho5',rho5, ...
        'lb',lb,'ub',ub,'ub_ritz_double',lam_ritz, ...
        'rel_gap',(ub-lb)./ub*100, ...
        'A0_iv',A0_iv,'A1_iv',A1_iv,'A2_iv',A2_iv, ...
        'A2_max_rad',max(rad(A2_iv(:))), ...
        'B_spd_lam1',r1.B_spd,'B_spd_cluster',rcl.B_spd);
    save(sprintf('regular_tet_one_element_p%d_lg_rigorous_results.mat', p), 'R', '-v7.3');
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
