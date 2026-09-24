function R = run_conforming_lg_rigorous(mesh, p, n_eig, lhat, groups, label)
%RUN_CONFORMING_LG_RIGOROUS  Rigorous CG^p + RT^p Lehmann-Goerisch engine.
%
% Shared core for conforming-mesh drivers (regular tet refinement, T_F).
%   * A0, A1 from interval CG assembly (assemble_cg_lg_small_matrices_iv)
%   * A2, w2 from interval RT/DG assembly (interval_lg_rt_terms)
%   * certified upper (Ritz) and lower (LG) bounds (lg_certify)
% w1 is an approximate double saddle-point field (theory allows any fixed w1).
%
% groups : cell array, each struct('idx', <subset of 1:n_eig>, 'rho', <double>)
%          giving the a priori separation lower bound rho for that group.
% Returns struct R with lb, ub, interval matrices and diagnostics.

    global INTERVAL_MODE; INTERVAL_MODE = 0;
    NodeList = mesh.NodeList; ElementList = mesh.ElementList;
    NumElt = mesh.NumElt;

    DegK = get_DOF(3, p);
    [L2G, DOF_BD, DimCG] = local_register_cg_dof(mesh, p);
    interior = setdiff(1:DimCG, DOF_BD);
    DimInt = numel(interior);
    fprintf('[%s] CG^%d: DimCG=%d interior=%d (%d tets)\n', label, p, DimCG, DimInt, NumElt);

    % ---- double CG stiffness/mass for the Ritz problem ----
    Amm = double(getInnerProdMatrix_Reference(p, p));
    Amm1 = double(getInnerProdMatrix_Reference(p-1, p-1));
    iK = zeros(NumElt*DegK^2,1); jK = iK; vK = iK; vM = iK;
    ptr = 0;
    for e = 1:NumElt
        LN = NodeList(ElementList(e,:), :);
        vol = get_volume(LN);
        dofs = L2G(e,:);
        MG = get_GradMat(p, LN);
        LK = zeros(DegK);
        for d = 1:3, LK = LK + double(MG(:,:,d))'*Amm1*double(MG(:,:,d)); end
        LK = LK*vol; LM = Amm*vol;
        [jl, il] = meshgrid(1:DegK, 1:DegK);
        rr = ptr + (1:DegK^2);
        iK(rr)=dofs(il(:)); jK(rr)=dofs(jl(:)); vK(rr)=LK(:); vM(rr)=LM(:);
        ptr = ptr + DegK^2;
    end
    K_g = sparse(iK,jK,vK,DimCG,DimCG);
    M_g = sparse(iK,jK,vM,DimCG,DimCG);
    K_int = (K_g(interior,interior)+K_g(interior,interior)')/2;
    M_int = (M_g(interior,interior)+M_g(interior,interior)')/2;

    if DimInt <= 2200
        [Ve, De] = eig(full(K_int), full(M_int));
    else
        opts.disp=0; [Ve, De] = eigs(K_int, M_int, n_eig+2, 1.0, opts);
    end
    [lam_ritz_all, si] = sort(real(diag(De)));
    Ve = real(Ve(:, si(1:n_eig)));
    lam_ritz = lam_ritz_all(1:n_eig);
    fprintf('[%s] Ritz upper bounds (double): ', label); fprintf('%.10f ', lam_ritz); fprintf('\n');

    v_all = zeros(DimCG, n_eig);
    for k = 1:n_eig
        vf = zeros(DimCG,1); vf(interior) = Ve(:,k);
        vf = vf / sqrt(vf(interior)'*M_int*vf(interior));
        v_all(:,k) = vf;
    end

    % ---- approximate RT auxiliary field w1 (double saddle-point) ----
    rt = build_scalar_rt_matrices(mesh, p);
    A_rt = rt.A_rt; B_rt = rt.B_rt; M_dg = rt.M_dg;
    DimRT = rt.DimRT; DimDG = rt.DimDG; DegK_dg = rt.DegK;
    D_V = lhat*M_dg; SP = [A_rt, B_rt'; B_rt, -D_V];
    fprintf('[%s] Factoring shifted SP (%d)...\n', label, DimRT+DimDG);
    [Lsp,Usp,Psp,Qsp] = lu(SP);
    w1_d = zeros(DimRT, n_eig); w2_d = zeros(DimDG, n_eig); Vdg = zeros(DimDG, n_eig);
    for i = 1:n_eig
        for e = 1:NumElt
            Vdg((e-1)*DegK_dg+(1:DegK_dg), i) = v_all(L2G(e,:), i);
        end
        rhs = [zeros(DimRT,1); -M_dg*Vdg(:,i)];
        sol = Qsp*(Usp\(Lsp\(Psp*rhs)));
        sol = sol + Qsp*(Usp\(Lsp\(Psp*(rhs - SP*sol))));
        w1_d(:,i) = sol(1:DimRT);
        w2_d(:,i) = sol(DimRT+1:end);
    end

    % ---- rigorous interval Goerisch matrices ----
    fprintf('[%s] Assembling rigorous A0,A1 (interval)...\n', label);
    [A0_iv, A1_iv] = assemble_cg_lg_small_matrices_iv(mesh, p, v_all, L2G, lhat);
    fprintf('[%s] Assembling rigorous A2 (full form, tight via b=MFB\\w1)...\n', label);
    rtt = interval_lg_rt_terms(mesh, p, w1_d, Vdg, lhat);
    A2_iv = rtt.A2_iv;
    fprintf('[%s] A2(k,k)*(lambda_k+lhat) [~1]: ', label);
    for k = 1:n_eig, fprintf('%.6f ', mid(A2_iv(k,k))*(lam_ritz(k)+lhat)); end
    fprintf('\n[%s] A2 max interval radius: %.3e\n', label, max(rad(A2_iv(:))));

    % ---- certified two-sided bounds ----
    INTERVAL_MODE = 1; lhat_iv = intval(lhat);
    [ub_enc, ~] = veig(A0_iv - lhat_iv*A1_iv, A1_iv, 1:n_eig);
    ub = sup(ub_enc(1:n_eig)); ub = ub(:);
    lam_ritz = lam_ritz(:);
    lb = nan(n_eig,1);
    cert = cell(numel(groups),1);
    for g = 1:numel(groups)
        cg = lg_certify(A0_iv, A1_iv, A2_iv, groups{g}.idx, groups{g}.rho, lhat);
        lb(groups{g}.idx) = cg.lb;
        cert{g} = cg;
        fprintf('[%s] group [%s] rho=%.6f: B_spd=%d margin=%.3e %s\n', label, ...
            num2str(groups{g}.idx), groups{g}.rho, cg.B_spd, cg.B_pd_margin, cg.error);
    end
    INTERVAL_MODE = 0;

    R = struct('label',label,'p',p,'n_eig',n_eig,'lhat',lhat, ...
        'n_tets',NumElt,'dim_cg',DimCG,'dim_int',DimInt,'dim_rt',DimRT,'dim_dg',DimDG, ...
        'ub',ub,'lb',lb,'ub_ritz_double',lam_ritz, ...
        'A0_iv',A0_iv,'A1_iv',A1_iv,'A2_iv',A2_iv, ...
        'A2_max_rad',max(rad(A2_iv(:))), ...
        'cert',{cert});
end
