function rt = build_interval_rt_matrices(mesh, M)
%BUILD_INTERVAL_RT_MATRICES  Rigorous interval RT mixed matrices.
%
% Interval-arithmetic counterpart of build_scalar_rt_matrices: assembles
%   rt.A_rt : (DimRT x DimRT) RT mass   (sigma_i, sigma_j)
%   rt.B_rt : (DimDG x DimRT) divergence (div sigma_i, v_j)
%   rt.M_dg : (DimDG x DimDG) DG mass    (u_i, u_j)
% as INTLAB interval sparse matrices, with every element block computed from
% interval element geometry via the VFEM Bernstein primitives (never by
% wrapping a double matrix with intval).  Sign corrections and DOF mappings are
% identical to build_scalar_rt_matrices so the result is the rigorous enclosure
% of the same matrices.  The DG mass uses getInnerProdMatrix(M,M,K_vol), matching
% build_interval_dg_mass_matrix.

    global INTERVAL_MODE
    old_mode = INTERVAL_MODE; INTERVAL_MODE = 1;
    cleanup = onCleanup(@() restore_mode(old_mode));

    ElementList = mesh.ElementList;
    NodeList    = mesh.NodeList;
    NumElt      = size(ElementList, 1);
    NumF        = mesh.NumF;

    [~, Element2Facet, ElementFacetDirectSign] = ...
        mesh_get_Facet2Element_with_sign_fast(mesh);

    DegK       = get_DOF(3, M);
    DegF       = get_DOF(2, M);
    DegRTInner = get_DOF(3, max(M-1,0)) * 3;
    if M == 0, DegRTInner = 0; end
    DegRTElt   = DegF*4 + DegRTInner;
    DOF_Facet  = NumF * DegF;
    DimRT      = DOF_Facet + NumElt*DegRTInner;
    DimDG      = NumElt * DegK;

    A_M_M_Ref = getInnerProdMatrix_Reference(M, M);
    if M > 0
        A_M_M_Plus_Ref      = getInnerProdMatrix_Reference(M, M+1);
        A_M_Plus_M_Plus_Ref = getInnerProdMatrix_Reference(M+1, M+1);
        ijkl_list = get_IJKL(M);
        ijkl_homo_idx = zeros(get_DOF(2,M), 1); hh = 1;
        for k = 1:size(ijkl_list,1)
            if ijkl_list(k,4) == 0, ijkl_homo_idx(hh) = k; hh = hh + 1; end
        end
    end

    A_rt = intval(sparse(DimRT, DimRT));
    B_rt = intval(sparse(DimDG, DimRT));
    M_dg = intval(sparse(DimDG, DimDG));

    for e = 1:NumElt
        LocalNodes = intval(NodeList(ElementList(e,:), :));
        K_vol = get_volume(LocalNodes);

        L2GMapping = zeros(DegRTElt, 1);
        L2GSignNegtive = zeros(1, DegF*4);
        for k = 1:4
            IdxF = Element2Facet(e, k);
            LocalDOF = (1:DegF) + (k-1)*DegF;
            L2GMapping(LocalDOF) = (1:DegF) + (IdxF-1)*DegF;
            if ElementFacetDirectSign(e, k) < 0, L2GSignNegtive(LocalDOF) = 1; end
        end
        negp = find(L2GSignNegtive > 0);
        if DegRTInner > 0
            LocalDOF_inner = (1:DegRTInner) + DegF*4;
            L2GMapping(LocalDOF_inner) = (1:DegRTInner) + (e-1)*DegRTInner + DOF_Facet;
        end
        L2G_DG = (e-1)*DegK + (1:DegK);

        A_M_M = A_M_M_Ref * K_vol;

        if M == 0
            TransMat = GetElementTransMat(LocalNodes, K_vol, M);
            LMA = intval(zeros(DegRTElt));
            ia=1:DegK; ib=DegK+(1:DegK); ic=2*DegK+(1:DegK);
            LMA(ia,ia)=A_M_M; LMA(ib,ib)=A_M_M; LMA(ic,ic)=A_M_M;
            LocalMat = TransMat' * LMA * TransMat;
        else
            A_M_M_Plus      = A_M_M_Plus_Ref * K_vol;
            A_M_Plus_M_Plus = A_M_Plus_M_Plus_Ref * K_vol;
            Up = get_MatDegreeUpByX(M, LocalNodes);
            UpH = Up(:, ijkl_homo_idx, :);
            LMA = intval(zeros(DegRTElt));
            ia=1:DegK; ib=DegK+(1:DegK); ic=2*DegK+(1:DegK); id=(3*DegK+1):DegRTElt;
            LMA(ia,ia)=A_M_M; LMA(ib,ib)=A_M_M; LMA(ic,ic)=A_M_M;
            LMd = intval(zeros(numel(id)));
            for k=1:3, LMd = LMd + UpH(:,:,k)'*A_M_Plus_M_Plus*UpH(:,:,k); end
            LMA(id,id)=LMd;
            LMA(ia,id)=A_M_M_Plus*UpH(:,:,1);
            LMA(ib,id)=A_M_M_Plus*UpH(:,:,2);
            LMA(ic,id)=A_M_M_Plus*UpH(:,:,3);
            LMA(id,ia)=LMA(ia,id)'; LMA(id,ib)=LMA(ib,id)'; LMA(id,ic)=LMA(ic,id)';
            TransMat = GetElementTransMat(LocalNodes, K_vol, M);
            LocalMat = TransMat' * LMA * TransMat;
        end
        LocalMat(negp,:) = -LocalMat(negp,:);
        LocalMat(:,negp) = -LocalMat(:,negp);

        DivMat = get_DivMat(M, LocalNodes);
        LocalMat_Div = TransMat' * DivMat' * A_M_M;   % RT rows x DG cols
        LocalMat_Div(negp,:) = -LocalMat_Div(negp,:);
        LocalB = LocalMat_Div';                        % DG rows x RT cols

        % DG mass via the direct physical formula (matches build_interval_dg_mass_matrix)
        M_dg_blk = getInnerProdMatrix(M, M, K_vol);

        A_rt(L2GMapping, L2GMapping) = A_rt(L2GMapping, L2GMapping) + LocalMat;
        B_rt(L2G_DG, L2GMapping)     = B_rt(L2G_DG, L2GMapping) + LocalB;
        M_dg(L2G_DG, L2G_DG)         = M_dg_blk;   % block diagonal, no overlap
    end

    rt.A_rt = A_rt; rt.B_rt = B_rt; rt.M_dg = M_dg;
    rt.DimRT = DimRT; rt.DimDG = DimDG; rt.DegK = DegK;
end

function restore_mode(old_mode)
    global INTERVAL_MODE
    INTERVAL_MODE = old_mode;
end
