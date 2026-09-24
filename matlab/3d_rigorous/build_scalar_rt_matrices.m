function rt_data = build_scalar_rt_matrices(mesh, M)
% BUILD_SCALAR_RT_MATRICES  Build scalar RT mixed matrices using triplet assembly.
%
% For the scalar mixed problem:
%   (sigma, tau) + (div tau, u) = 0   for all tau in RT_h
%   (div sigma, v) = (f, v)           for all v in DG_h
%
% Returns:
%   rt_data.A_rt  - (DimRT, DimRT) RT mass: (sigma_i, sigma_j)
%   rt_data.B_rt  - (DimDG, DimRT) divergence: (div sigma_i, v_j)  [note: DG rows, RT cols]
%   rt_data.M_dg  - (DimDG, DimDG) DG mass: (u_i, u_j)
%
% Uses triplet (COO) assembly for speed — avoids repeated sparse reallocation.

    global INTERVAL_MODE

    ElementList = mesh.ElementList;
    NodeList = mesh.NodeList;
    NumElt = size(ElementList, 1);
    NumF = mesh.NumF;

    [Facet2Element, Element2Facet, ElementFacetDirectSign] = ...
        mesh_get_Facet2Element_with_sign_fast(mesh);

    DegK = get_DOF(3, M);
    DegF = get_DOF(2, M);
    DegRTInner = get_DOF(3, max(M-1,0)) * 3;
    if M == 0
        DegRTInner = 0;
    end
    DegRTElt = DegF * 4 + DegRTInner;

    DOF_Facet = NumF * DegF;
    DOF_Inner = NumElt * DegRTInner;
    DimRT = DOF_Facet + DOF_Inner;
    DimDG = NumElt * DegK;

    fprintf('  Scalar RT%d: DimRT=%d, DimDG=%d\n', M, DimRT, DimDG);
    fprintf('  DegK=%d, DegF=%d, DegRTInner=%d, DegRTElt=%d\n', ...
            DegK, DegF, DegRTInner, DegRTElt);

    ijkl_list = get_IJKL(M);
    if M > 0
        ijkl_homo_list = zeros(get_DOF(2,M), 4);
        ijkl_homo_idx = zeros(get_DOF(2,M), 1);
        homo_idx = 1;
        for k = 1:size(ijkl_list, 1)
            if ijkl_list(k, 4) == 0
                ijkl_homo_list(homo_idx, :) = ijkl_list(k, :);
                ijkl_homo_idx(homo_idx) = k;
                homo_idx = homo_idx + 1;
            end
        end
    end

    % Pre-allocate triplet arrays
    nnz_A = NumElt * DegRTElt^2;
    nnz_B = NumElt * DegK * DegRTElt;
    nnz_C = NumElt * DegK^2;

    ii_A = zeros(nnz_A, 1); jj_A = zeros(nnz_A, 1); vv_A = zeros(nnz_A, 1);
    ii_B = zeros(nnz_B, 1); jj_B = zeros(nnz_B, 1); vv_B = zeros(nnz_B, 1);
    ii_C = zeros(nnz_C, 1); jj_C = zeros(nnz_C, 1); vv_C = zeros(nnz_C, 1);

    ptr_A = 0; ptr_B = 0; ptr_C = 0;

    A_M_M_Ref = getInnerProdMatrix_Reference(M, M);
    if M > 0
        A_M_M_Plus_Ref = getInnerProdMatrix_Reference(M, M+1);
        A_M_Plus_M_Plus_Ref = getInnerProdMatrix_Reference(M+1, M+1);
    end

    for IdxElt = 1:NumElt
        LocalIdxFList = Element2Facet(IdxElt, :);
        LocalNodes = NodeList(ElementList(IdxElt,:), :);
        K_vol = get_volume(LocalNodes);

        % DOF mapping
        L2GMapping = zeros(DegRTElt, 1);
        L2GSignNegtive = zeros(1, DegF*4);

        for k = 1:4
            IdxF = LocalIdxFList(k);
            LocalDOF = (1:DegF) + (k-1)*DegF;
            L2GMapping(LocalDOF) = (1:DegF) + (IdxF-1)*DegF;
            if ElementFacetDirectSign(IdxElt, k) < 0
                L2GSignNegtive(LocalDOF) = 1;
            end
        end
        L2GSignNegtivePart = find(L2GSignNegtive > 0);

        if DegRTInner > 0
            LocalDOF_inner = (1:DegRTInner) + DegF*4;
            L2GMapping(LocalDOF_inner) = (1:DegRTInner) + (IdxElt-1)*DegRTInner + DOF_Facet;
        end

        L2G_DG = (IdxElt-1)*DegK + (1:DegK);

        A_M_M = A_M_M_Ref * K_vol;

        % RT mass matrix (local)
        if M == 0
            TransMat = GetElementTransMat(LocalNodes, K_vol, M);
            LocalMatA = zeros(DegRTElt, DegRTElt);
            idx_a = 1:DegK; idx_b = DegK+(1:DegK); idx_c = 2*DegK+(1:DegK);
            LocalMatA(idx_a, idx_a) = double(A_M_M);
            LocalMatA(idx_b, idx_b) = double(A_M_M);
            LocalMatA(idx_c, idx_c) = double(A_M_M);
            LocalMat = double(TransMat)' * LocalMatA * double(TransMat);
        else
            A_M_M_Plus = A_M_M_Plus_Ref * K_vol;
            A_M_Plus_M_Plus = A_M_Plus_M_Plus_Ref * K_vol;
            MatDegreeUpByX = get_MatDegreeUpByX(M, LocalNodes);
            MatDegreeUpByX_Homo = MatDegreeUpByX(:, ijkl_homo_idx, :);

            LocalMatA = zeros(DegRTElt, DegRTElt);
            idx_a = 1:DegK; idx_b = DegK+(1:DegK); idx_c = 2*DegK+(1:DegK);
            idx_d = (3*DegK+1):DegRTElt;

            LocalMatA(idx_a, idx_a) = double(A_M_M);
            LocalMatA(idx_b, idx_b) = double(A_M_M);
            LocalMatA(idx_c, idx_c) = double(A_M_M);

            LocalMat_d = zeros(length(idx_d));
            for k = 1:3
                LocalMat_d = LocalMat_d + double(MatDegreeUpByX_Homo(:,:,k))' * double(A_M_Plus_M_Plus) * double(MatDegreeUpByX_Homo(:,:,k));
            end
            LocalMatA(idx_d, idx_d) = LocalMat_d;
            LocalMatA(idx_a, idx_d) = double(A_M_M_Plus) * double(MatDegreeUpByX_Homo(:,:,1));
            LocalMatA(idx_b, idx_d) = double(A_M_M_Plus) * double(MatDegreeUpByX_Homo(:,:,2));
            LocalMatA(idx_c, idx_d) = double(A_M_M_Plus) * double(MatDegreeUpByX_Homo(:,:,3));
            LocalMatA(idx_d, idx_a) = LocalMatA(idx_a, idx_d)';
            LocalMatA(idx_d, idx_b) = LocalMatA(idx_b, idx_d)';
            LocalMatA(idx_d, idx_c) = LocalMatA(idx_c, idx_d)';

            TransMat = GetElementTransMat(LocalNodes, K_vol, M);
            LocalMat = double(TransMat)' * LocalMatA * double(TransMat);
        end

        % Apply sign corrections
        LocalMat(L2GSignNegtivePart, :) = -1 * LocalMat(L2GSignNegtivePart, :);
        LocalMat(:, L2GSignNegtivePart) = -1 * LocalMat(:, L2GSignNegtivePart);

        % Store in triplet arrays — A_rt
        nA = DegRTElt^2;
        [jj_loc, ii_loc] = meshgrid(1:DegRTElt, 1:DegRTElt);
        ii_A(ptr_A+(1:nA)) = L2GMapping(ii_loc(:));
        jj_A(ptr_A+(1:nA)) = L2GMapping(jj_loc(:));
        vv_A(ptr_A+(1:nA)) = LocalMat(:);
        ptr_A = ptr_A + nA;

        % Divergence coupling: B_rt
        LocalMat_Div = double(TransMat)' * double(get_DivMat(M, LocalNodes))' * double(A_M_M);
        LocalMat_Div(L2GSignNegtivePart, :) = -1 * LocalMat_Div(L2GSignNegtivePart, :);
        % B_rt has DG rows, RT cols: B_rt(L2G_DG, L2GMapping) += LocalMat_Div'
        LocalB = LocalMat_Div';  % DegK x DegRTElt
        nB = DegK * DegRTElt;
        [jj_loc, ii_loc] = meshgrid(1:DegRTElt, 1:DegK);
        ii_B(ptr_B+(1:nB)) = L2G_DG(ii_loc(:));
        jj_B(ptr_B+(1:nB)) = L2GMapping(jj_loc(:));
        vv_B(ptr_B+(1:nB)) = LocalB(:);
        ptr_B = ptr_B + nB;

        % DG mass — diagonal for DG0, block-diagonal for higher
        nC = DegK^2;
        [jj_loc, ii_loc] = meshgrid(1:DegK, 1:DegK);
        ii_C(ptr_C+(1:nC)) = L2G_DG(ii_loc(:));
        jj_C(ptr_C+(1:nC)) = L2G_DG(jj_loc(:));
        vv_C(ptr_C+(1:nC)) = double(A_M_M(:));
        ptr_C = ptr_C + nC;

        if mod(IdxElt, 5000) == 0
            fprintf('    RT: %d / %d elements\n', IdxElt, NumElt);
        end
    end

    % Assemble sparse matrices from triplets
    rt_data.A_rt = sparse(ii_A(1:ptr_A), jj_A(1:ptr_A), vv_A(1:ptr_A), DimRT, DimRT);
    rt_data.B_rt = sparse(ii_B(1:ptr_B), jj_B(1:ptr_B), vv_B(1:ptr_B), DimDG, DimRT);
    rt_data.M_dg = sparse(ii_C(1:ptr_C), jj_C(1:ptr_C), vv_C(1:ptr_C), DimDG, DimDG);
    rt_data.DimRT = DimRT;
    rt_data.DimDG = DimDG;
    rt_data.DegK = DegK;

    fprintf('  Scalar RT done: A_rt %dx%d (nnz=%d), B_rt %dx%d (nnz=%d)\n', ...
            size(rt_data.A_rt,1), size(rt_data.A_rt,2), nnz(rt_data.A_rt), ...
            size(rt_data.B_rt,1), size(rt_data.B_rt,2), nnz(rt_data.B_rt));
end
