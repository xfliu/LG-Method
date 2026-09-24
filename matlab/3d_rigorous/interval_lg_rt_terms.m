function out = interval_lg_rt_terms(mesh, M, w1_d, Vdg, lhat)
%INTERVAL_LG_RT_TERMS  Rigorous Goerisch matrix A2 (full form), tight.
%
% Computes the Goerisch matrix
%   A2_{ij} = (w_i^{(1)},w_j^{(1)})_{A_rt} + lhat*(w_i^{(2)},w_j^{(2)})_{M_dg}
% entirely on the PHYSICAL element via closed-form Bernstein integrals -- no
% reference element, no loose domain-transform inverse (GetElementTransMat).
%
% Key idea (avoids the ||w1||^2 interval blow-up): the RT field w_i^{(1)} given by
% its facet-DOF coefficients w1 has a tight *raw Bernstein* representation
%   b_i = TransMat * w1_i = MatFunctionalOnBasis \ w1_i,
% obtained by a VERIFIED SOLVE whose enclosure width scales with the small raw
% norm, not with ||w1||.  Then everything is assembled from b:
%   (w_i^{(1)},w_j^{(1)})_{A_rt} = b_i' * LMA * b_j           (LMA = raw RT mass)
%   Div w_i^{(1)} = get_DivMat * b_i                          (DG coefficients)
%   w_i^{(2)} = (Div w_i^{(1)} + v_i)/lhat                    (constraint, exact)
%   lhat*(w_i^{(2)},w_j^{(2)})_{M_dg} = lhat * w2_i' * M_dg * w2_j
% LMA, get_DivMat and M_dg = getInnerProdMatrix(M,M,K_vol) are all tight (interval
% element geometry, closed-form Bernstein integrals).  This is the full Goerisch
% form (valid for an approximate w^{(1)}) and is rigorously tight.
%
% INPUT
%   mesh, M : mesh and RT/DG polynomial degree
%   w1_d    : (DimRT x n) approximate RT auxiliary fields, facet-DOF basis (double)
%   Vdg     : (DimDG x n) trial eigenfunctions in the DG (Bernstein) basis
%   lhat    : scalar shift
% OUTPUT (struct): A2_iv (interval, symmetrised), A2rt_iv, A2dg_iv, A2_max_rad.

    global INTERVAL_MODE
    old_mode = INTERVAL_MODE; INTERVAL_MODE = 1;
    cleanup = onCleanup(@() restore_mode(old_mode));

    ElementList = mesh.ElementList;
    NodeList    = mesh.NodeList;
    NumElt      = size(ElementList, 1);
    NumF        = mesh.NumF;
    n           = size(w1_d, 2);

    [~, Element2Facet, ElementFacetDirectSign] = ...
        mesh_get_Facet2Element_with_sign_fast(mesh);

    DegK       = get_DOF(3, M);
    DegF       = get_DOF(2, M);
    DegRTInner = get_DOF(3, max(M-1,0)) * 3;
    if M == 0, DegRTInner = 0; end
    DegRTElt   = DegF*4 + DegRTInner;
    DOF_Facet  = NumF * DegF;

    A_M_M_Ref      = getInnerProdMatrix_Reference(M, M);
    A_M_M_Plus_Ref = getInnerProdMatrix_Reference(M, M+1);
    A_M_Plus_Ref   = getInnerProdMatrix_Reference(M+1, M+1);
    ijkl_list = get_IJKL(M);
    hidx = find(ijkl_list(:,4) == 0);

    w1_iv  = intval(w1_d);
    Vdg_iv = intval(Vdg);
    lhat_iv = intval(lhat);

    A2rt = intval(zeros(n, n));
    A2dg = intval(zeros(n, n));

    for e = 1:NumElt
        LN_iv = intval(NodeList(ElementList(e,:), :));
        K_vol = get_volume(LN_iv);

        % ---- local RT DOF map + facet signs (as build_scalar_rt_matrices) ----
        L2GMapping = zeros(DegRTElt, 1);
        negmask = zeros(1, DegF*4);
        for k = 1:4
            IdxF = Element2Facet(e, k);
            LocalDOF = (1:DegF) + (k-1)*DegF;
            L2GMapping(LocalDOF) = (1:DegF) + (IdxF-1)*DegF;
            if ElementFacetDirectSign(e, k) < 0, negmask(LocalDOF) = 1; end
        end
        negp = find(negmask > 0);
        if DegRTInner > 0
            L2GMapping((1:DegRTInner)+DegF*4) = (1:DegRTInner) + (e-1)*DegRTInner + DOF_Facet;
        end
        L2G_DG = (e-1)*DegK + (1:DegK);

        % ---- raw Bernstein RT mass LMA (tight, no transform) ----
        A_M_M = A_M_M_Ref * K_vol;
        A_M_M_Plus = A_M_M_Plus_Ref * K_vol;
        A_M_Plus_M_Plus = A_M_Plus_Ref * K_vol;
        Up = get_MatDegreeUpByX(M, LN_iv);
        UpH = Up(:, hidx, :);
        LMA = intval(zeros(DegRTElt));
        ia = 1:DegK; ib = DegK+(1:DegK); ic = 2*DegK+(1:DegK); id = (3*DegK+1):DegRTElt;
        LMA(ia,ia)=A_M_M; LMA(ib,ib)=A_M_M; LMA(ic,ic)=A_M_M;
        LMd = intval(zeros(numel(id)));
        for k=1:3, LMd = LMd + UpH(:,:,k)'*A_M_Plus_M_Plus*UpH(:,:,k); end
        LMA(id,id)=LMd;
        LMA(ia,id)=A_M_M_Plus*UpH(:,:,1);
        LMA(ib,id)=A_M_M_Plus*UpH(:,:,2);
        LMA(ic,id)=A_M_M_Plus*UpH(:,:,3);
        LMA(id,ia)=LMA(ia,id)'; LMA(id,ib)=LMA(ib,id)'; LMA(id,ic)=LMA(ic,id)';

        MFB = get_rt_funcmat(NodeList(ElementList(e,:), :), K_vol, M);
        DivMat = get_DivMat(M, LN_iv);          % raw (a,b,c,d) -> Div DG coeffs
        M_dg_blk = getInnerProdMatrix(M, M, K_vol);

        % ---- gather signed local DOF coefficients, raw rep by verified solve ----
        w1_loc = w1_iv(L2GMapping, :);
        w1_loc(negp, :) = -w1_loc(negp, :);
        B = verified_solve(MFB, w1_loc);        % raw Bernstein coeffs (DegRTElt x n)
        v_loc = Vdg_iv(L2G_DG, :);

        % ---- A2 contributions ----
        A2rt = A2rt + B' * (LMA * B);
        DivW1 = DivMat * B;                     % (DegK x n)
        W2 = (DivW1 + v_loc) / lhat_iv;         % constraint: w2=(Div w1 + v)/lhat
        A2dg = A2dg + lhat_iv * (W2' * (M_dg_blk * W2));
    end

    A2rt = hull(A2rt, A2rt');
    A2dg = hull(A2dg, A2dg');
    A2_iv = hull(A2rt + A2dg, (A2rt + A2dg)');

    out.A2_iv = A2_iv;
    out.A2rt_iv = A2rt;
    out.A2dg_iv = A2dg;
    out.A2_max_rad = full(max(rad(A2_iv(:))));
end

function x = verified_solve(A, b)
    % Accurate ill-conditioned verified solve (Rump 'illco') keeps the enclosure
    % tight even when MatFunctionalOnBasis is very ill-conditioned (high degree);
    % fall back to the standard solver if unavailable.
    x = [];
    try
        x = verifylss(A, b, 'illco');
    catch
        x = [];
    end
    if isempty(x) || any(any(isnan(inf(x)))) || any(any(isnan(sup(x)))) || ...
            any(any(isinf(inf(x)))) || any(any(isinf(sup(x))))
        try
            x = verifylss(A, b);
        catch
            x = [];
        end
    end
    if isempty(x) || any(any(isnan(inf(x)))) || any(any(isnan(sup(x)))) || ...
            any(any(isinf(inf(x)))) || any(any(isinf(sup(x))))
        error('interval_lg_rt_terms: verified RT raw-rep solve failed.');
    end
end

function restore_mode(old_mode)
    global INTERVAL_MODE
    INTERVAL_MODE = old_mode;
end
