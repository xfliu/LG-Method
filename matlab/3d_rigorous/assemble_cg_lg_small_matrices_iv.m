function [A0_iv, A1_iv] = assemble_cg_lg_small_matrices_iv(mesh, p, v_all, L2G, lhat)
%ASSEMBLE_CG_LG_SMALL_MATRICES_IV  Rigorous A0, A1 for conforming CG^p trials.
%
% Computes, in INTLAB interval arithmetic and element by element,
%   A0(i,j) = a(v_i,v_j) + lhat*(v_i,v_j) = M(v_i,v_j)   [shifted energy]
%   A1(i,j) = (v_i,v_j)                  = N(v_i,v_j)    [L2 mass]
% with EVERY element stiffness/mass block assembled from interval geometry via
% the VFEM Bernstein primitives -- not by wrapping a floating-point K_int/M_int
% with intval().  Trial coefficient vectors v_all are FIXED floats (admissible).
%
% v_all : (DimCG x n) global CG coefficient vectors (zero on Dirichlet DOFs).
% L2G   : (NumElt x DegK) local-to-global CG DOF map.

    global INTERVAL_MODE
    old_mode = INTERVAL_MODE;
    INTERVAL_MODE = 1;
    cleanup = onCleanup(@() restore_mode(old_mode));

    ElementList = mesh.ElementList;
    NodeList    = mesh.NodeList;
    NumElt      = size(ElementList, 1);
    DegK        = get_DOF(3, p);
    n           = size(v_all, 2);

    Amm_iv  = getInnerProdMatrix_Reference(p, p);
    Amm1_iv = getInnerProdMatrix_Reference(p-1, p-1);
    v_iv    = intval(v_all);
    lhat_iv = intval(lhat);

    A0_iv = intval(zeros(n, n));
    A1_iv = intval(zeros(n, n));

    for e = 1:NumElt
        LN  = intval(NodeList(ElementList(e,:), :));
        vol = get_volume(LN);
        MG  = get_GradMat(p, LN);             % (DegK x DegK x 3), interval
        LK  = intval(zeros(DegK));
        for d = 1:3
            LK = LK + MG(:,:,d)' * Amm1_iv * MG(:,:,d);
        end
        LK = LK * vol;                        % local stiffness  a(.,.)
        LM = Amm_iv * vol;                    % local mass       (.,.)

        v_loc = v_iv(L2G(e,:), :);            % (DegK x n)
        A0_iv = A0_iv + v_loc' * ((LK + lhat_iv*LM) * v_loc);
        A1_iv = A1_iv + v_loc' * (LM * v_loc);
    end

    A0_iv = hull(A0_iv, A0_iv');
    A1_iv = hull(A1_iv, A1_iv');
end

function restore_mode(old_mode)
    global INTERVAL_MODE
    INTERVAL_MODE = old_mode;
end
