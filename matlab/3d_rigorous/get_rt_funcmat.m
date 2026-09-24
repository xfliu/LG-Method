function MFB = get_rt_funcmat(LocalNodes, Vol, M)
%GET_RT_FUNCMAT  The (row-scaled) RT DOF-functional-on-Bernstein-basis matrix.
%
% Identical to VFEM's GetElementTransMat, but returns MatFunctionalOnBasis
% (after row scaling) WITHOUT inverting it.  The transform used by VFEM is
% TransMat = inv(MFB); for the Goerisch quadratic form we instead want
% b = TransMat*w1 = MFB \ w1, obtained by a verified SOLVE whose enclosure width
% scales with the (small) raw-coefficient norm ||b||, not ||w1||.  This avoids the
% loose full interval inverse and realises "assemble directly on the physical
% element via closed-form Bernstein integrals".

    DegRTElt = get_DegRTElt(M);
    Nodes = I_intval(LocalNodes);
    NormVec = I_zeros(4,3);
    LocalFacetList = [2,3,4; 1,3,4; 1,2,4; 1,2,3];
    for k = 1:4
        F = LocalFacetList(k,:);
        a = Nodes(F(2),:) - Nodes(F(1),:);
        b = Nodes(F(3),:) - Nodes(F(1),:);
        last_idx = sum(1:4) - sum(F);
        e_14 = Nodes(last_idx,:) - Nodes(F(1),:);
        NormVec(k,:) = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)];
        if e_14*NormVec(k,:)' > 0
            NormVec(k,:) = - NormVec(k,:);
        end
        NormVec(k,:) = NormVec(k,:)/norm(NormVec(k,:));
    end
    K_vol = Vol;

    S1_local_index = get_dof_on_facet(M,1);
    S2_local_index = get_dof_on_facet(M,2);
    S3_local_index = get_dof_on_facet(M,3);
    S4_local_index = get_dof_on_facet(M,4);

    MatFunctionalOnBasis = I_zeros(DegRTElt,DegRTElt);
    DOF_F = get_DOF(2,M);
    DOF_Elt = get_DOF(3,M);
    DOF_Edge = M+1;

    F_idx = 0;
    S1_DOF_Indx = 1:DOF_F;
    MatFunctionalOnBasis(S1_DOF_Indx, S1_local_index) = NormVec(1,1)*eye(DOF_F,DOF_F);
    MatFunctionalOnBasis(S1_DOF_Indx, S1_local_index + DOF_Elt ) = NormVec(1,2)*eye(DOF_F,DOF_F);
    MatFunctionalOnBasis(S1_DOF_Indx, S1_local_index + DOF_Elt*2 ) = NormVec(1,3)*eye(DOF_F,DOF_F);
    [s1_idx, s4_idx] = get_common_dof_on_facets(M,1,4);
    MatFunctionalOnBasis(s1_idx, s4_idx + DOF_Elt*3 ) = (Nodes(2,:)*NormVec(1,:)')*eye(DOF_Edge,DOF_Edge);

    F_idx = F_idx + DOF_F;
    S2_DOF_Indx = (1:DOF_F)+DOF_F;
    MatFunctionalOnBasis(S2_DOF_Indx, S2_local_index )    = NormVec(2,1)*eye(DOF_F,DOF_F);
    MatFunctionalOnBasis(S2_DOF_Indx, S2_local_index+DOF_Elt ) = NormVec(2,2)*eye(DOF_F,DOF_F);
    MatFunctionalOnBasis(S2_DOF_Indx, S2_local_index+DOF_Elt*2 ) = NormVec(2,3)*eye(DOF_F,DOF_F);
    [s2_idx, s4_idx] = get_common_dof_on_facets(M,2,4);
    MatFunctionalOnBasis(s2_idx+F_idx, s4_idx + DOF_Elt*3 ) = (Nodes(3,:)*NormVec(2,:)')*eye(DOF_Edge,DOF_Edge);

    F_idx = F_idx + DOF_F;
    S3_DOF_Indx = (1:DOF_F)+DOF_F*2;
    MatFunctionalOnBasis(S3_DOF_Indx, S3_local_index )    = NormVec(3,1)*eye(DOF_F,DOF_F);
    MatFunctionalOnBasis(S3_DOF_Indx, S3_local_index+DOF_Elt ) = NormVec(3,2)*eye(DOF_F,DOF_F);
    MatFunctionalOnBasis(S3_DOF_Indx, S3_local_index+DOF_Elt*2 ) = NormVec(3,3)*eye(DOF_F,DOF_F);
    [s3_idx, s4_idx] = get_common_dof_on_facets(M,3,4);
    MatFunctionalOnBasis(s3_idx+F_idx, s4_idx + DOF_Elt*3 ) = (Nodes(4,:)*NormVec(3,:)')*eye(DOF_Edge,DOF_Edge);

    F_idx = F_idx + DOF_F;
    S4_DOF_Indx = (1:DOF_F)+DOF_F*3;
    MatFunctionalOnBasis(S4_DOF_Indx, S4_local_index )    = NormVec(4,1)*eye(DOF_F,DOF_F);
    MatFunctionalOnBasis(S4_DOF_Indx, S4_local_index+DOF_Elt ) = NormVec(4,2)*eye(DOF_F,DOF_F);
    MatFunctionalOnBasis(S4_DOF_Indx, S4_local_index+DOF_Elt*2 ) = NormVec(4,3)*eye(DOF_F,DOF_F);
    s4_idx = 1:DOF_F;
    MatFunctionalOnBasis(s4_idx+F_idx, s4_idx + DOF_Elt*3 ) = Nodes(1,:)*NormVec(4,:)'*eye(DOF_F,DOF_F);

    DOF_M_minus_1 = get_DOF(3,M-1);
    Interior_1_DOF_Indx = (1:DOF_M_minus_1)+DOF_F*4;
    Interior_2_DOF_Indx = (1:DOF_M_minus_1)+DOF_F*4+DOF_M_minus_1;
    Interior_3_DOF_Indx = (1:DOF_M_minus_1)+DOF_F*4+DOF_M_minus_1*2;

    dof_a_idx=1:DOF_Elt;
    dof_b_idx=dof_a_idx + DOF_Elt;
    dof_c_idx=dof_b_idx + DOF_Elt;
    dof_d_idx=(DOF_Elt*3+1):DegRTElt;

    A_M_MinusOne_M  = getInnerProdMatrix(M-1,M,K_vol);
    A_M_MinusOne_M_PlusOne  = getInnerProdMatrix(M-1,M+1,K_vol);

    MatDegreeUpByX = get_MatDegreeUpByX(M, LocalNodes);
    MatDegreeUpByX_Homo = MatDegreeUpByX(:, S4_local_index,:);

    MatFunctionalOnBasis(Interior_1_DOF_Indx, dof_a_idx) = A_M_MinusOne_M;
    MatFunctionalOnBasis(Interior_1_DOF_Indx, dof_d_idx) = A_M_MinusOne_M_PlusOne*MatDegreeUpByX_Homo(:,:,1);
    MatFunctionalOnBasis(Interior_2_DOF_Indx, dof_b_idx) = A_M_MinusOne_M;
    MatFunctionalOnBasis(Interior_2_DOF_Indx, dof_d_idx) = A_M_MinusOne_M_PlusOne*MatDegreeUpByX_Homo(:,:,2);
    MatFunctionalOnBasis(Interior_3_DOF_Indx, dof_c_idx) = A_M_MinusOne_M;
    MatFunctionalOnBasis(Interior_3_DOF_Indx, dof_d_idx) = A_M_MinusOne_M_PlusOne*MatDegreeUpByX_Homo(:,:,3);

    scaling = max(abs(MatFunctionalOnBasis'));
    MFB = diag( 1./scaling ) * MatFunctionalOnBasis;
end
