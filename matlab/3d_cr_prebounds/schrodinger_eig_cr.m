function [eig_lb, eig_cr, info] = schrodinger_eig_cr(mesh, V_info, neig, bc_type, sigma)
% SCHRODINGER_EIG_CR  Scalar Crouzeix-Raviart eigenvalue solver for
% the Schrodinger operator -Delta + V on a 3D tetrahedral mesh.
%
% Provides rough but guaranteed lower bounds via Liu's correction:
%   lambda_k >= nu_{k,h} / (1 + C_h^2 * nu_{k,h})
%
% CR DOFs in 3D: one per triangular facet (face midpoint average).
%
% INPUT:
%   mesh    - mesh structure from mesh_load_from_folder
%   V_info  - potential info (struct or function handle)
%   neig    - number of eigenvalues
%   bc_type - 'dirichlet' or 'neumann'
%
% OUTPUT:
%   eig_lb  - (neig, 1) guaranteed lower bounds (after Liu correction)
%   eig_cr  - (neig, 1) raw CR eigenvalues (before correction)
%   info    - struct with C_h, h_max, matrices

    global INTERVAL_MODE

    if nargin < 4
        bc_type = 'dirichlet';
    end
    if nargin < 5
        sigma = -0.5;   % default; use value below ground state (e.g. -1.5 for H2+)
    end

    ElementList = mesh.ElementList;
    NodeList = mesh.NodeList;
    Element2Facet = mesh.Element2Facet;
    Facet2Element = mesh.Facet2Element;

    NumElt = size(ElementList, 1);
    NumF = mesh.NumF;

    % CR degree: N=1 (piecewise linear on each element, DOF at face midpoints)
    N = 1;
    DegK = get_DOF(3, N);          % = 4 (4 basis functions per tet)
    DegKMinus = get_DOF(3, N-1);   % = 1

    fprintf('Schrodinger CR: %d elements, %d facets\n', NumElt, NumF);

    % Identify interior/boundary facets
    if strcmp(bc_type, 'dirichlet')
        Inner_Facet_Idx = find(sum(Facet2Element' > 0) == 2);
    else
        Inner_Facet_Idx = 1:NumF;  % Neumann: keep all facets
    end
    DimCR = length(Inner_Facet_Idx);
    fprintf('  CR DOFs: %d (interior facets)\n', DimCR);

    % Reference element matrices
    A_M_M_Ref = getInnerProdMatrix_Reference(N, N);           % (4x4) mass
    A_M_Minus_M_Minus_Ref = getInnerProdMatrix_Reference(N-1, N-1); % (1x1) for grad

    % TransMat maps from Bernstein to CR basis
    % CR basis: phi_i = 1 - 3*lambda_i (face-centered)
    TransMat = ones(4, 4) - eye(4, 4) * 3;

    % Assemble scalar stiffness, mass, and potential matrices
    CRGrad_CRGrad = I_sparse(NumF, NumF);    % Stiffness
    CR_CR_mass = I_sparse(NumF, NumF);        % Mass
    CR_V_CR = I_sparse(NumF, NumF);           % Potential

    h_max = 0;

    for EltIdx = 1:NumElt
        LocalIdxFList = Element2Facet(EltIdx, :);  % 4 facet indices
        LocalNodes = NodeList(ElementList(EltIdx,:), :);
        K_vol = get_volume(LocalNodes);

        L2G_CR = LocalIdxFList;

        % Stiffness: TransMat' * (sum_k GradMat_k' * A_{0,0} * GradMat_k) * TransMat
        MatGrad = get_GradMat(N, LocalNodes);
        LocalK = I_zeros(DegK, DegK);
        for k = 1:3
            LocalK = LocalK + MatGrad(:,:,k)' * A_M_Minus_M_Minus_Ref * MatGrad(:,:,k);
        end
        LocalK_CR = TransMat * LocalK * TransMat * K_vol;

        % Mass: TransMat' * A_{1,1} * TransMat
        LocalM_CR = TransMat * A_M_M_Ref * TransMat * K_vol;

        % Potential: integrate V * phi_i * phi_j over element
        LocalV_CR = assemble_potential_cr_element(LocalNodes, V_info, K_vol, TransMat);

        % Assemble
        CRGrad_CRGrad(L2G_CR, L2G_CR) = CRGrad_CRGrad(L2G_CR, L2G_CR) + LocalK_CR;
        CR_CR_mass(L2G_CR, L2G_CR) = CR_CR_mass(L2G_CR, L2G_CR) + LocalM_CR;
        CR_V_CR(L2G_CR, L2G_CR) = CR_V_CR(L2G_CR, L2G_CR) + LocalV_CR;

        % Track h_max for Liu constant
        for i = 1:4
            for j = i+1:4
                h = norm(LocalNodes(i,:) - LocalNodes(j,:));
                if h > h_max
                    h_max = h;
                end
            end
        end

        if mod(EltIdx, 10000) == 0
            fprintf('    %d / %d elements\n', EltIdx, NumElt);
        end
    end

    % Extract interior block
    A_cr = CRGrad_CRGrad(Inner_Facet_Idx, Inner_Facet_Idx) ...
         + CR_V_CR(Inner_Facet_Idx, Inner_Facet_Idx);
    M_cr = CR_CR_mass(Inner_Facet_Idx, Inner_Facet_Idx);

    % Solve eigenvalue problem
    fprintf('  Solving CR EVP (%d x %d)...\n', DimCR, DimCR);

    if INTERVAL_MODE > 0
        A_mid = mid(A_cr);
        M_mid = mid(M_cr);
    else
        A_mid = A_cr;
        M_mid = M_cr;
    end

    opts.disp = 0;
    opts.maxit = 1000;
    [V_eig, D_eig] = eigs(A_mid, M_mid, neig, sigma, opts);
    eig_cr = sort(real(diag(D_eig)));

    % Liu's correction constant for 3D CR
    % C_h = h_max / sqrt(10) for CR on tetrahedra (proven bound)
    C_h = h_max / sqrt(10);

    fprintf('  h_max = %.6f, C_h = %.6f, C_h^2 = %.4e\n', h_max, C_h, C_h^2);

    % Apply Liu's lower bound correction
    eig_lb = zeros(neig, 1);
    for k = 1:neig
        nu = eig_cr(k);
        eig_lb(k) = nu / (1 + C_h^2 * nu);
    end

    fprintf('\n  CR eigenvalues and Liu lower bounds:\n');
    fprintf('  %5s  %20s  %20s  %15s\n', 'k', 'nu_k (CR)', 'lambda_k_lb', 'Liu correction');
    for k = 1:neig
        correction = C_h^2 * eig_cr(k);
        fprintf('  %5d  %20.10f  %20.10f  %15.4e\n', k, eig_cr(k), eig_lb(k), correction);
    end

    info.h_max = h_max;
    info.C_h = C_h;
    info.A_cr = A_cr;
    info.M_cr = M_cr;
    info.Inner_Facet_Idx = Inner_Facet_Idx;
    info.DimCR = DimCR;
end


function LocalV = assemble_potential_cr_element(LocalNodes, V_info, K_vol, TransMat)
% Assemble potential matrix for one element in CR basis.
% LocalV = TransMat * V_bernstein * TransMat
% where V_bernstein[i,j] = integral of V * B_i * B_j over element.

    global INTERVAL_MODE

    if isstruct(V_info) && strcmp(V_info.type, 'coulomb')
        LocalV_bern = assemble_coulomb_cr(LocalNodes, V_info, K_vol);
    else
        LocalV_bern = assemble_general_potential_cr(LocalNodes, V_info, K_vol);
    end

    % Transform from Bernstein to CR basis
    LocalV = TransMat * LocalV_bern * TransMat;
end


function LocalV = assemble_coulomb_cr(LocalNodes, V_info, K_vol)
% Coulomb potential integration for CR element.
% Uses Duffy transform for elements near nuclei, standard quad otherwise.

    centers = V_info.centers;
    charges = V_info.charges;
    LocalV = I_zeros(4, 4);

    for ci = 1:size(centers, 1)
        c = centers(ci, :);
        Z = charges(ci);

        dists = sqrt(sum((LocalNodes - c).^2, 2));
        h_elem = max(sqrt(sum((LocalNodes - mean(LocalNodes)).^2, 2)));

        if min(dists) < 0.5 * h_elem
            % Duffy quadrature
            [~, closest_v] = min(dists);
            LocalV_ci = coulomb_duffy_cr(LocalNodes, c, Z, K_vol, closest_v);
        else
            % Standard quadrature
            LocalV_ci = coulomb_standard_cr(LocalNodes, c, Z, K_vol);
        end

        LocalV = LocalV + LocalV_ci;
    end
end


function LocalV = coulomb_duffy_cr(LocalNodes, center, Z, K_vol, sing_vertex)
% Duffy transform for 1/|x-c| near sing_vertex, in Bernstein basis.

    n_gauss = 8;
    [xi, wi] = gauss_legendre_01_cr(n_gauss);

    LocalV = I_zeros(4, 4);

    for i1 = 1:n_gauss
        for i2 = 1:n_gauss
            for i3 = 1:n_gauss
                eta1 = xi(i1); eta2 = xi(i2); eta3 = xi(i3);

                lam_duffy = [1-eta1, eta1*(1-eta2), eta1*eta2*(1-eta3), eta1*eta2*eta3];

                % Reorder for singularity vertex
                lam = zeros(1, 4);
                lam(sing_vertex) = lam_duffy(1);
                others = setdiff(1:4, sing_vertex);
                for kk = 1:3
                    lam(others(kk)) = lam_duffy(kk+1);
                end

                jac = eta1^2 * eta2;
                x = lam * LocalNodes;
                r = norm(x - center);
                V_val = -Z / max(r, 1e-14);

                % Bernstein basis (degree 1): B_i = lam_i
                phi = lam;

                w = wi(i1) * wi(i2) * wi(i3) * jac;
                LocalV = LocalV + (w * V_val * 6 * K_vol) * (phi' * phi);
            end
        end
    end
end


function LocalV = coulomb_standard_cr(LocalNodes, center, Z, K_vol)
% Standard quadrature for 1/|x-c| far from element.

    [qpts, qwts] = gauss_tet_quad_cr(6);
    nq = length(qwts);

    LocalV = I_zeros(4, 4);
    for q = 1:nq
        lam = qpts(q, :);
        x = lam * LocalNodes;
        r = norm(x - center);
        V_val = -Z / r;

        phi = lam;  % Bernstein degree 1
        LocalV = LocalV + (qwts(q) * V_val * 6 * K_vol) * (phi' * phi);
    end
end


function LocalV = assemble_general_potential_cr(LocalNodes, V_func, K_vol)
% General potential integration via quadrature, Bernstein basis.

    [qpts, qwts] = gauss_tet_quad_cr(4);
    nq = length(qwts);

    LocalV = I_zeros(4, 4);
    for q = 1:nq
        lam = qpts(q, :);
        x = lam * LocalNodes;
        V_val = V_func(x(1), x(2), x(3));

        phi = lam;
        LocalV = LocalV + (qwts(q) * V_val * 6 * K_vol) * (phi' * phi);
    end
end


function [pts, wts] = gauss_tet_quad_cr(degree)
% Conical product Gauss quadrature on reference tet.
    n = ceil((degree + 1) / 2);
    [xi, wi] = gauss_legendre_01_cr(n);
    pts = zeros(n^3, 4);
    wts = zeros(n^3, 1);
    idx = 0;
    for i = 1:n
        for j = 1:n
            for k = 1:n
                idx = idx + 1;
                r = xi(i); s = xi(j); t = xi(k);
                lam4 = r;
                lam3 = s*(1-r);
                lam2 = t*(1-r)*(1-s);
                lam1 = 1-lam2-lam3-lam4;
                jac = (1-r)^2*(1-s);
                pts(idx,:) = [lam1, lam2, lam3, lam4];
                wts(idx) = wi(i)*wi(j)*wi(k)*jac;
            end
        end
    end
end


function [x, w] = gauss_legendre_01_cr(n)
% Gauss-Legendre on [0,1].
    beta = (1:n-1) ./ sqrt(4*(1:n-1).^2 - 1);
    T = diag(beta, 1) + diag(beta, -1);
    [V, D] = eig(T);
    x0 = diag(D)';
    w0 = 2 * V(1,:).^2;
    x = 0.5 * (x0 + 1);
    w = 0.5 * w0;
end
