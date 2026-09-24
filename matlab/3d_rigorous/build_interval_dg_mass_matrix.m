function M_dg_iv = build_interval_dg_mass_matrix(mesh, M)
%BUILD_INTERVAL_DG_MASS_MATRIX  Rigorous DG mass assembly by Bernstein formula.
%
% Each local mass block is computed on the physical tetrahedron as
% getInnerProdMatrix(M,M,K_vol), where K_vol is obtained from interval
% element geometry.  No reference-element transform or numerical quadrature
% is used.

global INTERVAL_MODE
old_interval_mode = INTERVAL_MODE;
INTERVAL_MODE = 1;
cleanup = onCleanup(@() restore_interval_mode(old_interval_mode));

ElementList = mesh.ElementList;
NodeList = mesh.NodeList;
NumElt = size(ElementList, 1);
DegK = get_DOF(3, M);
DimDG = NumElt * DegK;

M_dg_iv = I_sparse(DimDG, DimDG);
for e = 1:NumElt
    loc = (e-1)*DegK + (1:DegK);
    LocalNodes = I_intval(NodeList(ElementList(e,:), :));
    K_vol = get_volume(LocalNodes);
    M_dg_iv(loc, loc) = getInnerProdMatrix(M, M, K_vol);
end

if ~isa(M_dg_iv, 'intval')
    error('Interval DG mass assembly failed: M_dg is not interval-valued.');
end
end

function restore_interval_mode(old_interval_mode)
global INTERVAL_MODE
INTERVAL_MODE = old_interval_mode;
end
