function mesh = special_tet_make_mesh(n_refine)
% Mesh for the fundamental tetrahedron T_F from special_tetrahedron.md.
%
% T_F = conv{(0,0,0), (0,0,1), (1/2,1/2,1/2), (-1/2,1/2,1/2)}.

if nargin < 1
    n_refine = 0;
end

vfem_root = '/home/xfliu/Workspace/New_Project_SchrodingerEVP/VFEM_LIB/VFEM3D';
addpath(vfem_root);
addpath(genpath(fullfile(vfem_root, 'VFEM_Lib')));

NodeList = [
     0.0, 0.0, 0.0;
     0.0, 0.0, 1.0;
     0.5, 0.5, 0.5;
    -0.5, 0.5, 0.5
];
ElementList = sort([1, 2, 3, 4], 2);

mesh = complete_mesh_struct(NodeList, ElementList);
for r = 1:n_refine
    mesh = fast_refine_tet_mesh(mesh);
end
end

function mesh_fine = fast_refine_tet_mesh(mesh)
EL = mesh.ElementList;
NL = mesh.NodeList;
NumElt = size(EL, 1);
NumNode = size(NL, 1);

edge_pairs = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
edges = zeros(NumElt * 6, 2);
for k = 1:6
    rows = (k:6:NumElt*6)';
    edges(rows,:) = sort(EL(:, edge_pairs(k,:)), 2);
end
[edges_unique, ~, edge_idx] = unique(edges, 'rows');
mid_nodes = 0.5 * (NL(edges_unique(:,1),:) + NL(edges_unique(:,2),:));
new_nodes = [NL; mid_nodes];
mid = reshape(NumNode + edge_idx, 6, NumElt)';

new_elts = zeros(NumElt * 8, 4);
for e = 1:NumElt
    v1=EL(e,1); v2=EL(e,2); v3=EL(e,3); v4=EL(e,4);
    m12=mid(e,1); m13=mid(e,2); m14=mid(e,3);
    m23=mid(e,4); m24=mid(e,5); m34=mid(e,6);
    base = (e-1)*8;
    new_elts(base+1,:) = sort([v1, m12, m13, m14]);
    new_elts(base+2,:) = sort([v2, m12, m23, m24]);
    new_elts(base+3,:) = sort([v3, m13, m23, m34]);
    new_elts(base+4,:) = sort([v4, m14, m24, m34]);
    new_elts(base+5,:) = sort([m12, m13, m14, m24]);
    new_elts(base+6,:) = sort([m12, m13, m23, m24]);
    new_elts(base+7,:) = sort([m13, m14, m24, m34]);
    new_elts(base+8,:) = sort([m13, m23, m24, m34]);
end

mesh_fine = complete_mesh_struct(new_nodes, new_elts);
fprintf('  Refined mesh: %d nodes, %d tets, hmax~%.6f\n', ...
        mesh_fine.NumNode, mesh_fine.NumElt, estimate_hmax(mesh_fine));
end

function hmax = estimate_hmax(mesh)
hmax = 0;
EL = mesh.ElementList;
NL = mesh.NodeList;
for e = 1:min(size(EL,1), 500)
    verts = NL(EL(e,:),:);
    for i = 1:4
        for j = i+1:4
            hmax = max(hmax, norm(verts(i,:) - verts(j,:)));
        end
    end
end
end

function mesh = complete_mesh_struct(NodeList, ElementList)
ElementList = sort(ElementList, 2);
FacetList = get_FacetList(ElementList);

mesh.ElementList = ElementList;
mesh.NodeList = NodeList;
mesh.FacetList = FacetList;

mesh.NumF = size(FacetList, 1);
mesh.NumNode = size(NodeList, 1);
mesh.NumElt = size(ElementList, 1);

edges = [ElementList(:,[1 2]); ElementList(:,[1 3]); ElementList(:,[1 4]); ...
         ElementList(:,[2 3]); ElementList(:,[2 4]); ElementList(:,[3 4])];
mesh.EdgeList = unique(sort(edges, 2), 'rows');
mesh.NumEdge = size(mesh.EdgeList, 1);

[mesh.Facet2Element, mesh.Element2Facet] = fast_facet2element(ElementList, FacetList);
end

function [Facet2Element, Element2Facet] = fast_facet2element(ElementList, FacetList)
NumElt = size(ElementList, 1);
NumF = size(FacetList, 1);
local_facets = [2 3 4; 1 3 4; 1 2 4; 1 2 3];
all_facets = zeros(NumElt * 4, 3);
elt_ids = repelem((1:NumElt)', 4);
loc_ids = repmat((1:4)', NumElt, 1);
for lf = 1:4
    rows = (lf:4:NumElt*4)';
    all_facets(rows,:) = sort(ElementList(:, local_facets(lf,:)), 2);
end
[~, idx] = ismember(all_facets, sort(FacetList, 2), 'rows');

Element2Facet = reshape(idx, 4, NumElt)';
Facet2Element = zeros(NumF, 2);
for r = 1:numel(idx)
    f = idx(r);
    if Facet2Element(f,1) == 0
        Facet2Element(f,1) = elt_ids(r);
    else
        Facet2Element(f,2) = elt_ids(r);
    end
end
end
