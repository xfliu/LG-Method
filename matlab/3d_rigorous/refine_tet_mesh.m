function mesh_fine = refine_tet_mesh(mesh)
% Red-refinement: each tet -> 8 sub-tets by adding 6 edge midpoints.
% h_max is halved.  Produces a nested mesh (no external files needed).

EL = mesh.ElementList;  % NumElt x 4, sorted ascending
NL = mesh.NodeList;     % NumNode x 3
NumElt  = size(EL,1);
NumNode = size(NL,1);

% ---- Build edge-to-midpoint index map ----
% Encode edge (a,b) with a<b as a*big+b for fast lookup via sparse matrix
big   = NumNode + 1;
% Collect all 6 edges per tet
edges = zeros(NumElt*6, 2);
ptr   = 0;
edge_pairs = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
for e = 1:NumElt
    v = EL(e,:);
    for k = 1:6
        a = v(edge_pairs(k,1)); b = v(edge_pairs(k,2));
        ptr = ptr+1;
        edges(ptr,:) = [min(a,b), max(a,b)];
    end
end
edges = unique(edges,'rows');
NumEdge = size(edges,1);

% Midpoint coordinates
mid_nodes = 0.5*(NL(edges(:,1),:) + NL(edges(:,2),:));
% Midpoint global index = NumNode + local_edge_index
new_nodes = [NL; mid_nodes];
NewNumNode = NumNode + NumEdge;

% Sparse lookup: edge (a,b) -> midpoint global index
% Store as sparse: row=a, col=b, val = NumNode + edge_idx
S = sparse(edges(:,1), edges(:,2), NumNode + (1:NumEdge)', big, big);

% ---- Build 8 sub-tets per element ----
new_elts = zeros(NumElt*8, 4);
eptr = 0;
for e = 1:NumElt
    v1=EL(e,1); v2=EL(e,2); v3=EL(e,3); v4=EL(e,4);
    % midpoints (a<b guaranteed since EL is sorted ascending)
    m12 = S(v1,v2); m13 = S(v1,v3); m14 = S(v1,v4);
    m23 = S(v2,v3); m24 = S(v2,v4); m34 = S(v3,v4);
    eptr=eptr+1; new_elts(eptr,:) = sort([v1, m12, m13, m14]);
    eptr=eptr+1; new_elts(eptr,:) = sort([v2, m12, m23, m24]);
    eptr=eptr+1; new_elts(eptr,:) = sort([v3, m13, m23, m34]);
    eptr=eptr+1; new_elts(eptr,:) = sort([v4, m14, m24, m34]);
    eptr=eptr+1; new_elts(eptr,:) = sort([m12, m13, m14, m24]);
    eptr=eptr+1; new_elts(eptr,:) = sort([m12, m13, m23, m24]);
    eptr=eptr+1; new_elts(eptr,:) = sort([m13, m14, m24, m34]);
    eptr=eptr+1; new_elts(eptr,:) = sort([m13, m23, m24, m34]);
end

% ---- Assemble fine mesh struct ----
tmp.ElementList = new_elts(1:eptr,:);
tmp.NodeList    = new_nodes(1:NewNumNode,:);

FacetList = get_FacetList(tmp.ElementList);
tmp.FacetList   = FacetList;
NumF_new        = size(FacetList,1);
MaxEdgeNum      = NewNumNode + NumF_new + 10;
[EdgeList_new, NumEdge_new] = mesh_get_EdgeList(tmp, MaxEdgeNum);
tmp.EdgeList    = EdgeList_new;
tmp.NumF        = NumF_new;
tmp.NumNode     = NewNumNode;
tmp.NumEdge     = NumEdge_new;
tmp.NumElt      = eptr;
[F2E, E2F]      = mesh_get_Facet2Element(tmp);
tmp.Facet2Element = F2E;
tmp.Element2Facet = E2F;

mesh_fine = tmp;

% Compute and report hmax
hmax = 0;
NL2  = tmp.NodeList;
EL2  = tmp.ElementList;
for e = 1:min(eptr, 500)   % sample for speed
    verts = NL2(EL2(e,:),:);
    for i=1:4; for j=i+1:4
        hmax = max(hmax, norm(verts(i,:)-verts(j,:)));
    end; end
end
fprintf('  Refined mesh: %d nodes, %d tets, hmax~%.6f\n', NewNumNode, eptr, hmax);
end
