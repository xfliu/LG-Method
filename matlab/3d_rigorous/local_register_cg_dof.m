function [L2G, DOF_BD, DimCG] = local_register_cg_dof(mesh, p)
%LOCAL_REGISTER_CG_DOF  Conforming CG^p DOF registration on a tet mesh.
% Returns the local-to-global map L2G (NumElt x DegK), the list of Dirichlet
% boundary DOFs DOF_BD, and the global dimension DimCG.  Entity-based keys make
% it robust across meshes (used by both regular and T_F rigorous drivers).

    ElementList = mesh.ElementList;
    FacetList = mesh.FacetList;
    Facet2Element = mesh.Facet2Element;
    NumElt = size(ElementList, 1);
    ijkl = get_IJKL(p);
    DegK = size(ijkl, 1);
    boundary_facets = sort(FacetList(Facet2Element(:,2) == 0, :), 2);

    map = containers.Map('KeyType', 'char', 'ValueType', 'int32');
    bd_cache = false(NumElt * DegK, 1);
    L2G = zeros(NumElt, DegK);
    next_id = 0;

    for e = 1:NumElt
        verts_e = ElementList(e, :);
        for a = 1:DegK
            alpha = ijkl(a, :);
            support = find(alpha > 0);
            verts = verts_e(support);
            weights = alpha(support);
            key = make_dof_key(e, verts, weights);
            if isKey(map, key)
                gid = map(key);
            else
                next_id = next_id + 1;
                gid = next_id;
                map(key) = gid;
                if is_boundary_entity(verts, boundary_facets)
                    bd_cache(gid) = true;
                end
            end
            L2G(e, a) = gid;
        end
    end

    DimCG = next_id;
    bd_cache = bd_cache(1:DimCG);
    DOF_BD = find(bd_cache);
end

function key = make_dof_key(e, verts, weights)
    if numel(verts) == 4
        key = sprintf('c:%d:%d:%d:%d:%d', e, weights(1), weights(2), weights(3), weights(4));
    else
        data = [verts(:)'; weights(:)'];
        key = sprintf('s:%d:', numel(verts));
        key = [key, sprintf('%d,', data(:))];
    end
end

function tf = is_boundary_entity(verts, boundary_facets)
    if numel(verts) == 4
        tf = false; return;
    end
    tf = false;
    for i = 1:size(boundary_facets, 1)
        if all(ismember(verts, boundary_facets(i, :)))
            tf = true; return;
        end
    end
end
