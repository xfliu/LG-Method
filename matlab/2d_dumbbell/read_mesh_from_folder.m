function mesh = read_mesh_from_folder(folder_name)
%% Load data from specified folder to create a structure mesh.
%   mesh = read_mesh_from_folder(folder_name)
%
% 入力:
%   folder_name - The folder for mesh data.
%
% 出力:
%   mesh - 構造体（fields: nodes, edges, elements, domain, min_edge_length）
%          nodes: 各行が [x, y] 座標（XML中の z 成分は無視）
%          elements: 各行が 3 頂点インデックス（XML の triangle 要素；注意: XML は 0 始まりなので MATLAB 用に +1 している）
%          edges: 各行が要素から抽出したユニークな辺（各辺は昇順にソート済み）
%          domain: サンプルとして最初の三角形の頂点を [1, 3, 2] の順に格納
%          bd_edges: Boudary edges represented by two end node indices.
%          %min_edge_length: 各三角形要素の3辺のうち最大の辺長の値の中で最小のもの
%
% Example:
%   mesh = read_dolfin_mesh('./mesh8/');

    mesh_path = folder_name;
    vert = load([mesh_path, 'vert.dat']);
    edges = load([mesh_path, 'edge.dat']);
    tri  = load([mesh_path, 'tri.dat']);
    bd_edges   = load([mesh_path, 'bd.dat']);
    
    domain = [];

    is_edge_bd = find_is_edge_bd(edges, bd_edges, size(edges,1), size(bd_edges,1));
    bd_edge_ids =find(is_edge_bd>0)';

    %% 結果の構造体を生成
    mesh = struct('nodes', vert, 'edges', double(edges), 'elements', double(tri), 'bd_edge_ids', bd_edge_ids, ...
                  'domain', domain);
end


function is_edge_bd = find_is_edge_bd(edge, bd_edge, ne, nb)
    edge_idx                = sort(edge,    2) * [ne; 1];
    bd_edge_idx             = sort(bd_edge, 2) * [ne; 1];
    [~, bd_edge_idx]        = ismember(bd_edge_idx, edge_idx);
    is_edge_bd              = zeros(ne, 1);
    is_edge_bd(bd_edge_idx) = ones(nb, 1);
end
