function mesh = read_dolfin_mesh(filename)
% READ_DOLFIN_MESH 読み込んだ Dolfin XML メッシュファイルから mesh 構造体を生成する
%
%   mesh = read_dolfin_mesh(filename)
%
% 入力:
%   filename - 読み込む XML ファイルのパス（例: 'mesh.xml'）
%
% 出力:
%   mesh - 構造体（fields: nodes, edges, elements, domain, min_edge_length）
%          nodes: 各行が [x, y] 座標（XML中の z 成分は無視）
%          elements: 各行が 3 頂点インデックス（XML の triangle 要素；注意: XML は 0 始まりなので MATLAB 用に +1 している）
%          edges: 各行が要素から抽出したユニークな辺（各辺は昇順にソート済み）
%          domain: サンプルとして最初の三角形の頂点を [1, 3, 2] の順に格納
%          min_edge_length: 各三角形要素の3辺のうち最大の辺長の値の中で最小のもの
%
% 例:
%   mesh = read_dolfin_mesh('mesh.xml');

    % XML ファイルの読み込み
    xDoc = xmlread(filename);
    
    %% 頂点情報の抽出
    verticesList = xDoc.getElementsByTagName('vertex');
    numVertices = verticesList.getLength;
    nodes = zeros(numVertices, 2);
    
    for k = 0:numVertices-1
        vertex = verticesList.item(k);
        % x, y の属性値を文字列として取得し数値に変換
        xStr = char(vertex.getAttribute('x'));
        yStr = char(vertex.getAttribute('y'));
        nodes(k+1, :) = [str2double(xStr), str2double(yStr)];
    end
    
    %% 三角形（セル）情報の抽出
    trianglesList = xDoc.getElementsByTagName('triangle');
    numTriangles = trianglesList.getLength;
    elements = zeros(numTriangles, 3);
    
    for k = 0:numTriangles-1
        triangle = trianglesList.item(k);
        v0 = str2double(char(triangle.getAttribute('v0')));
        v1 = str2double(char(triangle.getAttribute('v1')));
        v2 = str2double(char(triangle.getAttribute('v2')));
        % XML 中のインデックスは 0 始まりなので MATLAB 用に +1
        elements(k+1, :) = [v0+1, v1+1, v2+1];
    end
    
    %% 辺情報の作成（各三角形の辺を抽出してユニークにする）
    edge_list = [];
    for k = 1:size(elements,1)
        tri = elements(k,:);
        % 各三角形は 3 辺を持つ: (v1,v2), (v2,v3), (v3,v1)
        edges_tri = [tri([1,2]); tri([2,3]); tri([3,1])];
        % 各辺は小さい方のインデックスが前になるようにソート
        edges_tri = sort(edges_tri,2);
        edge_list = [edge_list; edges_tri];
    end
    edges = unique(edge_list, 'rows');
    
    %% 各三角形要素の最大辺長を計算し，その中で最小の長さを求める
    triangle_max_lengths = zeros(numTriangles, 1);
    for k = 1:numTriangles
        tri = elements(k,:);
        % 各辺の長さを計算
        d1 = norm(nodes(tri(1),:) - nodes(tri(2),:));
        d2 = norm(nodes(tri(2),:) - nodes(tri(3),:));
        d3 = norm(nodes(tri(3),:) - nodes(tri(1),:));
        triangle_max_lengths(k) = max([d1, d2, d3]);
    end
    min_edge_length = min(triangle_max_lengths);
    
    %% domain の定義
    % ここでは例として，最初の三角形の頂点を [1, 3, 2] の順で domain とする
    %    domain = nodes([1, 3, 2], :);
    domain = [];

    bd_edge_ids = get_boundary_edge(elements, edges);
    
    %% 結果の構造体を生成
    mesh = struct('nodes', nodes, 'edges', double(edges), 'elements', double(elements), ...
                  'bd_edge_ids', bd_edge_ids, ... 
                  'domain', domain, 'min_edge_length', min_edge_length);
end

function boundary_edge_ids = get_boundary_edge(elements, edges)
% GET_BOUNDARY_EDGE_INDICES
%   boundary_edge_ids = get_boundary_edge_indices(mesh)
%
% 入力:
%   mesh - 構造体 (fields: nodes, edges, elements, domain)
%          elements: [M x 3] 行列（三角形の頂点インデックス, 1始まり）
%          edges: [E x 2] 行列（ユニークな辺, 各辺は昇順にソート済み）
%
% 出力:
%   boundary_edge_ids - mesh.edges 内で境界に属する辺の番号リスト
%
% 処理の流れ:
%   1. 各三角形から3辺を抽出し、全エッジリスト allEdges を作成
%   2. unique, accumarray を用いて各ユニーク辺の出現回数を求める
%   3. 出現回数が1の辺を境界辺と判定
%   4. mesh.edges 内の各辺が境界辺に含まれているかを調べ、その番号を返す

    % --- 1. 各三角形の辺を列挙 ---
    numElements = size(elements, 1);
    allEdges = zeros(numElements * 3, 2);
    idx = 1;
    for i = 1:numElements
        tri = elements(i,:);
        % 各三角形の辺（昇順にソート）
        edge1 = sort([tri(1), tri(2)]);
        edge2 = sort([tri(2), tri(3)]);
        edge3 = sort([tri(3), tri(1)]);
        allEdges(idx, :)   = edge1;
        allEdges(idx+1, :) = edge2;
        allEdges(idx+2, :) = edge3;
        idx = idx + 3;
    end

    % --- 2. ユニークな辺とその出現回数を取得 ---
    [uEdges, ~, ic] = unique(allEdges, 'rows');
    counts = accumarray(ic, 1);
    
    % --- 3. 出現回数が1回の辺は境界辺 ---
    boundaryEdges = uEdges(counts == 1, :);
    
    % --- 4. mesh.edges 内で境界辺に該当する番号を抽出 ---
    numMeshEdges = size(edges, 1);
    boundary_edge_ids = [];
    for j = 1:numMeshEdges
        % mesh.edges(j,:) はすでに昇順にソートされている前提
        current_edge = edges(j,:);
        if ismember(current_edge, boundaryEdges, 'rows')
            boundary_edge_ids(end+1,1) = j; %#ok<AGROW>
        end
    end
end
