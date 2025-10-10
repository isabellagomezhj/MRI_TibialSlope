function shapeIndex = computeKoenderinkShapeIndex(FV, numIterations)
% Compute the Koenderink shape index from a mesh structure FV
% using the Rucinskiewicz curvature algorithm.
%
% INPUTS:
%   FV             - Mesh struct with fields: FV.vertices and FV.faces
%   numIterations  - (Optional) Number of smoothing iterations [default = 5]
%
% OUTPUT:
%   shapeIndex     - Nx1 array of Koenderink shape index per vertex
    if nargin < 2
        numIterations = 5;  % Default smoothing iterations
    end

    % Step 1: Face normals
    N = CalcFaceNormals(FV);

    % Step 2: Vertex normals
    [VertexNormals, Avertex, Acorner, up, vp] = CalcVertexNormals(FV, N);

    % Step 3: Curvature tensors
    [~, VertexSFM, ~] = CalcCurvature(FV, VertexNormals, N, Avertex, Acorner, up, vp);

    % Step 4: Principal curvatures
    [PrincipalCurvatures, ~, ~] = getPrincipalCurvatures(FV, VertexSFM, up, vp);

    % Step 5: Koenderink shape index formula
    k1 = PrincipalCurvatures(1, :);
    k2 = PrincipalCurvatures(2, :);
    shapeIndex = (atan(k2 ./ k1) + pi) / (2 * pi);

    % Step 6: Smoothing shape index across neighbors
    smoothedValues = shapeIndex;
    for i = 1:numIterations
        for v = 1:size(FV.vertices, 1)
            neighbors = find_neighbor_vertices(FV, v); % Returns vertex indices
            smoothedValues(v) = mean(shapeIndex(neighbors));
        end
    end

    shapeIndex = smoothedValues;  % Final smoothed shape index
end
