%% Bilateral mesh denoising
% Based on Fleishman et al., 2003
% Fleishman, Shachar & Drori, Iddo & Cohen-Or, Daniel. (2003). 
% Bilateral Mesh Denoising. ACM Transactions on Graphics. 22. 
% 10.1145/1201775.882368. 

function smoothed_pts = bilateral_mesh_denoising(FV, sigma_c, sigma_s, iterations)
    % FV: Mesh structure with fields FV.vertices (Nx3) and FV.faces (Mx3)
    % sigma_c: Spatial smoothing parameter
    % sigma_s: Normal-based smoothing parameter
    % iterations: Number of iterations for smoothing

    smoothed_pts = FV.vertices;  % Copy initial vertex positions

    % Compute per-vertex normals
    normals = vertex_normals(FV);

    alpha = 0.1; % volume-preservation factor

    for iter = 1:iterations
        new_pts = smoothed_pts; % Store updated vertex positions
        
        for v = 1:size(FV.vertices, 1)
            % Get neighboring vertices
            neighbors = find_neighbor_vertices(FV, v);
            
            if isempty(neighbors)
                continue;
            end

            % Current vertex and normal
            v_pos = smoothed_pts(v, :);
            n = normals(v, :);

            % Initialize weight sum
            sum_h = 0;
            normalizer = 0;

            % Process each neighboring vertex
            for i = 1:length(neighbors)
                q = smoothed_pts(neighbors(i), :);

                % Compute spatial weight (distance-based)
                t = norm(v_pos - q);
                wc = exp(-t^2 / (2 * sigma_c^2));

                % Compute similarity weight (normal projection)
                h = dot(n, (q - v_pos));
                ws = exp(-h^2 / (2 * sigma_s^2));

                % Accumulate weighted displacement
                sum_h = sum_h + (wc * ws) * h;
                normalizer = normalizer + (wc * ws);
            end

            % Update vertex position
            if normalizer > 0
                new_pts(v, :) = v_pos + n * (sum_h / normalizer);
            end
        end

        % Volume preservation
        displacements = new_pts - smoothed_pts;
        for v = 1:size(FV.vertices, 1)
            new_pts(v, :) = new_pts(v, :) - alpha * dot(displacements(v, :), normals(v, :)) * normals(v, :);
        end       

        % Update vertices for next iteration
        smoothed_pts = new_pts;
    end
end

%% Helpers

function normals = vertex_normals(FV)
    % Computes per-vertex normals using face normals
    normals = zeros(size(FV.vertices));

    % Compute face normals
    face_normals = cross(FV.vertices(FV.faces(:,2), :) - FV.vertices(FV.faces(:,1), :), ...
                         FV.vertices(FV.faces(:,3), :) - FV.vertices(FV.faces(:,1), :));
    face_normals = normalize_vectors(face_normals);
    
    % Accumulate face normals for each vertex
    for i = 1:size(FV.faces, 1)
        for j = 1:3
            normals(FV.faces(i, j), :) = normals(FV.faces(i, j), :) + face_normals(i, :);
        end
    end

    % Normalize vertex normals
    normals = normalize_vectors(normals);
end

function neighbors = find_neighbor_vertices(FV, vertexIndex)
    % Extract faces containing the vertex
    connectedFaces = FV.faces(any(FV.faces == vertexIndex, 2), :);
    
    % Flatten and get unique vertex indices
    neighbors = unique(connectedFaces(:));
    
    % Remove the input vertex itself
    neighbors(neighbors == vertexIndex) = [];
end

function normalized_vectors = normalize_vectors(vectors)
    norms = sqrt(sum(vectors.^2, 2));
    norms(norms == 0) = 1;  % Avoid division by zero
    normalized_vectors = vectors ./ norms;
end

