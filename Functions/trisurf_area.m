%% Function to Calculate Face Areas
function areas = trisurf_area(faces, vertices)
    % Calculates the area of triangular faces
    v1 = vertices(faces(:, 1), :);
    v2 = vertices(faces(:, 2), :);
    v3 = vertices(faces(:, 3), :);
    areas = 0.5 * sqrt(sum(cross(v2 - v1, v3 - v1, 2).^2, 2));
end


