%% Function to Find Neighboring Vertices
function neighbors = find_neighbor_vertices(FV, vertexIndex)
    % Extract faces that contain the given vertex
    connectedFaces = FV.faces(any(FV.faces == vertexIndex, 2), :);
    % Flatten the list of vertices in these faces
    neighbors = unique(connectedFaces(:));
    % Remove the input vertex itself from the neighbor list
    neighbors(neighbors == vertexIndex) = [];
end