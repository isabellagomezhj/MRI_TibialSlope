function [n, plane_center] = find_rims_SI(pts, cns)

% Exclude pts more than 4 cm distal to top of tibia in MRI
plat_idx = pts(:,3) > max(pts(:,3))-40;  
cns_plat = cns(all(plat_idx(cns), 2), :);
pts_plat = pts(plat_idx,:);

map = zeros(size(pts, 1), 1);
map(plat_idx) = 1:sum(plat_idx);
cns_plat = map(cns_plat);

FV.vertices = pts_plat;
FV.faces = cns_plat;
smoothing_iterations = 15;

% Calculate curvature index
shapeIndex = computeKoenderinkShapeIndex(FV, smoothing_iterations);

% % Plot
% figure; hold on
% trisurf(FV.faces, FV.vertices(:, 1), FV.vertices(:, 2), FV.vertices(:, 3), shapeIndex, 'EdgeColor', 'none'); % Mesh visualization
% colormap([linspace(0, 0, 256)', linspace(0, 1, 256)', linspace(1, 0, 256)']); % Blue to Red colormap
% colorbar;
% title('Shape Index Heatmap')
% xlabel('X'); ylabel('Y'); zlabel('Z');
% axis equal;
% grid on;
% lighting phong; 
% view(3);

%% Shape Index Curvature Gradient
% Step 1: Calculate average shape index for each face
faceShapeIndex = mean(shapeIndex(FV.faces), 2);

% Step 2: Build adjacency graph
tr = triangulation(FV.faces, FV.vertices);
neighbors = tr.neighbors();

% Step 3: Calculate gradient in shape index between adjacent faces
gradients = zeros(length(FV.faces), 1);

for i = 1:length(FV.faces)
    % current face shape index
    faceSI = faceShapeIndex(i);

    % find maximum difference among its 3 neighbors and add to array
    neighborIdxs = neighbors(i,:);
    neighborIdxs = neighborIdxs(~isnan(neighborIdxs));

    if ~isempty(neighborIdxs)
        adjSI = faceShapeIndex(neighborIdxs);
        gradients(i) = max(abs(faceSI - adjSI));
    end
end

% % Plot gradient heatmap
% figure;
% title('Shape Index Gradient Heatmap');
% patch('Faces', FV.faces, 'Vertices', FV.vertices, ...
%       'FaceVertexCData', gradients, ...
%       'FaceColor', 'flat', ...
%       'EdgeColor', 'none');
% axis equal; grid on; lighting phong; view(3);
% colorbar;
% caxis([0, prctile(gradients, 99)]);
% xlabel('X'); ylabel('Y'); zlabel('Z');

%% Step 4: Identify rim based on gradient map
% gradientThreshold = prctile(gradients, 35);
gradientThreshold = 0.01;
lowGradientFaces = FV.faces(gradients < gradientThreshold, :);

% % Plot faces below gradient threshold
% figure;
% hold on;
% title('Faces Isolated Based on SI Gradient');
% trisurf(cns_plat, pts_plat(:,1), pts_plat(:,2), pts_plat(:,3), 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', [.7 .7 .7], 'FaceAlpha', 0.5); hold on
% trisurf(lowGradientFaces, FV.vertices(:, 1), FV.vertices(:, 2), FV.vertices(:, 3), ...
%         'FaceColor', 'magenta', 'EdgeColor', 'black');
% axis equal; grid on; lighting phong; view(3);
% xlabel('X'); ylabel('Y'); zlabel('Z');

% Threshold based on curvature index
% Step 1: Calculate thresholds for shape index
% shapeIndexThreshold = prctile(shapeIndex, 45);
% shapeIndexThreshold2 = prctile(shapeIndex, 75);
shapeIndexThreshold = 0.4833;
shapeIndexThreshold2 = 0.535;

% Step 2: Filter faces with at least one vertex between thresholds
midShapeIndexFaces = any(shapeIndex(FV.faces) > shapeIndexThreshold & shapeIndex(FV.faces) < shapeIndexThreshold2, 2);
% midShapeIndexFaces = any(shapeIndex(FV.faces) > shapeIndexThreshold, 2);
filteredFaces = FV.faces(midShapeIndexFaces, :);

if isempty(filteredFaces)
    error('No faces found between the 45th and 75th percentiles of shape index values.');
end

% % Plot
% figure;
% trisurf(cns_plat, pts_plat(:,1), pts_plat(:,2), pts_plat(:,3), 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', [.7 .7 .7], 'FaceAlpha', 0.5); hold on
% trisurf(filteredFaces, FV.vertices(:, 1), FV.vertices(:, 2), FV.vertices(:, 3), ...
%         'FaceColor', 'cyan', 'EdgeColor', 'black');
% % trisurf(lowGradientFaces, FV.vertices(:, 1), FV.vertices(:, 2), FV.vertices(:, 3), ...
% %         'FaceColor', 'magenta', 'EdgeColor', 'black');
% title('Faces Isolated Based on SI');
% legend({'tibia', '45-75th percentile SI faces'})
% xlabel('X');
% ylabel('Y');
% zlabel('Z');
% axis equal;
% grid on;
% lighting phong;
% view(3);

% Isolate faces that are both low-gradient and mid-shape-index
commonFaces = intersect(filteredFaces, lowGradientFaces, 'rows');

% % Plot
% figure;
% trisurf(cns_plat, pts_plat(:,1), pts_plat(:,2), pts_plat(:,3), 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', [.7 .7 .7], 'FaceAlpha', 0.5); hold on
% trisurf(commonFaces, FV.vertices(:, 1), FV.vertices(:, 2), FV.vertices(:, 3), ...
%         'FaceColor', 'blue', 'EdgeColor', 'black');
% title('Isolated Faces');
% legend({'tibia', 'combined faces'})
% xlabel('X');
% ylabel('Y');
% zlabel('Z');
% axis equal;
% grid on;
% lighting phong;
% view(3);

%% Filter out small and inlier components

% Step 1: Build adjacency graph
edges = [commonFaces(:, [1, 2]); commonFaces(:, [2, 3]); commonFaces(:, [3, 1])];
edges = sort(edges, 2);
edges = unique(edges, 'rows');

if isempty(edges)
    error('No edges could be created for the filtered faces.');
end

meshGraph = graph(edges(:, 1), edges(:, 2));
connectedComponents = conncomp(meshGraph);
numComponents = max(connectedComponents);

% Step 2: Map faces to their component
faceComponents = connectedComponents(commonFaces(:, 1));

% Step 3: Compute component areas and centroid
componentAreas = zeros(numComponents, 1);
componentCentroids = zeros(numComponents, 3);

for comp = 1:numComponents
    compFaces = commonFaces(faceComponents == comp, :);
    if ~isempty(compFaces)
        % Calculate area
        componentAreas(comp) = sum(trisurf_area(compFaces, FV.vertices));

        % Calculate centroid
        compVertexIndices = unique(compFaces(:));
        componentCentroids(comp, :) = mean(FV.vertices(compVertexIndices, :), 1);
    end
end

% Filter out small components
size_idxs = componentAreas > 30;

% Filter out components close to the center of the plateau
% nonzero_mask = any(componentCentroids, 2);
% cent = mean(componentCentroids(nonzero_mask, :));
% axial_dists = vecnorm(componentCentroids(:, 1:2) - cent(1,1:2),2,2);
% dist_idxs = axial_dists > 20;

% Filter out components close to the distal "cut"
cut_idxs = abs(componentCentroids(:,3) - min(FV.vertices(:,3))) > 15;

% final_comp_mask = size_idxs & dist_idxs & cut_idxs;
final_comp_mask = size_idxs & cut_idxs;
final_faces_mask = final_comp_mask(faceComponents);
final_rim_faces = commonFaces(final_faces_mask, :);

% % Plot
% figure;
% trisurf(cns_plat, pts_plat(:,1), pts_plat(:,2), pts_plat(:,3), 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', [.7 .7 .7], 'FaceAlpha', 0.5); hold on
% trisurf(final_rim_faces, FV.vertices(:, 1), FV.vertices(:, 2), FV.vertices(:, 3), ...
%         'FaceColor', 'green', 'EdgeColor', 'black');
% title('Intermediate Filtering of Rim Faces');

%% Fit plane to the isolated vertices
% colors = {'cyan', 'magenta', 'yellow'};

for i = 1:2 
isolated_vertices_indices = unique(final_rim_faces(:));
isolated_vertices = FV.vertices(isolated_vertices_indices, :);

[n1, ~, plane_center1] = affine_fit(isolated_vertices);

if n1(3) < 0
    n1 = -n1;
end

vectors_to_centroids = componentCentroids - plane_center1;
signed_distances = vectors_to_centroids * n1;
% center_distances = vecnorm(vectors_to_centroids,1,2);

% Remove components >4 mm below/above the plane
% plane_idxs = signed_distances > -5 & signed_distances < 2 & center_distances > 25;
plane_idxs = signed_distances > -4 & signed_distances < 4;

% final_comp_mask = size_idxs & dist_idxs & cut_idxs & plane_idxs;
final_comp_mask = size_idxs & cut_idxs & plane_idxs;
final_faces_mask = final_comp_mask(faceComponents);
final_rim_faces = commonFaces(final_faces_mask, :);

% trisurf(final_rim_faces, FV.vertices(:, 1), FV.vertices(:, 2), FV.vertices(:, 3), ...
%         'FaceColor', colors{i}, 'EdgeColor', 'black');
end

%% Fit plane to final vertices
isolated_vertices_indices = unique(final_rim_faces(:));
isolated_vertices = FV.vertices(isolated_vertices_indices, :);

[n, plane, plane_center] = affine_fit(isolated_vertices);

% %% Find tubercle cutoff point
% % project points onto plane ap axis
% axis1 = plane(:, 1);
% axis2 = plane(:, 2);
% 
% ang_axis1 = abs(dot(axis1, [0, 1, 0]));
% ang_axis2 = abs(dot(axis2, [0, 1, 0]));
% 
% % find ap axis (closest to global Y axis)
% centered_vertices = isolated_vertices - plane_center;
% if ang_axis1 > ang_axis2
%     ap_axis_intrinsic = axis1;
% else
%     ap_axis_intrinsic = axis2;   
% end
% 
% ap_coords = centered_vertices * ap_axis_intrinsic;
% 
% % check direction of axis
% if dot(ap_axis_intrinsic, [0, 1, 0]) < 0
%     ap_coords = -ap_coords;
% end
% 
% ant_prct_val = prctile(ap_coords, 20);
% [~, ant_idx] = min(abs(ap_coords - ant_prct_val));
% ant_pt = isolated_vertices(ant_idx, :);

% %% Plot
% figure;
% trisurf(cns_plat, pts_plat(:,1), pts_plat(:,2), pts_plat(:,3), 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', [.7 .7 .7], 'FaceAlpha', 0.5); hold on
% trisurf(final_rim_faces, FV.vertices(:, 1), FV.vertices(:, 2), FV.vertices(:, 3), ...
%         'FaceColor', 'cyan', 'EdgeColor', 'black');
% 
% % Plane
% [S1,S2] = meshgrid([-30 0 30]);
% X = plane_center(1)+[S1(:) S2(:)]*plane(1,:)';
% Y = plane_center(2)+[S1(:) S2(:)]*plane(2,:)';
% Z = plane_center(3)+[S1(:) S2(:)]*plane(3,:)';
% surf(reshape(X,3,3),reshape(Y,3,3),reshape(Z,3,3),'facecolor','blue','facealpha',0.5);
% plot3(plane_center(1), plane_center(2), plane_center(3), 'r .', 'MarkerSize', 30);
% quiver3(plane_center(1), plane_center(2), plane_center(3), n(1), n(2), n(3), 'r', 'AutoScaleFactor', 20);
% % plot3(ant_pt(1), ant_pt(2), ant_pt(3), 'y .', 'MarkerSize', 30);
% 
% title('Final Isolated Faces and Best-Fit Plane');
% legend({'tibia', 'filtered faces'})
% xlabel('X');
% ylabel('Y');
% zlabel('Z');
% axis equal;
% grid on;
% lighting phong;
% view(3);

end