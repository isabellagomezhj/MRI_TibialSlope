%% Calculate lateral and medial tibial plateau slope from short knee MRI
% Inputs:
% fld = path to tibia .iv file
% NOTE: the tibia should roughly be aligned with the Z-axis 
% superior(+)-inferior(-) and the Y-axis anterior(-)-posterior(+)
% kneeside = subject laterality
% plt_figs = whether to plot tibia points and plateau normal vectors (bool)

% Outputs:
% slopes = 1 x 3 table with the full, medial, and lateral plateau slopes in
% the sagittal plane
% NOTE: negative indicates a posterior-inferior directed slope ("higher
% slope")

%% User inputs
% path to tibia .iv file
fld = "S:\BiomechanicsResearch\groupImhauser\Modeling\UVM_test\Data_Reduced\FemalePair7\Case\model_inputs\Geometries\Tib_high_density.iv";
kneeside = 'R';
plt_figs = 0;

%% Main script
slopes = zeros(1,3);

% auto ACS calculation
tibacs = getTibiaACSv7(fld, kneeside);

% import tibia points
pts = read_vrml_fast(fld);

% find plateau normal
[mri_sem_ax, t, ~, ~, ~] = ellipse_ax(pts, 5, 100);

% transform so tibia axes are roughly oriented along MRI axes and we can
% find the rims
temp_locpts = inv(t)*[pts'; ones(1,size(pts,1))];
temp_locpts = temp_locpts(1:3,:)';

% find vector normal to full plateau, medial and lateral plateaus
[~, ~, n, ~, med_n, lat_n] = find_rims(temp_locpts, kneeside, mri_sem_ax);

% transform back to MRI coordinates
n_mri = t(1:3,1:3) * n;
med_n_mri = t(1:3,1:3) * med_n;
lat_n_mri = t(1:3,1:3) * lat_n;

% plot CS and plateau normals
if plt_figs == 1
    figure;
    plot3(pts(:,1), pts(:,2), pts(:,3), '.', 'MarkerSize', .3); hold on
    quiver3(tibacs(1,4), tibacs(2,4), tibacs(3,4), n_mri(1), n_mri(2), n_mri(3), 100, 'LineWidth', 5);
    quiver3(tibacs(1,4), tibacs(2,4), tibacs(3,4), med_n_mri(1), med_n_mri(2), med_n_mri(3), 100, 'LineWidth', 5);
    quiver3(tibacs(1,4), tibacs(2,4), tibacs(3,4), lat_n_mri(1), lat_n_mri(2), lat_n_mri(3), 100, 'LineWidth', 5);
    plotcs_colour(tibacs,50,[1 0 0;0 1 0;0 0 1]); % acs
    legend({'tibia', 'full plateau normal', 'medial plateau normal', 'lateral plateau normal'});
    xlabel('X'); ylabel('Y'); zlabel('Z');
end

% calc sagittal slopes
n_anatomic = inv(tibacs(1:3,1:3)) * n_mri;
full_slope = atan2d(n_anatomic(3), n_anatomic(1));
slopes(1,1) = full_slope;

n_med_anatomic = inv(tibacs(1:3,1:3)) * med_n_mri;
med_slope = atan2d(n_med_anatomic(3), n_med_anatomic(1));
slopes(1,2) = med_slope;

n_lat_anatomic = inv(tibacs(1:3,1:3)) * lat_n_mri;
lat_slope = atan2d(n_lat_anatomic(3), n_lat_anatomic(1));
slopes(1,3) = lat_slope;

slopes = table(slopes(:,1), slopes(:,2), slopes(:,3), 'VariableNames', {'full_plateau', 'medial', 'lateral'});

fprintf('Full plateau slope: %.2fd | Medial plateau slope: %.2fd | Lateral plateau slope: %.2fd \n', full_slope, med_slope, lat_slope);
