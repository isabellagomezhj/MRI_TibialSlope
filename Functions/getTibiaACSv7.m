%% Tibia auto ACS function
% 7/31/2025
% Isabella Gomez Hjerthen

function [tibacs, antdist, tib_length] = getTibiaACSv7(fld, kneeside, antcorr)
tic
scale = 0.82; % UVM
mid_scale = 0.7; % UVM
% scale = 0.55;
% mid_scale = 0.43;
slice = 5;
plot_figs = 0;

% tibia in MRI space
[pts,cns] = read_vrml_fast(fld);
cns = cns(:, 1:3) + 1;

if plot_figs == 1
    % plotting original pts and MRI coordinate system
    figure;
    plotcs_colour([1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1],50,[1 0 0;0 1 0;0 0 1]); hold on
    plot3(pts(:,1), pts(:,2), pts(:,3), '.', 'MarkerSize', .3);
    xlabel('X'); ylabel('Y'); zlabel('Z');
    view(3);
end

%% Find vector normal to plateau using curvature index
[n, plane_center] = find_rims_SI(pts, cns);

% Transform tibia points so plateau normal is the z-axis
temp_y = unit(cross(n, [1 0 0]));
temp_x = unit(cross(temp_y, n));

t_med = [temp_x' temp_y' n];
med_pts = pts * t_med;

ctr = plane_center * t_med;
rim_z = ctr(3);

if plot_figs == 1
    % plot points and plateau normal vector
    figure;
    % plotcs_colour(t,50,[1 0 0;0 1 0;0 0 1]); hold on
    plot3(pts(:,1), pts(:,2), pts(:,3), '.', 'MarkerSize', .3); hold on
    plot3(med_pts(:,1), med_pts(:,2), med_pts(:,3), '.', 'MarkerSize', .3);
    legend({'original', 'plateau-normal'});
    quiver3(0,0,0, n(1), n(2), n(3), 50);
    xlabel('X'); ylabel('Y'); zlabel('Z');
end

% figure;
% plot3(med_pts(:,1), med_pts(:,2), med_pts(:,3), '.', 'MarkerSize', .3); hold on

%% Find most contoured ellipse
[~, t, antdist1, ~, ~] = ellipse_ax(med_pts, slice, rim_z);

% transform so ellipse axes are AP/ML
med_pts = inv(t)*[med_pts'; ones(1,size(med_pts,1))];
med_pts = med_pts(1:3,:)';


%% Mesh smoothing for noise reduction
FV.vertices = med_pts;
FV.faces = cns;

% smoothed_pts = med_pts; % copy med_pts coordinates
% numIterations = 10;
% for i = 1:numIterations
%     for v = 1:size(med_pts, 1)
%         neighbors = find_neighbor_vertices(FV, v); % Function to find neighbors
%         smoothed_pts(v) = mean(smoothed_pts(neighbors));
%     end
% end

% smoothed_pts = taubin_smooth(FV, 0.5, -0.5, 3);

sigma_c = 50;  % Controls spatial smoothness
sigma_s = 0.5;  % Controls feature preservation
iterations = 15;  % Number of smoothing passes

smoothed_pts = bilateral_mesh_denoising(FV, sigma_c, sigma_s, iterations);

% scale pts back to original size
% orig_range = max(FV.vertices) - min(FV.vertices);
% final_range = max(smoothed_pts) - min(smoothed_pts);
% scale_factors = orig_range ./ final_range;
% smoothed_pts = smoothed_pts .* scale_factors;

if plot_figs == 1
    idx = smoothed_pts(:,3) < 0;
    figure;
    plot3(smoothed_pts(idx,1), smoothed_pts(idx,2), smoothed_pts(idx,3), '.', 'MarkerSize', .3); hold on
    plot3(med_pts(:,1), med_pts(:,2), med_pts(:,3), '.', 'MarkerSize', .3);
    legend({'smoothed', 'original'});
    xlabel('X'); ylabel('Y'); zlabel('Z');

    figure;
    subplot(1,2,1);
    plot3(smoothed_pts(:,1), smoothed_pts(:,2), smoothed_pts(:,3), '.', 'MarkerSize', .3); hold on
    title('smoothed');
    subplot(1,2,2);
    plot3(med_pts(:,1), med_pts(:,2), med_pts(:,3), '.', 'MarkerSize', .3); hold on
    title('original');
end

% re-define med_pts so that the shaft is smoothed but plateau is preserved
cut = antdist1 * -0.5;
med_pts = [med_pts(med_pts(:,3) > cut,:); smoothed_pts(smoothed_pts(:,3) < cut,:)];

if plot_figs == 1
    figure;
    plot3(med_pts(:,1), med_pts(:,2), med_pts(:,3), '.', 'MarkerSize', .3);
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title('Tibia points with smoothed shaft');
end

%% Sagittal plane angle of posterior trace of shaft wrt plateau normal
if nargin < 3
    [ang, medn_ax, antdist, tpt] = post_trace_angle(med_pts, slice, scale, 0, plot_figs);
else
    [ang, medn_ax, antdist, tpt] = post_trace_angle_antcorr(med_pts, slice, scale, 0, plot_figs, antcorr);
end

Rx = [1 0 0; 0 cos(deg2rad(ang)) -sin(deg2rad(ang)); 0 sin(deg2rad(ang)) cos(deg2rad(ang))];

rot_pts = med_pts * Rx; % rotate points by theta

if plot_figs == 1
    figure;
    plot3(rot_pts(:,1), rot_pts(:,2), rot_pts(:,3), '. k', 'MarkerSize', 0.3); hold on
    plot3(med_pts(:,1), med_pts(:,2), med_pts(:,3), '.', 'MarkerSize', 0.3);
end

%% Top point for tubercle exclusion
% find most anterior point on each plateau rim
[med_pt, lat_pt, ~, ~, ~, ~] = find_rims(med_pts, kneeside, medn_ax);
tub_tpt = mean([med_pt; lat_pt]);
tub_tpt = tub_tpt * Rx;

% % rotate point found with shape index
% tub_tpt_med = ant_pt * t_med;
% tub_tpt_temp = inv(t) * [tub_tpt_med'; 1];
% tub_tpt_temp = tub_tpt_temp(1:3)';
% tub_tpt = tub_tpt_temp * Rx;

% tub_tpt = [0,antdist*0.4,0];

%% Tubercle exclusion

% translate tibia so all pts anterior to tub_tpt have positive y-values
tub_pts = rot_pts - tub_tpt;

notub_idx = tub_pts(:,2) < 0;
notub = rot_pts(notub_idx, :);

if plot_figs == 1
    % plot axis used to exclude tubercule
    % figure;
    plot3(rot_pts(:,1), rot_pts(:,2), rot_pts(:,3), '.', 'MarkerSize', .3); hold on
    % plot3(vertring(:,1), vertring(:,2), vertring(:,3), '.', 'MarkerSize', .3);
    % plot3([tub_tpt(:,1) tub_bpt(:,1)], [tub_tpt(:,2) tub_bpt(:,2)], [tub_tpt(:,3) tub_bpt(:,3)], 'r');
    % plot3(tub_bpt(:,1), tub_bpt(:,2), tub_bpt(:,3), '. r', 'MarkerSize', 30);
    plot3(tub_tpt(:,1), tub_tpt(:,2), tub_tpt(:,3), '. r', 'MarkerSize', 30);

    % figure;
    plot3(rot_pts(:,1), rot_pts(:,2), rot_pts(:,3), '.', 'MarkerSize', .3); hold on
    plot3(notub(:,1), notub(:,2), notub(:,3), '.', 'MarkerSize', .3);
    xlabel('x'); ylabel('y');
    % plotcs_colour([1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1],50,[1 0 0;0 1 0;0 0 1]);
end

%% Re-define SI axis
tpt_rot = tpt * Rx;

% bottom point of SI axis
if nargin < 3
    vertpt = tpt_rot(:,3) - scale*antdist;
else
    vertpt = tpt_rot(:,3) - scale*antdist + antcorr;
end
vertring_idx = (notub(:,3) > (vertpt - slice / 2)) & (notub(:,3) < (vertpt + slice / 2));
vertring = notub(vertring_idx, :);

bpt_rot = [(min(vertring(:,1)) + max(vertring(:,1)))/2 (min(vertring(:,2)) + max(vertring(:,2)))/2 (min(vertring(:,3)) + max(vertring(:,3)))/2];


% mid point for SI axis def
if nargin < 3
    vertpt = tpt_rot(:,3) - mid_scale*antdist;
else
    vertpt = tpt_rot(:,3) - mid_scale*antdist + antcorr;
end
vertring_idx = (notub(:,3) > (vertpt - slice / 2)) & (notub(:,3) < (vertpt + slice / 2));
vertring2 = notub(vertring_idx, :);

mpt_rot = [(min(vertring2(:,1)) + max(vertring2(:,1)))/2 (min(vertring2(:,2)) + max(vertring2(:,2)))/2 (min(vertring2(:,3)) + max(vertring2(:,3)))/2];

v3 = unit(mpt_rot - bpt_rot);

if plot_figs == 1
    figure;
    plot3(rot_pts(:,1), rot_pts(:,2), rot_pts(:,3), '.', 'MarkerSize', 0.3); hold on
    plot3(notub(:,1), notub(:,2), notub(:,3), '.', 'MarkerSize', .3);
    plot3([mpt_rot(1) bpt_rot(1)], [mpt_rot(2) bpt_rot(2)], [mpt_rot(3) bpt_rot(3)], '- r');
    plot3(mpt_rot(1), mpt_rot(2), mpt_rot(3), '. g', 'MarkerSize', 30);
    plot3(bpt_rot(1), bpt_rot(2), bpt_rot(3), '. r', 'MarkerSize', 30);
    plot3(tpt_rot(1), tpt_rot(2), tpt_rot(3), '. y', 'MarkerSize', 30);
    plot3(vertring(:,1), vertring(:,2), vertring(:,3), '. k');
    plot3(vertring2(:,1), vertring2(:,2), vertring2(:,3), '. k');
end

% fprintf ('Plateau normal: %f, %f, %f | Post. trace ang: %f | Tub pt: %f, %f, %f | SI lower pt: %f, %f, %f | SI mid pt: %f, %f, %f \n', n(1), n(2), n(3), ang, tub_tpt(1), tub_tpt(2), tub_tpt(3), bpt_rot(1), bpt_rot(2), bpt_rot(3), mpt_rot(1), mpt_rot(2), mpt_rot(3));

%% Measure length of tibia
tib_length = abs(tpt_rot(3) - mean(notub(notub(:,3) < min(notub(:,3)) + 5,3)));

%% Find remaining axes
v2 = unit(cross(v3,[1,0,0]));
% front = findAnteriorTibia(v2,strrep(fld,'Geometries\Tib.iv',''));
% if front == 0 
%     v2 = -v2;
% end

if v2(2) < 0
    v2 = -v2;
end

v1 = unit(cross(v2,v3));
v3 = unit(cross(v1,v2));
ninertial = [v3' v1' v2'];

% figure;
% plot3(rot_pts(:,1), rot_pts(:,2), rot_pts(:,3), '.', 'MarkerSize', 0.3); hold on
% plotcs_colour([ninertial [0;0;0]; 0 0 0 1],50,[1 0 0;0 1 0;0 0 1]);

%% Find spines to place origin
% most proximal point - should be medial spine peak
[~,A] = max(med_pts(:,3));
pk1 = med_pts(A,:);

% finding lateral spine peak
pk_range = 10; % mm
[idx,~] = find(((pk1(1) - med_pts(:,1)).^2 + (pk1(2) - med_pts(:,2)).^2) > pk_range^2); % index for points outside a 10mm radius from medial spine
pk_pts = med_pts(idx,:); % possible peak points
[~,B] = max(pk_pts(:,3)); % most proximal point
pk2 = pk_pts(B,:);

% center of peaks
tpt2 = mean([pk1;pk2]);

if plot_figs == 1
    figure;
    plot3(med_pts(:,1), med_pts(:,2), med_pts(:,3), '.', 'MarkerSize', .3); hold on
    plot3(pk1(1), pk1(2), pk1(3), '. r', 'MarkerSize', 20);
    plot3(pk2(1), pk2(2), pk2(3), '. y', 'MarkerSize', 20);
    plot3(tpt2(1), tpt2(2), tpt2(3), '. k', 'MarkerSize', 20);
    legend({'tibia', 'peak 1', 'peak 2', 'origin'});
    title('Spines and CS origin');
end

% cappts = med_pts(med_pts(:,3)>(top-3),:); %get top 3 mm of tibia and find centerpoint 
% tpt2 = mean(cappts);

% plot
% figure;
% plot3(pts(:,1), pts(:,2), pts(:,3), '.', 'MarkerSize', .3); hold on; % tibia
% plot3(cappts(:,1), cappts(:,2), cappts(:,3), '. b'); 
% plot3(tpt2(1), tpt2(2), tpt2(3), '. g');
% plot3(tpt(1), tpt(2), tpt(3), '. r'); % for comparison

% Notch  
notch = [tpt2';1];

% rotate ninertial back to plateau-normal coords
ninertial_rot = ninertial' * Rx';
ninertial = ninertial_rot';

% TIBACS in ellipse-defined coordinates
ell_acs = [ninertial notch(1:3); 0 0 0 1];

% plot
% figure;
% % tibia surface points and CS
% plot3(med_pts(:,1), med_pts(:,2), med_pts(:,3), '.', 'MarkerSize', .3); hold on
% plotcs_colour(med_acs,50,[1 0 0;0 1 0;0 0 1]); % acs

%% Transform med_acs back to MRI coordinates
temp_acs = t * ell_acs; % medial normal to temp local

% plot
% figure;
% % tibia surface points and CS
% plot3(temp_locpts(:,1), temp_locpts(:,2), temp_locpts(:,3), '.', 'MarkerSize', .3); hold on
% plotcs_colour(temp_acs,50,[1 0 0;0 1 0;0 0 1]); % acs

tibacs = [t_med [0;0;0]; 0 0 0 1] * temp_acs; % temp local to MRI

if plot_figs == 1
    % plot final CS
    figure;
    % tibia surface points and CS
    plot3(pts(:,1), pts(:,2), pts(:,3), '.', 'MarkerSize', .3); hold on
    plotcs_colour(tibacs,50,[1 0 0;0 1 0;0 0 1]); % acs
end

fprintf('ACS calculation runtime: %.2f seconds\n', toc);

end
