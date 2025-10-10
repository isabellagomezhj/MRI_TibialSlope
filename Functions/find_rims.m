function [med_pt, lat_pt, n, rim_z, med_n, lat_n, idxs] = find_rims(tib, kneeside, long_axis)
%% Finding plateau rims
% Adapted from slopeMedPlat.m and slopeLatPlat.m
% "S:\BiomechanicsResearch\groupImhauser\TibialSlopesKneeKinematics\Code\Consolidate"
% Inputs:
% tib = tibia mesh points in local coordinates established by ellipse fit (Amirtharaj method)
% kneeside = laterality
% long_axis = length of major axis of ellipse fit to most contoured slice

og_tib = tib; % for plotting

% flip right knees
if strcmp(kneeside, 'R')
    tib(:,1) = -tib(:,1);
end

tib_initial = tib;

% find lateral plateau bounds in initial local coordinates
latbound_2 = bounding(tib, long_axis, 'lateral');

% find lateral plateau bounds with rotation algorithm
iterations = 20;
[lat_pts, ~] = iter_rot(iterations, tib_initial, latbound_2, 'lateral', long_axis);

% find medial plateau bounds in initial local coordinates
medbound_2 = bounding(tib, long_axis, 'medial');

% find medial plateau bounds with rotation algorithm
[med_pts, ~] = iter_rot(iterations, tib_initial, medbound_2, 'medial', long_axis);

% flip rim points for right knees
if strcmp(kneeside, 'R')
    lat_pts(:,1) = -lat_pts(:,1);
    med_pts(:,1) = -med_pts(:,1);
end

% best-fit plane to medial pts
% [med_n, med_plane, mp] = affine_fit(med_pts);

% best-fit plane to all rim pts
[n, med_plane, mp] = affine_fit([med_pts; lat_pts]);
rim_z = mp(3);

% return most anterior point on each rim
lat_prct_val = prctile(lat_pts(:,2), 90);
[~, lat_idx] = min(abs(lat_pts(:,2) - lat_prct_val));
lat_pt = lat_pts(lat_idx, :);

med_prct_val = prctile(med_pts(:,2), 90);
[~, med_idx] = min(abs(med_pts(:,2) - med_prct_val));
med_pt = med_pts(med_idx, :);

% return normal to each plateau
[med_n, ~, ~] = affine_fit(med_pts);
[lat_n, ~, ~] = affine_fit(lat_pts);

% return indices for each plateau's points (need to account for floating
% point precision errors)
full_rim = [med_pts; lat_pts];

M = size(full_rim, 1);
idxs = zeros(M, 1);
tol = 1e-6;

for i = 1:M
    dists = sqrt(sum((og_tib - full_rim(i,:)).^2, 2));
    [min_dist, min_loc] = min(dists);

    if min_dist < tol
        idxs(i) = min_loc;
    else
        idxs(i) = 0;
    end
end


% % plotting final set of rim points found
% figure;
% plot3(og_tib(:,1), og_tib(:,2), og_tib(:,3), '.', 'MarkerSize', 0.3); hold on
% plot3(lat_pts(:,1), lat_pts(:,2), lat_pts(:,3), 'r.', 'MarkerSize', 20);
% plot3(med_pts(:,1), med_pts(:,2), med_pts(:,3), 'r.', 'MarkerSize', 20);
% plot3(lat_pt(:,1), lat_pt(:,2), lat_pt(:,3), 'g.', 'MarkerSize', 20);
% plot3(med_pt(:,1), med_pt(:,2), med_pt(:,3), 'g.', 'MarkerSize', 20);
% xlabel('X'); ylabel('Y'); zlabel('Z');
% title('Final Rim Points');

% plotting rim points and plane fit to them w/ normal vector
% figure;
% plot3(og_tib(:,1), og_tib(:,2), og_tib(:,3), '.', 'MarkerSize', 0.3); hold on
% plot3(med_pts(:,1), med_pts(:,2), med_pts(:,3), 'r.', 'MarkerSize', 20);
% plot3(lat_pts(:,1), lat_pts(:,2), lat_pts(:,3), 'r.', 'MarkerSize', 20);
% plot3(lat_pt(:,1), lat_pt(:,2), lat_pt(:,3), 'g.', 'MarkerSize', 20);
% plot3(med_pt(:,1), med_pt(:,2), med_pt(:,3), 'g.', 'MarkerSize', 20);
% plot3(mp(1), mp(2), mp(3), 'b.'); %center of plane
% quiver3(mp(1), mp(2), mp(3), n(1), n(2), n(3), 20, 'b', 'LineWidth', 3); %normal vector
% 
% % plot best fit plane
% [S1,S2] = meshgrid([-30 0 30]);
% 
% X = mp(1)+[S1(:) S2(:)]*med_plane(1,:)';
% Y = mp(2)+[S1(:) S2(:)]*med_plane(2,:)';
% Z = mp(3)+[S1(:) S2(:)]*med_plane(3,:)';
% 
% surf(reshape(X,3,3),reshape(Y,3,3),reshape(Z,3,3),'facecolor','blue','facealpha',0.5);
% xlabel('X'); ylabel('Y'); zlabel('Z');
% % title('Medial Rim and Best-Fit Plane');
% title('Plateau Rims and Best-Fit Plane');

end

%%%% Functions
%% Iterative rotations function
function [bounds, n] = iter_rot(iterations, tib_initial, latbound_2, side, long_axis)
%----------------------------ITERATE ROTATION STEPS------------------------------------------------
for m = 1:iterations

    clear tib_r ell_idx ellipse long_axis_r

    % calculate ap angle based on initial run of kneedle
    [n_1,~,~] = affine_fit(latbound_2); 

    if n_1 == [0; 0; 1]
        tib_r = tib_initial;
    else
        % rotate WHOLE TIB POINT CLOUD SO PLATEAU = XY PLANE

        if n_1(3) < 0
            n_1 = -n_1;
        end
        
        a = [1 0 0];
        y = cross(n_1, a);
        y_norm = y/norm(y);
        x = cross(y_norm, n_1);
        x_norm = x/norm(x);
        r = cat(2, x_norm', y_norm');
        r = cat(2, r, n_1);
        
        tib_r = tib_initial*r;
    end
    
    % fit ellipse to most contoured slice
    ell_idx = (tib_r(:,3) <= 2.5) & (tib_r(:,3) > -2.5);
    ellipse = fit_ellipse(tib_r(ell_idx,1), tib_r(ell_idx,2));

    if isempty(ellipse)
        long_axis_r = long_axis;
    else
        long_axis_r = ellipse.long_axis;
    end

    % find rim points for rotated tibia
    latbound_flat = bounding(tib_r, long_axis_r, side);
  
    % remove outliers
    quantiles = quantile(latbound_flat(:,3), [.25 .5 .75]);
    q1 = quantiles(1);
    q3 = quantiles(3);
    iqr = q3-q1;

    new_latbound_methodA = zeros(length(latbound_flat)-5,3);
    new_count_methodA = 1;

    for x=1:length(latbound_flat)
        if latbound_flat(x,3) < (q3+3*iqr) && latbound_flat(x,3) > (q1-3*iqr)
            new_latbound_methodA(new_count_methodA,1) = latbound_flat(x,1);
            new_latbound_methodA(new_count_methodA,2) = latbound_flat(x,2);
            new_latbound_methodA(new_count_methodA,3) = latbound_flat(x,3);
            new_count_methodA = new_count_methodA + 1;
        end
    end

    latbound_flat = new_latbound_methodA;
            
    latbound_2 = latbound_flat*r';

    % figure;
    % plot3(latbound_2(:,1), latbound_2(:,2), latbound_2(:,3), '. r'); hold on
    % plot3(tib_initial(:,1), tib_initial(:,2), tib_initial(:,3), '.', 'MarkerSize', 0.3);
    
end

% ------------END ITERATIVE ROTATIONS STEPS--------------------------

bounds = latbound_flat*r'; % final version to return
[n,~,~] = affine_fit(bounds);

end

%% Function to find rim points
% Inputs:
% int_tib = intact tibia pts
% long_axis = length of long axis from ellipse fit
% Outputs:
% latbound_2 = points bounding plateau rim

function latbound_2 = bounding(int_tib, long_axis, side)

[tib, PM] = find_maxes(int_tib, side);

% identify the plateau boundary using radial slices originating at approx.
% plateau center
range2 = 0.2; %range used to identify points for a radial slice

if strcmp(side, 'medial')
    % identify rough center of plateau
    center = [long_axis/4 PM(2) 0];

    theta = 0.1; %starting angle (original = -0.1*pi)
    % thetamax = 0.5*pi; %finishing angle for anterior half of rim - ruins
    % plane fitting
    thetamax = 1.0*pi; %finishing angle for full rim (original = 1.1*pi)
else
    center = [-long_axis/4 PM(2) 0];

    theta = 1.0*pi; %starting angle (original = 0.9*pi)
    % thetamax = 1.5*pi; %finishing angle for anterior half of rim
    thetamax = 1.9*pi; %finishing angle for full rim (original = 2.1*pi)
end

rsi = 1; %initiate counter

incr = deg2rad(3); %interval; iterate the radial slices by this amount from 0 to 360
num_rim_points = floor((thetamax-theta)/incr);
latbound_2 = nan(num_rim_points,3);
axrot = [0 0 1]; %axis for rotation

tib2 = tib(:,1:2); 
tib2(:,3) = 0;

% plot plateau pts in xy plane and estimated center
% figure(13);
% plot(tib2(:,1), tib2(:,2), '.'); hold on
% plot(center(:,1), center(:,2), 'r.');
% point = center + [100*sin(theta) 100*cos(theta) 0];
% plot(point(:,1), point(:,2), 'g.');

% fig_ct = 7;

while abs(theta) < abs(thetamax)
    clear radslice radslice2 halfradslice a b g h cutdist point

    % identify a point that, along with center, defines the direction of slice
    point = center + [100*sin(theta) 100*cos(theta) 0];

    % figure(13);
    % plot([center(:,1); point(:,1)], [center(:,2); point(:,2)]); hold on

    %Cycle through each point in tib2 and calculate the perpendicular
    %distance from the slice plane.  Store this distance in the fourth
    %column of tib2
    for col = 1:size(tib2, 1)
        ab = point - center;
        bc = tib2(col,1:3) - center;
        tib2(col,4) = norm(cross(ab,bc)) / norm(ab);        
    end


    %Identify all the points in tib2 that are within "range2" mm of the
    %current slice
    [a, ~] = find(tib2(:,4)<range2 & tib2(:,4)>-range2);
    radslice = tib(a,:);

    %For each slice, only points moving radially outward in a single
    %direction are desired.  Points moving radially outward in the opposite
    %direction are removed.
    radslice2 = radslice;
    radslice2(:,3) = 0;
    cutdist = norm(point-center);

    %Cycle through each point in the current slice.  Find the distance
    %of each point from "point"
    for col2 = 1:size(radslice2,1)
        radslice2(col2,3) = norm(point-radslice(col2,:));
    end

    %Identify points that are farther from "center" than "point" and remove
    %them from the matrix.  This eliminates points in the current slice
    %that are in the opposite direction (i.e. in the pi direction  instead
    %of in the zero direction)
    [g, ~] = find(radslice2(:,3)<cutdist);
    halfradslice = radslice(g,:);

    %Identify the "knee" in the data using a script called
    %find_transition_step_kneedle created by Bobby Kent
    %(bobbykent14@gmail.com).  This script requires that the data be in
    %quadrant I, and y must increase with x.  The data is transformed
    %in this section in order to ensure these requirements are met.
    rotangle = theta - pi;

    hrstrans = halfradslice - repmat(center, size(halfradslice,1), 1); %translate to origin

    try
        v_rot = rodrigues_rot(hrstrans, axrot, rotangle); %rotate so the slice is always at 180 degrees
        hrsorigin = v_rot + repmat(center, size(v_rot,1), 1); %translate back to initial origin

        [Q, ~] = min(hrsorigin); %find minimum values within translated and rotated matrix
        quad1 = hrsorigin - repmat(Q, size(hrsorigin,1),1); %translate to ensure all points are in quadrant 1
    
        [t_step] = find_transition_step_kneedle(quad1(:,3), quad1(:,2)); %identify "knee"
    
        latbound_2(rsi, 1:3) = halfradslice(t_step,:);

        % figure(fig_ct + 200);
        % plot3(quad1(:,1), quad1(:,2), quad1(:,3), '.'); hold on
        % plot3(latbound_2(rsi,1), latbound_2(rsi,2), latbound_2(rsi,3), 'r.');
        % plot(quad1(:,3), quad1(:,2), '.'); hold on
        % plot(quad1(t_step,3), quad1(t_step,2), 'r.');

    catch
        latbound_2(rsi, :) = nan;
    end

    % figure(fig_ct);
    % plot(tib2(:,1), tib2(:,2), '.'); hold on
    % plot(radslice2(:,1), radslice2(:,2), 'g.');
    % plot(halfradslice(:,1), halfradslice(:,2), 'y.');
    % plot(latbound_2(rsi, 1), latbound_2(rsi, 2), 'r.');
    % fig_ct = fig_ct + 1;

    % increment counter variables
    theta = theta + incr;
    rsi = rsi + 1;

end

latbound_2 = latbound_2(~any(isnan(latbound_2), 2), :);

% plot set of rim points found
% figure;
% plot3(tib(:,1), tib(:,2), tib(:,3), '.', 'MarkerSize', 0.3); hold on
% plot3(latbound_2(:,1), latbound_2(:,2), latbound_2(:,3), 'r.');
% xlabel('X'); ylabel('Y'); zlabel('Z');

end

%% Function to find maximum points on a region of the tibia

function [tib, PM] = find_maxes(pts, side)

% find most proximal point of the tibia
% [~,I] = max(pts);
% PM = pts(I(3),:);

% trying to deal with case where the front of the plateau or the rim 
% is higher than the spines bc of MRI position
y_mid = mean(pts(pts(:,3) >= 0, 2));
% y_mid = mean(pts(pts(3,:) >= 0, 2));
% post_pts = pts(pts(:,2) <= y_mid, :);

% constrain to posterior half and middle 30-70% (ML) of points
post_pts = pts(pts(:,2) <= y_mid & pts(:,1) > prctile(pts(pts(:,3) > 0, 1),30) & pts(:,1) < prctile(pts(pts(:,3) > 0, 1),70) & pts(:,3) > 0, :);
% post_pts = pts(1:3, pts(2,:) <= y_mid & pts(1,:) > prctile(pts(pts(3,:) > 0, 1),30) & pts(1,:) < prctile(pts(pts(3,:) > 0, 1),70) & pts(3,:) > 0);
[~, I] = max(post_pts(:,3));
PM = post_pts(I, :);

% figure;
% plot3(pts(:,1), pts(:,2), pts(:,3), '.', 'MarkerSize', 0.3); hold on
% plot3(PM(1), PM(2), PM(3), '.r', 'MarkerSize',15);
% xlabel('X'); ylabel('Y'); zlabel('Z');

% exclude pts >10mm medial and >30mm inferior to peak of spine
cutz = PM(3)-30;

if strcmp(side, 'medial')
    cutx = PM(1)+5;
    tib = pts(pts(:,1)>cutx & pts(:,3)>cutz,:);
else
    cutx = PM(1)-5;
    tib = pts(pts(:,1)<cutx & pts(:,3)>cutz,:);
end

% % plot pts that will be used to find rim and measure slope
% figure;
% plot3(tib(:,1), tib(:,2), tib(:,3), '.', 'MarkerSize', 0.3); hold on
% plot3(pts(:,1), pts(:,2), pts(:,3), '.');
% 
% xv = linspace(min(tib(:,1)), max(tib(:,1)), 100);
% yv = linspace(min(tib(:,2)), max(tib(:,2)), 100);
% [X,Y] = meshgrid(xv, yv);
% Z = griddata(tib(:,1), tib(:,2), tib(:,3), X, Y);
% surf(X, Y, Z);
% title(side);
end


