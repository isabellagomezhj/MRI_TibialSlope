%% Determine major axis of ellipse fit to axial slice, make transformation matrix, and measure ap depth

function [fin_sem_ax, t, antdist, v1, tpt] = ellipse_ax(pts, slice, rim_z)
% z_min = min(pts(:,3));

if rim_z == 100
    z_max = max(pts(:,3));
    z_min = z_max - 20;
else
    z_max = rim_z;
    z_min = rim_z - 5;
end

max_ell = 0;

% figure;
% plot3(pts(:,1), pts(:,2), pts(:,3), '.', 'MarkerSize', .3); hold on; % tibia

for z_start = z_max:-slice:z_min
    z_end = z_start - slice;
    
    band_indices = (pts(:,3) <= z_start) & (pts(:,3) > z_end);
    % plot3(pts(band_indices,1), pts(band_indices,2), pts(band_indices,3), '.');
    
    if sum(band_indices) > 5
        band_x = pts(band_indices,1);
        band_y = pts(band_indices,2);
    
        % fit ellipse to band
        ellipse = fit_ellipse(band_x, band_y);
    
        if ~isempty(ellipse)
            sem_ax1 = ellipse.long_axis;
            sem_ax2 = ellipse.short_axis;
        
            if sem_ax1+sem_ax2 > max_ell
        
                max_ell = sem_ax1+sem_ax2;
                max_ell_idx = band_indices;

                x0 = ellipse.X0;
                y0 = ellipse.Y0;
                x0_in = ellipse.X0_in;
                y0_in = ellipse.Y0_in;
                sa1 = ellipse.a;
                sa2 = ellipse.b;
                phi = ellipse.phi;
                fin_sem_ax = ellipse.long_axis;

                % AP distance
                antdist = ellipse.short_axis;
            end
        end
    end
end

% top point (ellipse center)
tpt = [x0_in, y0_in, mean(pts(max_ell_idx,3))];

% % plotting the largest ellipse
% figure;
% plot3(pts(:,1), pts(:,2), pts(:,3), '.', 'MarkerSize', .3); hold on; % tibia
% plot3(pts(max_ell_idx,1), pts(max_ell_idx,2), pts(max_ell_idx,3), '.g'); % pts used for fit

R = [cos(phi), sin(phi); -sin(phi), cos(phi)];

% theta = linspace(0, 2*pi);
% X = x0 + sa1*cos(theta);
% Y = y0 + sa2*sin(theta);
% ellipse = R * [X; Y];
% xe = ellipse(1,:);
% ye = ellipse(2,:);
% ze = ones(length(theta)) * mean(pts(max_ell_idx,3));
% plot3(xe, ye, ze, 'r');

% ML axis
ml_ax = R * [x0 + sa1*[-1 1]; [y0 y0]];
vert1 = [ml_ax(1,1), ml_ax(2,1), mean(pts(max_ell_idx,3))];
vert2 = [ml_ax(1,2), ml_ax(2,2), mean(pts(max_ell_idx,3))];
v1 = unit(vert1-vert2); % axis unit vector

% AP axis
ap_ax = R * [[x0 x0]; y0 + sa2*[-1 1]];
vert3 = [ap_ax(1,1), ap_ax(2,1), mean(pts(max_ell_idx,3))];
vert4 = [ap_ax(1,2), ap_ax(2,2), mean(pts(max_ell_idx,3))];
v2 = unit(vert3-vert4); % axis unit vector

% temporary axes
v3_temp = unit(cross(v1,v2));
t = [v1' v2' v3_temp' tpt'; 0 0 0 1];

end

