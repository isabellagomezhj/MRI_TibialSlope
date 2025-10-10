%% Using posterior trace to find angle of shaft wrt plateau normal

function [ang, medn_ax, antdist, tpt] = post_trace_angle(med_pts, slice, scale, rim_z, plot_figs)
    % find most contoured ellipse to set top point
    [medn_ax, ~, antdist, ~, tpt] = ellipse_ax(med_pts, slice, rim_z);

    % calculate tibia length to set absolute distal boundary
    tib_length = abs(tpt(3) - mean(med_pts(med_pts(:,3) < min(med_pts(:,3)) + 5,3)));
    abs_z_bound = tpt(3) - tib_length;
    
    % set center point at z = antdist down shaft
    vertpt = tpt(:,3) - scale*antdist;
    vertring_idx = (med_pts(:,3) > (vertpt - slice / 2)) & (med_pts(:,3) < (vertpt + slice / 2));
    % vertring_idx = (med_pts(:,3) > (vertpt - slice)) & (med_pts(:,3) < (vertpt + slice));
    vertring = med_pts(vertring_idx, :);
    
    ring_ctr = [mean(vertring(:,1)) mean(vertring(:,2)) mean(vertring(:,3))];
    
    % set top and bottom point of range
    % shaft_ap_range = range(med_pts(med_pts(:,3) < ring_ctr(3), 2)); % AP depth of the bottom points of the shaft
    % z_adj = shaft_ap_range * sin(deg2rad(15)); % very high slope knee (15d) will need this extra distance up/down shaft to find smallest ring that is same as 0d knee
    
    % shaft_ap_range = range(vertring(:,2));
    % z_adj = shaft_ap_range * 0.2;

    % z_adj = antdist * 0.1;
    
    % z_adj_min = 2.5;
    % z_adj = 8;
    % band_w = z_adj_min + z_adj;

    z_adj_min = antdist * 0.05;
    z_adj = antdist * 0.15;

    % min_range = ring_ctr(3) - z_adj;
    min_range = ring_ctr(3) - z_adj_min;
    max_range = ring_ctr(3) + z_adj;
    
    % set of points to look for ring
    pos_ring_pts1 = med_pts(med_pts(:,3) > min_range & med_pts(:,3) < max_range, :);
    
    if plot_figs == 1
        figure;
        plot3(med_pts(:,1), med_pts(:,2), med_pts(:,3), '.', 'MarkerSize', 0.3); hold on
        plot3(pos_ring_pts1(:,1), pos_ring_pts1(:,2), pos_ring_pts1(:,3), '.', 'MarkerSize', 0.3);
        plot3(tpt(1), tpt(2), tpt(3), 'r .', 'MarkerSize', 10);
        plot3(ring_ctr(1), ring_ctr(2), ring_ctr(3), 'r .', 'MarkerSize', 10);
        xlabel('X'); ylabel('Y'); zlabel('Z');
    end

    % project onto sagittal plane and find angle of trace w/ z-axis (normal to
    % pleateau)

    pos_ring_pts_sag = [pos_ring_pts1(pos_ring_pts1(:,2) < ring_ctr(2), 2) pos_ring_pts1(pos_ring_pts1(:,2) < ring_ctr(2), 3)];
    
    % not actually necessary to find all of the trace points but good for now
    z_min = min(pos_ring_pts_sag(:,2));
    z_max = max(pos_ring_pts_sag(:,2));
    
    tracepts = [];
    ct = 0;

    for z_start = z_max:-1:z_min
        z_end = z_start - 1;

        band_indices = (pos_ring_pts_sag(:,2) <= z_start) & (pos_ring_pts_sag(:,2) > z_end);

        % select points in the most posterior 0.25 mm of the band
        band_pts = pos_ring_pts_sag(band_indices,:); % all points in the 5mm band
        post_val = min(band_pts(:,1)) + 0.25;
        post_pts = band_pts(band_pts(:,1) < post_val, :);

        if height(post_pts) > 0
            % pt = mean(post_pts,1);
            % 
            % if ct == 0 || abs(pt(1) - tracepts(ct,1)) < 0.2
            %     tracepts = [tracepts; mean(post_pts,1)];
            %     ct = ct + 1;
            % end    
            tracepts = [tracepts; post_pts];
        end

        % % if there's more than 10 points in the band, use the most posterior
        % % ones
        % if sum(band_indices) > 15
        %     band_pts = pos_ring_pts_sag(band_indices,:); % all points in the 5mm band
        %     post_pts = sortrows(band_pts,1); % anterior to posterior
        % 
        %     tracepts = [tracepts; post_pts(1:10,:)];
        % end
    end

    % % most distal and posterior point
    % dp_pt = tracepts(tracepts(:,2) == min(tracepts(:,2)), :);
    % dp_pt = tracepts(tracepts(:,2) < min(tracepts(:,2)) + 5, :);
    % 
    % if height(dp_pt) > 1
    %     dp_pt = mean(dp_pt);
    % end
    % 
    % % most proximal and posterior point
    % % pp_pt = tracepts(tracepts(:,2) == max(tracepts(:,2)), :);
    % pp_pt = tracepts(tracepts(:,2) > max(tracepts(:,2)) - 5, :);
    % 
    % if height(pp_pt) > 1
    %     pp_pt = mean(pp_pt);
    % end
    % 
    % if plot_figs == 1
    %     figure;
    %     plot(pos_ring_pts_sag(:,1), pos_ring_pts_sag(:,2), '.'); hold on
    %     plot(tracepts(:,1), tracepts(:,2), 'r.');
    %     plot(dp_pt(:,1), dp_pt(:,2), 'g.', 'MarkerSize', 20);
    %     plot(pp_pt(:,1), pp_pt(:,2), 'k.', 'MarkerSize', 20);
    %     xlabel('Y'); ylabel('Z');
    % end
    % 
    % ang1 = atan2d(pp_pt(2)-dp_pt(2), pp_pt(1)-dp_pt(1)) - 90;

    % linear regression of posterior points
    % reg_tracepts = tracepts(tracepts(:,2) < max(tracepts(:,2)) - 1 & tracepts(:,2) > min(tracepts(:,2)) + 1, :);
    % reg_tracepts = tracepts;
    % reg_tracepts = tracepts(tracepts(:,1) < prctile(tracepts(:,1), 85), :);
    % reg_tracepts = tracepts(tracepts(:,1) < min(tracepts(:,1)) + band_w*tand(20), :);
    % reg_tracepts = [movmean(tracepts(:,1), 30) tracepts(:,2)];
    % reg_tracepts = tracepts(tracepts(:,2) > max(tracepts(:,2)) - band_w & tracepts(:,1) < min(tracepts(:,1)) + band_w*tand(15), :);
    reg_tracepts = tracepts(tracepts(:,2) > abs_z_bound & tracepts(:,1) < prctile(tracepts(:,1),90), :);

    mdl = fitlm(reg_tracepts(:,1), reg_tracepts(:,2));

    ang1 = (abs(mdl.Coefficients.Estimate(2)) / mdl.Coefficients.Estimate(2)) * (atan2d(abs(mdl.Coefficients.Estimate(2)), 1) - 90);

    % check if points aren't linear
    if mdl.Rsquared.Ordinary < 0.2 || abs(ang1) > 20
        if abs(ang1) > 40
            % posterior trace is probably close to vertical (slope of 0)
            ang1 = 0;
        else
            % posterior trace is probably near the bottom of the tibia and
            % curved - need to find edge of the bone
            reg_tracepts_diffs = diff(reg_tracepts(:,2));
            [~, edge_idx] = min(reg_tracepts_diffs);
            edge_z = reg_tracepts(edge_idx+1,2);
            reg_tracepts = reg_tracepts(reg_tracepts(:,2) > edge_z, :);
            mdl = fitlm(reg_tracepts(:,1), reg_tracepts(:,2));
            ang1 = (abs(mdl.Coefficients.Estimate(2)) / mdl.Coefficients.Estimate(2)) * (atan2d(abs(mdl.Coefficients.Estimate(2)), 1) - 90);
        end
    end

    if abs(ang1) > 30
        % posterior trace is probably close to vertical (slope of 0)
        ang1 = 0;
    end

    % fprintf('Ang1: %s \n', ang1);

    if plot_figs == 1
        figure;
        plot(mdl);
        title('Posterior trace regression #1');
    end

    % repeat posterior trace fitting, accounting for preliminary angle
    % found
    corr = antdist * (1 - cosd(ang1));
    % if ang1 < 1
    %     corr = -corr;
    % end

    % set center point at z = antdist down shaft
    vertpt = tpt(:,3) - scale*antdist + corr;
    vertring_idx = (med_pts(:,3) > (vertpt - slice / 2)) & (med_pts(:,3) < (vertpt + slice / 2));
    % vertring_idx = (med_pts(:,3) > (vertpt - slice)) & (med_pts(:,3) < (vertpt + slice));
    vertring = med_pts(vertring_idx, :);

    ring_ctr = [mean(vertring(:,1)) mean(vertring(:,2)) mean(vertring(:,3))];

    % set top and bottom point of range
    % shaft_ap_range = range(vertring(:,2));
    % z_adj = shaft_ap_range * 0.2;

    % min_range = ring_ctr(3) - z_adj;
    min_range = ring_ctr(3) - z_adj_min;
    max_range = ring_ctr(3) + z_adj;

    % set of points to look for ring
    pos_ring_pts = med_pts(med_pts(:,3) > min_range & med_pts(:,3) < max_range, :);

    if plot_figs == 1
        figure;
        plot3(med_pts(:,1), med_pts(:,2), med_pts(:,3), '.', 'MarkerSize', 0.3); hold on
        plot3(pos_ring_pts(:,1), pos_ring_pts(:,2), pos_ring_pts(:,3), '.');
        plot3(pos_ring_pts1(:,1), pos_ring_pts1(:,2), pos_ring_pts1(:,3), '.');
        plot3(tpt(1), tpt(2), tpt(3), 'r .', 'MarkerSize', 10);
        plot3(ring_ctr(1), ring_ctr(2), ring_ctr(3), 'r .', 'MarkerSize', 10);
        xlabel('X'); ylabel('Y'); zlabel('Z');
        legend({'tibia', 'rd 2', 'rd 1'});
        title('Posterior Trace #2');
    end

    % project onto sagittal plane and find angle of trace w/ z-axis (normal to
    % pleateau)

    pos_ring_pts_sag = [pos_ring_pts(pos_ring_pts(:,2) < ring_ctr(2), 2) pos_ring_pts(pos_ring_pts(:,2) < ring_ctr(2), 3)];

    % not actually necessary to find all of the trace points but good for now
    z_min = min(pos_ring_pts_sag(:,2));
    z_max = max(pos_ring_pts_sag(:,2));

    tracepts = [];
    ct = 0;

    for z_start = z_max:-1:z_min
        z_end = z_start - 1;

        band_indices = (pos_ring_pts_sag(:,2) <= z_start) & (pos_ring_pts_sag(:,2) > z_end);

        % select points in the most posterior 0.25 mm of the band
        band_pts = pos_ring_pts_sag(band_indices,:); % all points in the 5mm band
        post_val = min(band_pts(:,1)) + 0.25;
        post_pts = band_pts(band_pts(:,1) < post_val, :);

        if height(post_pts) > 0
            % pt = mean(post_pts,1);
            % 
            % if ct == 0 || abs(pt(1) - tracepts(ct,1)) < 0.2
            %     tracepts = [tracepts; mean(post_pts,1)];
            %     ct = ct + 1;
            % end    
            tracepts = [tracepts; post_pts];
        end

        % % if there's more than 10 points in the band, use the most posterior
        % % ones
        % if sum(band_indices) > 15
        %     band_pts = pos_ring_pts_sag(band_indices,:); % all points in the 5mm band
        %     post_pts = sortrows(band_pts,1); % anterior to posterior
        % 
        %     tracepts = [tracepts; post_pts(1:10,:)];
        % end
    end

    % most distal and posterior point
    % dp_pt = tracepts(tracepts(:,2) == min(tracepts(:,2)), :);
    % dp_pt = tracepts(tracepts(:,2) < min(tracepts(:,2)) + 5, :);
    % 
    % if height(dp_pt) > 1
    %     dp_pt = mean(dp_pt);
    % end
    % 
    % % most proximal and posterior point
    % % pp_pt = tracepts(tracepts(:,2) == max(tracepts(:,2)), :);
    % pp_pt = tracepts(tracepts(:,2) > max(tracepts(:,2)) - 5, :);
    % 
    % if height(pp_pt) > 1
    %     pp_pt = mean(pp_pt);
    % end
    % 
    % if plot_figs == 1
    %     figure;
    %     plot(pos_ring_pts_sag(:,1), pos_ring_pts_sag(:,2), '.'); hold on
    %     plot(tracepts(:,1), tracepts(:,2), 'r.');
    %     plot(dp_pt(:,1), dp_pt(:,2), 'g.', 'MarkerSize', 20);
    %     plot(pp_pt(:,1), pp_pt(:,2), 'k.', 'MarkerSize', 20);
    %     xlabel('Y'); ylabel('Z');
    % end
    % 
    % ang = atan2d(pp_pt(2)-dp_pt(2), pp_pt(1)-dp_pt(1)) - 90;

    % linear regression of posterior points - rd2
    % reg_tracepts = tracepts(tracepts(:,2) < max(tracepts(:,2)) - 1 & tracepts(:,2) > min(tracepts(:,2)) + 1, :);
    % reg_tracepts = tracepts;
    % reg_tracepts = tracepts(tracepts(:,1) < prctile(tracepts(:,1), 85), :);
    % reg_tracepts = tracepts(tracepts(:,1) < min(tracepts(:,1)) + band_w*tand(20), :);
    % reg_tracepts = tracepts(tracepts(:,2) > max(tracepts(:,2)) - band_w & tracepts(:,1) < min(tracepts(:,1)) + band_w*tand(15), :);
    reg_tracepts = tracepts(tracepts(:,2) > abs_z_bound & tracepts(:,1) < prctile(tracepts(:,1),90), :);

    mdl = fitlm(reg_tracepts(:,1), reg_tracepts(:,2));

    ang = (abs(mdl.Coefficients.Estimate(2)) / mdl.Coefficients.Estimate(2)) * (atan2d(abs(mdl.Coefficients.Estimate(2)), 1) - 90);

    % check if points aren't linear
    if mdl.Rsquared.Ordinary < 0.2 || abs(ang) > 20
        if abs(ang) > 40 || height(reg_tracepts) < 2
            ang = ang1;
        else
            % posterior trace is probably near the bottom of the tibia and
            % curved - need to find edge of the bone
            reg_tracepts_diffs = diff(reg_tracepts(:,2));
            [~, edge_idx] = min(reg_tracepts_diffs);
            edge_z = reg_tracepts(edge_idx+1,2);
            reg_tracepts = reg_tracepts(reg_tracepts(:,2) > edge_z, :);
            mdl = fitlm(reg_tracepts(:,1), reg_tracepts(:,2));
            ang = (abs(mdl.Coefficients.Estimate(2)) / mdl.Coefficients.Estimate(2)) * (atan2d(abs(mdl.Coefficients.Estimate(2)), 1) - 90);
        end
    end

    % figure;
    % plot3(med_pts(:,1), med_pts(:,2), med_pts(:,3), '.', 'MarkerSize', 0.3); hold on
    % plot3(mean(pos_ring_pts(:,1)) * ones(length(tracepts)), tracepts(:,1), tracepts(:,2), '. r', 'MarkerSize', 20);
    % % plot3(tpt(:,1), tpt(:,2), tpt(:,3), '. y', 'MarkerSize', 30);
    % title('Plateau-normal-aligned points and posterior trace');

    if abs(ang) > 30
        % posterior trace is probably close to vertical (slope of 0)
        ang = ang1;
    end

    % fprintf('Ang: %s \n', ang);

    if plot_figs == 1
        figure;
        plot(mdl);
        title('Posterior trace regression #2');
    end
end