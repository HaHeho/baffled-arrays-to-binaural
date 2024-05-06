function plot_ETC_azims_IRs(az, times, irs, ...
    title_str, ear_id, do_normalize_peak)

    if nargin < 4; title_str = ''; end
    if nargin < 5; ear_id = 1; end % 1 or 2
    if nargin < 6; do_normalize_peak = true; end

    ETC_DR = [-80, 0]; % in dB
    
    if all(isnan(irs), 'all')
        axis off;
        return;
    end

    if ~iscell(title_str); title_str = {title_str}; end

    if isvector(az)
        is_az = true;
    else
        % Split off azimuth if elevation data is given as well
        el = az(2, :);
%         az = az(1, :);
        az = -az(1, :); % transform azimuth into yaw
        is_az = all(el == 0 | el == 180 | el == -180);
    end

    % rearange data so there is no plotting artefacts
    if is_az
        [az, idx] = sort(az);
    else
        [el, idx] = sort(el);
    end
    irs = irs(:, idx, :);

    etc = db(abs(irs(:, :, ear_id)));
    if do_normalize_peak
        norm_lvl = max(etc, [], 'all');
        etc = etc - norm_lvl;
    else
        norm_lvl = 0;
    end
    etc(isinf(etc)) = -300; % replace -Inf
    
    % repeat last element for it to be plotted and shift data labels by half 
    % a step size for them to be in the center of the surface bar
    etc(:, end+1) = etc(:, end);
    if is_az
        md = mean(diff(az));
        az(end+1) = az(end) + md;
        az = az - md/2; % maybe this should be 0.5 ?
    else
        md = mean(diff(el));
        el(end+1) = el(end) + md;
        el = el - md/2; % maybe this should be 0.5 ?
    end

    ax = gca;
    if is_az
        % % this shows some strange artefacts when using interpolation i.e.,
        % % values below the set color range of the plot become visible
        % surface(az, times, etc, 'EdgeColor', 'None', 'FaceColor', 'Interp');
        surface(az, times, etc, 'EdgeColor', 'None', 'FaceColor', 'Flat');
        set(ax, 'YDir', 'Reverse'); % in case of yaw
    else
        surface(times, el, etc.', 'EdgeColor', 'None', 'FaceColor', 'Flat');
    end
    set(ax, 'Layer', 'Top', 'TickDir', 'Both', 'TickLength', [0.01, 0]);

    if ~isempty(title_str{1})
        if ear_id == 1
            title_str{end} = [title_str{end}, '   (left ear)'];
        elseif ear_id == 2
            title_str{end} = [title_str{end}, '   (right ear)'];
        end
        title(title_str, 'Interpreter', 'None');
    end
    
    if is_az
        xticks(-360 : 30 : 360);
        ax.XAxis.MinorTick = 'on';
        ax.XAxis.MinorTickValues = -360 : 5 : 360;
    else
        yticks(-360 : 30 : 360);
        ax.YAxis.MinorTick = 'on';
        ax.YAxis.MinorTickValues = -360 : 5 : 360;
    end
    
    grid off;
    box on;
    axis tight;
    caxis(ax, ETC_DR);
    
    if is_az
        xlabel('Head yaw rotation in degree');
        ylabel('Time in ms');
    else
        xlabel('Time in ms');
        ylabel('Head pitch rotation in degree');
    end

    cb_str = 'Magnitude in dB';
    if norm_lvl
        cb_str = sprintf('%s   (normalized by %+.1f dB)', cb_str, -norm_lvl);
    end
    cb = colorbar('EastOutside');
    cb.Label.String = cb_str;
end
