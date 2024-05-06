function plot_spec_azims(az, freqs, specs, ...
    title_str, ear_id, y_freq_rng, do_normalize_mean, norm_freq_rng, c_mag_rng)

    if nargin < 4; title_str = ''; end
    if nargin < 5; ear_id = 1; end % 1 or 2
    if nargin < 6; y_freq_rng = [freqs(2), 20e3]; end % in Hz
    if nargin < 7; do_normalize_mean = true; end % true or false
    if nargin < 8; norm_freq_rng = y_freq_rng; end % in Hz
    if nargin < 9; c_mag_rng = [-30, 30]; end % in dB

    if all(isnan(specs), 'all')
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

    % if ~is_az
    %     % If it wasn't done before the rendering:
    %     % Translate SH spherical coordinate convention into continous elevation
    %     % rotation (azimuth 0 deg, elevation 0deg to -359 deg)
    %     idx = az == 180 | az == -180;
    %     el(idx) = -el(idx) - 180;
    %     az(idx) = 0;
    %     idx = el > 0;
    %     el(idx) = el(idx) - 360;
    %     clear idx;
    % 
    %     el = sph2nav(el);
    % end

    % rearange data so there is no plotting artefacts
    if is_az
        [az, idx] = sort(az);
    else
        [el, idx] = sort(el);
    end
    specs = specs(:, idx, :);

    spec = db(abs(specs(:, :, ear_id)));
    if do_normalize_mean
        norm_freqs = get_freqs_inclusive(freqs, norm_freq_rng);
        norm_lvl = mean(spec(norm_freqs, :), 'all');
        spec = spec - norm_lvl;
    else
        norm_lvl = 0;
    end
    spec(isinf(-spec)) = -300; % replace -Inf

    % repeat last element for it to be plotted and shift data labels by half 
    % a step size for them to be in the center of the surface bar
    spec(:, end+1) = spec(:, end);
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
        surface(az, freqs, spec, 'EdgeColor', 'None', 'FaceColor', 'Interp');
        set(ax, 'YScale', 'Log'); % in case of yaw
    else
        surface(freqs, el, spec.', 'EdgeColor', 'None', 'FaceColor', 'Interp');
        set(ax, 'XScale', 'Log');
    end
    set(ax, 'Layer', 'Top');

    if ~isempty(title_str{1})
        if ear_id == 1
            title_str{end} = [title_str{end}, '   (left ear)'];
        elseif ear_id == 2
            title_str{end} = [title_str{end}, '   (right ear)'];
        end
        title(title_str, 'Interpreter', 'None');
    end
    
    [ticks, ticklabels] = get_freqs_ticks();
    if is_az
        xticks(-360 : 30 : 360);
        ax.XAxis.MinorTick = 'on';
        ax.XAxis.MinorTickValues = -360 : 5 : 360;
        yticks(ticks);
        yticklabels(ticklabels);
    else
        yticks(-360 : 30 : 360);
        ax.YAxis.MinorTick = 'on';
        ax.YAxis.MinorTickValues = -360 : 5 : 360;
        xticks(ticks);
        xticklabels(ticklabels);
    end

    grid off;
    box on;
    axis tight;
    if is_az
        ylim(y_freq_rng);
    else
        xlim(y_freq_rng);
    end
    caxis(ax, c_mag_rng);

    if is_az
        xlabel('Head yaw rotation in degree');
        ylabel('Frequency in Hz');
    else
        xlabel('Frequency in Hz');
        ylabel('Head pitch rotation in degree');
    end

    cb_str = 'Magnitude in dB';
    if norm_lvl
        cb_str = sprintf('%s   (normalized by %+.1f dB)', cb_str, -norm_lvl);
    end
    cb = colorbar('EastOutside');
    cb.Label.String = cb_str;
end
