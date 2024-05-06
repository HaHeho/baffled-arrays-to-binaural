function plot_spec_azims_percentiles(freqs, specs, percentiles, ...
    title_str, ear_id, y_freq_rng, do_normalize_mean, norm_freq_rng, color)

    if nargin < 4; title_str = ''; end
    if nargin < 5; ear_id = 1; end % 1 or 2
    if nargin < 6; y_freq_rng = [freqs(2), 20e3]; end % in Hz
    if nargin < 7; do_normalize_mean = true; end % true or false
    if nargin < 8; norm_freq_rng = y_freq_rng; end % in Hz
    if nargin < 9; color = [0, 0, 0]; end
    
    SPEC_DR = [-20, 20]; % in dB

    if ~iscell(title_str); title_str = {title_str}; end

    if ear_id
        specs = db(abs(specs(:, :, ear_id)));
    else
        specs = db(abs(specs));
    end
    if do_normalize_mean
        norm_freqs = get_freqs_inclusive(freqs, norm_freq_rng);
        norm_lvl = mean(specs(norm_freqs, 1), 'all'); % based on median
        specs = specs - norm_lvl;
    else
        norm_lvl = 0;
    end

    specs_lower = [freqs.', specs(:, 2)];
    specs_upper = [freqs.', specs(:, end)];
    % ignore 0 Hz bin since log-scale does not work otherwise
    specs_comb = [specs_upper(2:end, :); flipud(specs_lower(2:end, :))];

    plot(freqs, specs(:, 1), 'LineWidth', 2, 'Color', color);
    hold on;
    fill(specs_comb(:, 1), specs_comb(:, 2), 0, ...
        'EdgeColor', color, 'FaceColor', color, 'FaceAlpha', 0.3);
    set(gca, 'XScale', 'log');

    if ~isempty(title_str{1})
        if ear_id == 1
            title_str{end} = [title_str{end}, '   (left ear)'];
        elseif ear_id == 2
            title_str{end} = [title_str{end}, '   (right ear)'];
        end
        title(title_str, 'Interpreter', 'None');
    end

    [ticks, ticklabels] = get_freqs_ticks();
    xticks(ticks);
    xticklabels(ticklabels);
    
    grid on;
    box on;
    set(gca, 'Layer', 'Top');

    xlim(y_freq_rng);
    ylim(SPEC_DR);

    y_lbl = 'Magnitude in dB';
    if norm_lvl
        y_lbl = sprintf('%s   (normalized by %+.1f dB)', y_lbl, -norm_lvl);
    end
    xlabel('Frequency in Hz');
    ylabel(y_lbl);
    legend({'Median', sprintf('%0.f^{th} to %0.f^{th} percentile', ...
        percentiles(1), percentiles(end))}, 'Location', 'NorthWest');
end
