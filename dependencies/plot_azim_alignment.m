function plot_azim_alignment(az, min_shift, min_shift_ind, ...
    specs_diff, specs_diff_mean, fs, title_str, diff_freq_rng, ...
    do_normalize)

    if nargin < 8; do_normalize = true; end

    if ~iscell(title_str); title_str = {title_str}; end

    tl = tiledlayout(2, size(specs_diff, 3), 'TileSpacing', 'tight', 'Padding', 'tight');
    drawnow;
    colormap(AKcolormaps('Greys', 128, true));

    if ~isvector(az)
        % Split off azimuth if elevation data is given as well
%         el = az(2, :);
        az = az(1, :);
    end

    % rearange data so there is no plotting artefacts
    [az, az_idx] = sort(az);
    specs_diff = specs_diff(:, az_idx, :);
    specs_diff_mean = specs_diff_mean(az_idx, :);

    if do_normalize
        norm_diff_sum = min(specs_diff_mean, [], 'all');
        specs_diff_mean = specs_diff_mean / norm_diff_sum;
        norm_diff = min(specs_diff, [], [1, 2]);
        specs_diff = specs_diff ./ norm_diff;
        norm_diff = squeeze(norm_diff);
    end

    specs_diff_mean = db(abs(specs_diff_mean));
    specs_diff_lim = db(abs(specs_diff));
    specs_diff_lim(isinf(-specs_diff_lim)) = nan; % replace -Inf
%     specs_diff_lim = [ ...
%         ceil(min(specs_diff_lim, [], 'all', 'omitnan') / 5) * 5, ...
%         ceil(max(specs_diff_lim, [], 'all', 'omitnan') / 5) * 5];
    specs_diff_lim = [ ...
        ceil(min(specs_diff_lim, [], 'all', 'omitnan')), ...
        floor(max(specs_diff_lim, [], 'all', 'omitnan'))];
    freqs = linspace(0, fs/2, size(specs_diff, 1));
    [ticks, ticklabels] = get_freqs_ticks(diff_freq_rng);

    for e = 1 : size(specs_diff, 3) % each ear
        ax = nexttile(tl);
        plot_spec_azims(az, freqs, specs_diff, title_str, e, diff_freq_rng, false);
        set(ax, 'XDir', 'Normal', 'Layer', 'Top', ...
            'TickDir', 'Both', 'TickLength', [0.01, 0]);
        ax.YTick = ticks;
        ax.YTickLabel = ticklabels;
        xlabel(ax, 'Azimuth shift in degree');
        cb = get(ax, 'ColorBar');
        if do_normalize
            cb.Label.String = sprintf('%s   (normalized by %+.1f dB)', ...
                cb.Label.String, -db(norm_diff(e)));
        end
        if e == 1
            cb.Location = 'WestOutside';
        else
            yticklabels(ax, ''); ylabel(ax, '');
        end
        caxis(ax, specs_diff_lim);
        grid off;
        drawnow;
    end

    nexttile(tl, [1, 2]);
    plot(az, specs_diff_mean, 'LineWidth', 2.5, 'LineStyle', ':');
    hold on;
    plot(az, sum(specs_diff_mean, 2), 'k', 'LineWidth', 2);
    xline(min_shift_ind, 'LineWidth', 2.5, 'LineStyle', ':');
    xline(min_shift, 'LineWidth', 2);
    ax = gca;
    xticks(-180 : 30 : 180);
    ax.XAxis.MinorTick = 'on';
    ax.XAxis.MinorTickValues = -180 : 5 : 180;

    grid on;
    box on;
    axis tight;
    ylim([floor(min(specs_diff_mean, [], 'all')), ...
        ceil(max(specs_diff_mean, [], 'all'))]);

    xlabel('Azimuth shift in degree');
    label_str = 'Magnitude in dB';
    if do_normalize
        label_str = sprintf('%s   (normalized by %+.1f dB)', ...
            label_str, -db(norm_diff_sum));
    end
    ylabel(label_str);
    title_str{end} = [title_str{end}, '   and weighted-averaged over frequency'];
    title(title_str);

    lgd_str = {sprintf('Left ear   (minimum at %.0f deg)', min_shift_ind(1)), ...
        sprintf('Right ear   (minimum at %.0f deg)', min_shift_ind(2)), ...
        sprintf('Both ears   (minimum at %.0f deg)', min_shift)};
    legend(lgd_str, 'Location', 'Best');
    drawnow;
end
