function plot_spec_compare_percentiles(lgd_str, lgd_title_str, ...
    irs_fro_array, spec_hor_perc_array, percentiles, fs_vec, ...
    apply_frac_sm, do_normalize_mean, norm_freq_rng, plot_freq_rng, plot_mag_dr)

    ETC_DR = 100; % in dB
    MAG_LIM_FREQ_RNG = [40, 20e3]; % in Hz

    if nargin < 6; apply_frac_sm = 24; end % in fractional octaves
    if nargin < 7; do_normalize_mean = true; end % true or false
    if nargin < 8; norm_freq_rng = [20, 1e3]; end % in Hz
    if nargin < 9; plot_freq_rng = [20, 20e3]; end % in Hz
    if nargin < 10; plot_mag_dr = 45; end % in dB

    n_data = length(lgd_str);
    if n_data ~= length(irs_fro_array) || n_data ~= length(spec_hor_perc_array)
        error('Mismatch in provided cell sizes.');
    end
    if length(fs_vec) == 1
        fs_vec = repmat(fs_vec, 1, n_data);
    end

    colors = AKcolormaps('Dark2', n_data - 1);
    colors = [0, 0, 0; colors]; % add black as first color

    for d = 1 : n_data
        etc_fro_array{d} = db(abs(irs_fro_array{d})); %#ok<*AGROW> 

        spec_fro_array{d} = AKboth2singleSidedSpectrum(fft(irs_fro_array{d}));
        if apply_frac_sm
            spec_fro_array{d} = AKfractOctSmooth( ...
                spec_fro_array{d}, 'amp', fs_vec(d), apply_frac_sm);
        end
        spec_fro_array{d} = db(abs(spec_fro_array{d}));
        
        times_array{d} = 1:size(etc_fro_array{d}, 1);
        freqs_array{d} = linspace(0, fs_vec(d)/2, size(spec_fro_array{d}, 1));
        
        if do_normalize_mean
            norm_freqs = get_freqs_inclusive(freqs_array{d}, norm_freq_rng);
            
            norm_lvl(d) = mean(spec_fro_array{d}(norm_freqs, :), 'all');
            spec_fro_array{d} = spec_fro_array{d} - norm_lvl(d);

%             norm_lvl(d) = mean(spec_hor_perc_array{d}(norm_freqs, 1), 'all'); % based on median
%             spec_hor_perc_array{d} = spec_hor_perc_array{d} - norm_lvl(d);
        end
    end
    
%     % adjust shown minimum frequency
%     norm_freq_rng(1) = max([norm_freq_rng(1), min(cellfun(@(x) x(2), freqs_array))]);

    tl = tiledlayout(2, 2, 'TileSpacing', 'tight', 'Padding', 'tight');
    drawnow();

    %%
    ax = nexttile(tl);
    hold on;
    data_x_lim = 0; data_y_lim = -inf;
    for d = 1 : n_data
        h = plot(times_array{d}, etc_fro_array{d}, ...
            'LineWidth', 1, 'LineStyle', '-', 'Color', colors(d, :));
        h(1).DisplayName = lgd_str{d};
        data_x_lim = max([data_x_lim, h(1).XData(end)]);
        data_y_lim = max([data_y_lim, max(h(1).YData, [], 'all')]);
        if size(irs_fro_array{d}, 2) > 1
            h(2).LineStyle = ':';
            h(2).LineWidth = 1.75;
            h(2).HandleVisibility = 'off';
        end
        data_x_lim = max([data_x_lim; arrayfun(@(x) x.XData(end), h)]);
        data_y_lim = max([data_y_lim; arrayfun(@(y) max(y.YData), h)]);
    end
    xlim([0, data_x_lim]);
    data_y_lim = ceil(data_y_lim / 5) * 5;
    ylim([data_y_lim - ETC_DR, data_y_lim]);

    grid on;
    box on;
    set(ax, 'Layer', 'Top');

    xlabel('Time in samples');
    ylabel('Magnitude in dB');
    title('Energy Time Curve', 'Interpreter', 'None');

    lgd = legend('Interpreter', 'None', 'Location', 'Best');
    title(lgd, [lgd_title_str, '  (frontal)']);
    clear ax data_x_lim data_y_lim;
    drawnow();

    %%
    ax = nexttile(tl);
    data_y_lim = -inf;
    data_op = @max;
    for d = 1 : n_data
        h = semilogx(freqs_array{d}, spec_fro_array{d}, ...
            'LineWidth', 2, 'LineStyle', '-', 'Color', colors(d, :));
        h(1).DisplayName = lgd_str{d};
        data_y_lim = data_op([data_y_lim, data_op(h(1).YData, [], 'all')]);
        if size(irs_fro_array{d}, 2) > 1
            h(2).LineStyle = ':';
            h(2).LineWidth = 2.5;
            h(2).HandleVisibility = 'off';
        end
        freq_rng = get_freqs_inclusive(freqs_array{d}, MAG_LIM_FREQ_RNG);
        data_y_lim = data_op([data_y_lim; arrayfun(@(y) data_op(y.YData(freq_rng)), h)]);
        hold on;
    end

%     xlim([norm_freq_rng(1), freqs_array{end}(end)]);
    xlim(plot_freq_rng);
    data_y_lim = ceil(data_y_lim / 5) * 5;
    ylim([data_y_lim - plot_mag_dr, data_y_lim]);

    [ticks, ticklabels] = get_freqs_ticks();
    xticks(ticks);
    xticklabels(ticklabels);
    
    grid on;
    box on;
    set(ax, 'Layer', 'Top', 'YAxisLocation', 'Right');

    xlabel('Frequency in Hz');
    ylabel('Magnitude in dB');
    title_str = 'Magnitude Spectrum';
    if apply_frac_sm
        title_str = sprintf('%s   (1/%d oct. smoothing)', title_str, apply_frac_sm);
    end
    if do_normalize_mean
        title_str = [title_str, '   (normalized)'];
    end
    title(title_str, 'Interpreter', 'None');

    lgd = legend('Interpreter', 'None', 'Location', 'Best');
    title(lgd, [lgd_title_str, '  (frontal)']);
    clear ax data_y_lim data_op;
    drawnow();
    
    %%
    nexttile(tl);
    data_y_lim = -inf;
    data_op = @max;
    for d = 1 : n_data
        hold on;
        plot_spec_azims_percentiles(freqs_array{d}, spec_hor_perc_array{d}, ...
            percentiles, '', 1, plot_freq_rng, ...
            do_normalize_mean, norm_freq_rng, colors(d, :));
        ax = gca;
        freq_rng = get_freqs_inclusive(freqs_array{d}, MAG_LIM_FREQ_RNG);
        data_y_lim = data_op([data_y_lim; data_op(ax.Children(1).YData(freq_rng))]);
        % remove patch component and only keep line
        delete(ax.Children(2));
    end

    data_y_lim = ceil(data_y_lim / 5) * 5;
    ylim([data_y_lim - plot_mag_dr, data_y_lim]);

    title(title_str, 'Interpreter', 'None');
    lgd = legend(lgd_str, 'Interpreter', 'None', 'Location', 'Best');
    title(lgd, sprintf('%s %0.f^{th} to %0.f^{th} percentiles  (horizontal)  (ears average)', ...
        lgd_title_str, percentiles(1), percentiles(end)), ...
        'Interpreter', 'Tex');
    clear ax data_y_lim data_op;
    drawnow();

    %%
    nexttile(tl);
    data_y_lim = -inf;
    data_op = @max;
    for d = 1 : n_data
        hold on;
        plot_spec_azims_percentiles(freqs_array{d}, spec_hor_perc_array{d}, ...
            percentiles, '', 1, plot_freq_rng, ...
            do_normalize_mean, norm_freq_rng, colors(d, :));
        ax = gca;
        freq_rng = get_freqs_inclusive(freqs_array{d}, MAG_LIM_FREQ_RNG);
        data_y_lim = data_op([data_y_lim; data_op(ax.Children(1).YData(freq_rng))]);
        % remove patch component and only keep line
        delete(ax.Children(1));
        set(ax, 'YAxisLocation', 'Right');
    end

    data_y_lim = ceil(data_y_lim / 5) * 5;
    ylim([data_y_lim - plot_mag_dr, data_y_lim]);

    title(title_str, 'Interpreter', 'None');
    lgd = legend(lgd_str, 'Interpreter', 'None', 'Location', 'Best');
    title(lgd, [lgd_title_str, ' median  (horizontal)  (ears average)'], ...
        'Interpreter', 'Tex');
    clear ax data_y_lim data_op;
    drawnow();
end
