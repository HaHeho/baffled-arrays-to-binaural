function plot_spec_compare(lgd_str, ...
    irs_fro_array, spec_hor_array, spec_all_array, fs_vec, ...
    is_difference, apply_frac_sm, do_normalize_mean, norm_freq_rng)
    
    ETC_DR = 100; % in dB
    MAG_DR = 40; % in dB
    MAG_LIM_FREQ_RNG = [10, 20e3]; % in Hz

    if nargin < 6; is_difference = false; end % true or false
    if nargin < 7; apply_frac_sm = 12; end % in fractional octaves
    if nargin < 8; do_normalize_mean = true; end % true or false
    if nargin < 9; norm_freq_rng = [20, 1e3]; end % in Hz

    tl = tiledlayout(2, 2, 'TileSpacing', 'tight', 'Padding', 'tight');
    drawnow;

    is_show_all = ~isempty(spec_all_array);
    if ~is_show_all
        spec_all_array = spec_hor_array;
    end
    n_data = length(lgd_str);
    if n_data ~= length(irs_fro_array) || n_data ~= length(spec_hor_array) ...
            || n_data ~= length(spec_all_array)
        error('Mismatch in provided cell sizes.');
    end
    if length(fs_vec) == 1
        fs_vec = repmat(fs_vec, 1, n_data);
    end

    colors = AKcolormaps('Dark2', n_data);
    if ~is_difference
        colors = [0, 0, 0; colors]; % add black as first color
        colors = colors(1:end-1, :); % remove last color
    end

    for d = 1 : n_data
        if is_difference
            % in this case one-sided spectra should be provided here, where
            % smoothing has already been applied before
            spec_fro_array{d} = irs_fro_array{d};
        else
            etc_fro_array{d} = db(abs(irs_fro_array{d})); %#ok<*AGROW> 
            spec_fro_array{d} = AKboth2singleSidedSpectrum(fft(irs_fro_array{d}));
            if apply_frac_sm
                spec_fro_array{d} = AKfractOctSmooth(spec_fro_array{d}, 'amp', fs_vec(d), apply_frac_sm);
            end
        end
        spec_fro_array{d} = db(abs(spec_fro_array{d}));
        spec_hor_array{d} = db(abs(spec_hor_array{d}));
        spec_all_array{d} = db(abs(spec_all_array{d}));
        
        if ~is_difference
            times_array{d} = 1:size(etc_fro_array{d}, 1);
        end
        freqs_array{d} = linspace(0, fs_vec(d)/2, size(spec_fro_array{d}, 1));
        
        if do_normalize_mean
            norm_freqs = get_freqs_inclusive(freqs_array{d}, norm_freq_rng);
            
            norm_lvl(d) = mean(spec_fro_array{d}(norm_freqs, :), 'all');
            spec_fro_array{d} = spec_fro_array{d} - norm_lvl(d);

            norm_lvl(d) = mean(spec_hor_array{d}(norm_freqs, :), 'all');
            spec_hor_array{d} = spec_hor_array{d} - norm_lvl(d);

            norm_lvl(d) = mean(spec_all_array{d}(norm_freqs, :), 'all');
            spec_all_array{d} = spec_all_array{d} - norm_lvl(d);
        end
    end

    % adjust shown minimum frequency
    norm_freq_rng(1) = max([norm_freq_rng(1), min(cellfun(@(x) x(2), freqs_array))]);

    %%
    nexttile(tl);
    hold on;
    if is_difference
        axis off;
    else
        data_x_lim = 0; data_y_lim = -inf;
        for d = 1 : n_data
            h = plot(times_array{d}, etc_fro_array{d}, ...
                'LineWidth', 1, 'LineStyle', '-', 'Color', colors(d, :));
            h(1).DisplayName = lgd_str{d};
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
        set(gca, 'Layer', 'Top');

        xlabel('Time in samples');
        ylabel('Magnitude in dB');
        title('Energy Time Curve', 'Interpreter', 'None');

        lgd = legend('Interpreter', 'None', 'Location', 'NorthEast');
        title(lgd, 'HRIRs (frontal)');
        clear data_x_lim data_y_lim;
        drawnow;
    end

    %%
    nexttile(tl);
    if is_difference
        data_y_lim = inf;
        data_op = @min;
    else
        data_y_lim = -inf;
        data_op = @max;
    end
    for d = 1 : n_data
        h = semilogx(freqs_array{d}, spec_fro_array{d}, ...
            'LineWidth', 2, 'LineStyle', '-', 'Color', colors(d, :));
        h(1).DisplayName = lgd_str{d};
        if size(irs_fro_array{d}, 2) > 1
            h(2).LineStyle = ':';
            h(2).LineWidth = 2.5;
            h(2).HandleVisibility = 'off';
        end
        freq_rng = get_freqs_inclusive(freqs_array{d}, MAG_LIM_FREQ_RNG);
        data_y_lim = data_op([data_y_lim; arrayfun(@(y) data_op(y.YData(freq_rng)), h)]);
        hold on;
    end

    xlim([norm_freq_rng(1), freqs_array{end}(end)]);
    if is_difference
        data_y_lim = floor(data_y_lim / 5) * 5;
        ylim([data_y_lim, data_y_lim + MAG_DR]);
    else
        data_y_lim = ceil(data_y_lim / 5) * 5;
        ylim([data_y_lim - MAG_DR, data_y_lim]);
    end

    [ticks, ticklabels] = get_freqs_ticks();
    xticks(ticks);
    xticklabels(ticklabels);
    
    grid on;
    box on;
    set(gca, 'Layer', 'Top');

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

    lgd = legend('Interpreter', 'None', 'Location', 'SouthWest');
    lgd_title_str = 'HRIRs ';
    if is_difference
        lgd_title_str = [lgd_title_str, 'difference '];
    end
    title(lgd, [lgd_title_str, '(frontal)']);

    clear data_y_lim data_op;
    drawnow;
    
    %%
    nexttile(tl);
    if is_show_all
        if is_difference
            data_y_lim = inf;
            data_op = @min;
        else
            data_y_lim = -inf;
            data_op = @max;
        end
        for d = 1 : n_data
            h = semilogx(freqs_array{d}, spec_all_array{d}, ...
                'LineWidth', 2, 'LineStyle', '-', 'Color', colors(d, :));
            h(1).DisplayName = lgd_str{d};
            if size(spec_all_array{d}, 2) > 1
                h(2).LineStyle = ':';
                h(2).LineWidth = 2.5;
                h(2).HandleVisibility = 'off';
            end
            freq_rng = get_freqs_inclusive(freqs_array{d}, MAG_LIM_FREQ_RNG);
            data_y_lim = data_op([data_y_lim; arrayfun(@(y) data_op(y.YData(freq_rng)), h)]);
            hold on;
        end

        xlim([norm_freq_rng(1), freqs_array{end}(end)]);
        if is_difference
            data_y_lim = floor(data_y_lim / 5) * 5;
            ylim([data_y_lim, data_y_lim + MAG_DR]);
        else
            data_y_lim = ceil(data_y_lim / 5) * 5;
            ylim([data_y_lim - MAG_DR, data_y_lim]);
        end

        [ticks, ticklabels] = get_freqs_ticks();
        xticks(ticks);
        xticklabels(ticklabels);

        grid on;
        box on;
        set(gca, 'Layer', 'Top');
    
        xlabel('Frequency in Hz');
        ylabel('Magnitude in dB');
        title(title_str, 'Interpreter', 'None');

        lgd = legend('Interpreter', 'None', 'Location', 'SouthWest');
        title(lgd, [lgd_title_str, 'average (spherical)']);
        clear data_y_lim data_op;
        drawnow;
    else
        axis off;
    end

    %%
    nexttile(tl);
    if is_difference
        data_y_lim = inf;
        data_op = @min;
    else
        data_y_lim = -inf;
        data_op = @max;
    end
    for d = 1 : n_data
        h = semilogx(freqs_array{d}, spec_hor_array{d}, ...
            'LineWidth', 2, 'LineStyle', '-', 'Color', colors(d, :));
        h(1).DisplayName = lgd_str{d};
        if size(spec_hor_array{d}, 2) > 1
            h(2).LineStyle = ':';
            h(2).LineWidth = 2.5;
            h(2).HandleVisibility = 'off';
        end
        freq_rng = get_freqs_inclusive(freqs_array{d}, MAG_LIM_FREQ_RNG);
        data_y_lim = data_op([data_y_lim; arrayfun(@(y) data_op(y.YData(freq_rng)), h)]);
        hold on;
    end

    xlim([norm_freq_rng(1), freqs_array{end}(end)]);
    if is_difference
        data_y_lim = floor(data_y_lim / 5) * 5;
        ylim([data_y_lim, data_y_lim + MAG_DR]);
    else
        data_y_lim = ceil(data_y_lim / 5) * 5;
        ylim([data_y_lim - MAG_DR, data_y_lim]);
    end

    [ticks, ticklabels] = get_freqs_ticks();
    xticks(ticks);
    xticklabels(ticklabels);
    
    grid on;
    box on;
    set(gca, 'Layer', 'Top');

    xlabel('Frequency in Hz');
    ylabel('Magnitude in dB');
    title(title_str, 'Interpreter', 'None');

    lgd = legend('Interpreter', 'None', 'Location', 'SouthWest');
    title(lgd, [lgd_title_str, 'average (horizontal)']);
    clear data_y_lim data_op;
    drawnow;
end
