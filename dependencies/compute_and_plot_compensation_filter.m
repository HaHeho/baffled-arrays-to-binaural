function [ir_min, ir_lin, ir_inv_min, ir_inv_lin] = ...
    compute_and_plot_compensation_filter(fig_name, irs, fs, target_len, ...
    comp_freq_rng, comp_frac_sm)

    global DO_EXPORT_PLOT

    if nargin < 6; comp_frac_sm = 24; end % in fractional octaves
    if nargin < 5 || isempty(comp_freq_rng); comp_freq_rng = [0, fs/2]; end % in Hz

%     LF_EXT_FREQS = [0, 130]; % in Hz
    LF_EXT_FREQS = []; % no low-frequency extension
%     RES_WIN_LEN = 0; % in samples
    RES_WIN_LEN = 16; % in samples
    MAG_DR = [-20, 20]; % in dB
    ETC_DR = 100; % in dB

    is_comp_inv = nargout > 2;
    is_multi_ch = size(irs, 2) > 1;

    % transform into spectrum (invert for compensation filter)
    [spec, is_even] = AKboth2singleSidedSpectrum(fft(irs));
    spec = 1 ./ spec;

    freqs = linspace(0, fs/2, size(irs, 1)/2 + 1).';
    spec_bin = freqs < comp_freq_rng(1) | freqs > comp_freq_rng(2);
    
    % average ears
    spec_avg = compute_spectral_average(spec);

    % normalize mean
    norm_fact = mean(spec_avg(~spec_bin));
    spec_avg = spec_avg / norm_fact;
    spec = spec / norm_fact;

    % apply regularization
    spec_lim = spec_avg;
    spec_lim(spec_bin, :) = 1 ... % substitute magnitude and
        .* exp(1j * angle(spec_lim(spec_bin, :))); % leave phase untouched

    % apply smoothing
    spec_sm = AKfractOctSmooth(spec_lim, 'amp', fs, comp_frac_sm);
    
    % interpolate to target length
    freqs_ip = linspace(0, fs/2, target_len/2 + 1);
    spec_ip = interp1(freqs, spec_sm, freqs_ip, 'linear', 'extrap').';

    if isempty(LF_EXT_FREQS)
        % adjust DC bin
        % spec_ip(1, :) = 1; % does not work well since amplitde is adjusted afterwards
        spec_ip(1, :) = abs(spec_ip(2, :));
    else
        % low-frequency extension
        freqs_bin = get_freqs_exclusive(freqs_ip, LF_EXT_FREQS);
        spec_ip(freqs_bin, :) = spec_ip(find(~freqs_bin, 1), :);
    end

    % invert spectrum
    spec_ip_inv = 1 ./ spec_ip;

    % transform into IR
    target_len_is_even = ~mod(target_len, 2);
    ir_filt = ifft(AKsingle2bothSidedSpectrum(spec_ip, target_len_is_even));
    ir_filt_inv = ifft(AKsingle2bothSidedSpectrum(spec_ip_inv, target_len_is_even));

    % make linear phase with fade-in and fade-out window
    ir_lin = circshift(ir_filt, length(ir_filt) / 2);
    ir_inv_lin = circshift(ir_filt_inv, length(ir_filt_inv) / 2);
    if RES_WIN_LEN
        win = cos(linspace(pi/2, 0, RES_WIN_LEN).').^2;
        ir_lin(1:RES_WIN_LEN) = ir_lin(1:RES_WIN_LEN) .* win;
        ir_inv_lin(1:RES_WIN_LEN) = ir_inv_lin(1:RES_WIN_LEN) .* win;

        win = cos(linspace(0, pi/2, RES_WIN_LEN).').^2;
        ir_lin(end-RES_WIN_LEN+1:end) = ir_lin(end-RES_WIN_LEN+1:end) .* win;
        ir_inv_lin(end-RES_WIN_LEN+1:end) = ir_inv_lin(end-RES_WIN_LEN+1:end) .* win;
    end

    % make minimum phase with fade-out window
    ir_min = AKphaseManipulation(ir_filt, fs, 'min', 2, false);
    ir_inv_min = AKphaseManipulation(ir_filt_inv, fs, 'min', 2, false);
    % half length of minimum phase filter
    ir_min = ir_min(1:ceil(target_len / 2), :);
    ir_inv_min = ir_inv_min(1:ceil(target_len / 2), :);
    if RES_WIN_LEN
        win = cos(linspace(0, pi/2, RES_WIN_LEN).').^2;
        ir_min(end-RES_WIN_LEN+1:end) = ir_min(end-RES_WIN_LEN+1:end) .* win;
        ir_inv_min(end-RES_WIN_LEN+1:end) = ir_inv_min(end-RES_WIN_LEN+1:end) .* win;
    end

    % adjust amplitudes to not result in clipping when exporting filter IRs
    max_amp = max(abs([ir_min; ir_lin]));
    if is_comp_inv
        % also consider maximum amplitudes of inverse filters if requested
        max_amp = max([max_amp, max(abs([ir_inv_min; ir_inv_lin]))]);
    end
    fprintf('normalizing filter amplitude by %.1f dB to prevent export clipping ... ', ...
        -db(max_amp));
    ir_min = ir_min / max_amp;
    ir_lin = ir_lin / max_amp;
    ir_inv_min = ir_inv_min / max_amp;
    ir_inv_lin = ir_inv_lin / max_amp;
    % calculate spectra of final filters
    spec_lin = AKboth2singleSidedSpectrum(fft(ir_lin));
    spec_inv_lin = AKboth2singleSidedSpectrum(fft(ir_inv_lin));
    spec_min = AKboth2singleSidedSpectrum(fft(ir_min));
    spec_inv_min = AKboth2singleSidedSpectrum(fft(ir_inv_min));
    % also adjust spectra amplitudes for the plot to be consistent
    spec = spec / max_amp;
    spec_avg = spec_avg / max_amp;
    spec_lim = spec_lim / max_amp;
    spec_sm = spec_sm / max_amp;  
    spec_ip = spec_ip / max_amp;
    spec_ip_inv = spec_ip_inv / max_amp;
    fprintf('done.\n');

    % generate plot
    [fig_dir, fig_name, ~] = fileparts(fig_name);
    fprintf('Generating plot "%s" ... ', fig_name);
    fig = AKf();
    set(fig, 'NumberTitle', 'Off', 'Name', fig_name);
    tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'tight', 'Padding', 'tight');
    drawnow;

    colors = colororder();
    filter_str = {'Filter (lin. phase)', 'Filter (min. phase)'};
    filter_inv_str = {'Inverse Filter (lin. phase)', 'Inverse Filter (min. phase)'};
    
    nexttile(tl, [1, 2]);
    plot_style = {'m2d', 'fs', fs, 'lw', 2, 'dr', MAG_DR, 'c'};
    AKp(get_padded_IR(spec, is_even), plot_style{:}, 'k');
    if is_multi_ch
        AKp(get_padded_IR(spec_avg, is_even), plot_style{:}, colors(1, :));
        lgd_str = {'Input (left)', 'Input (right)', 'Input averaged'};
    else
        lgd_str = {'Input'}; 
    end
    lgd_str = {lgd_str{:}, 'Input regularized', 'Output (interpolated)', filter_str{:}}; %#ok<CCAT>
    AKp(get_padded_IR(spec_lim, is_even), plot_style{:}, colors(2, :));
    if comp_frac_sm
        AKp(get_padded_IR(spec_sm, is_even), plot_style{:}, colors(3, :));
        lgd_str = {lgd_str{1:4}, sprintf( ...
            'Input smoothed (1/%d frac. oct.)', comp_frac_sm), lgd_str{5:end}};
    end
    AKp(get_padded_IR(spec_ip, is_even), plot_style{:}, colors(6, :));
    AKp(get_padded_IR(spec_lin, is_even), plot_style{:}, colors(4, :), 'lw', 4);
    AKp(get_padded_IR(spec_min, is_even), plot_style{:}, colors(5, :), 'lw', 4);
    if is_comp_inv
        AKp(get_padded_IR(spec_ip_inv, is_even), plot_style{:}, ...
            colors(6, :), 'lw', 3, 'ls', ':');
        AKp(get_padded_IR(spec_inv_lin, is_even), plot_style{:}, ...
            colors(4, :), 'lw', 5, 'ls', ':');
        AKp(get_padded_IR(spec_inv_min, is_even), plot_style{:}, ...
            colors(5, :), 'lw', 5, 'ls', ':');
        lgd_str = [lgd_str, 'Inverse Output (interpolated)', filter_inv_str{:}];
    end
    xline(comp_freq_rng, 'LineWidth', 2, 'LineStyle', ':', 'Color', colors(2, :));
    legend(lgd_str, 'Location', 'Best');
    drawnow;

    nexttile(tl);
    plot_style = {'et2d', 'fs', fs, 'lw', 2, 'xu', 'n', 'dr', ETC_DR, 'c'};
    AKp(ir_lin, plot_style{:}, colors(4, :));
    AKp(ir_min, plot_style{:}, colors(5, :));
    lgd_str = filter_str;
    if is_comp_inv
        AKp(ir_inv_lin, plot_style{:}, colors(4, :), 'lw', 3, 'ls', ':');
        AKp(ir_inv_min, plot_style{:}, colors(5, :), 'lw', 3, 'ls', ':');
        lgd_str = [lgd_str, filter_inv_str{:}];
    end
    legend(lgd_str, 'Location', 'Best');
    xlim([1, size(ir_inv_lin, 1)]);
    drawnow;

    nexttile(tl);
    plot_style = {'gd2d', 'fs', fs, 'lw', 2, 'du', 'ms', 'c'};
    AKp(ir_lin, plot_style{:}, colors(4, :));
    AKp(ir_min, plot_style{:}, colors(5, :));
    lgd_str = filter_str;
    if is_comp_inv
        AKp(ir_inv_lin, plot_style{:}, colors(4, :), 'lw', 3, 'ls', ':');
        AKp(ir_inv_min, plot_style{:}, colors(5, :), 'lw', 3, 'ls', ':');
        lgd_str = [lgd_str, filter_inv_str{:}];
    end
    legend(lgd_str, 'Location', 'Best');
    drawnow;

    if DO_EXPORT_PLOT
        fprintf('exporting ... ');
        file_name = fullfile(fig_dir, [fig_name, '.pdf']);
        exportgraphics(fig, file_name);
    end
end

function ir = get_padded_IR(spec, is_even)
    IR_LEN = 1024; % in samples

    % make zero-phase
    if ~isreal(spec)
        spec = abs(spec);
    end

%     % get IR from single-sided spectrum
%     ir = ifft(AKsingle2bothSidedSpectrum(spec, is_even));
% 
%     % make linear phase (approximately)
%     ir = circshift(ir, round(size(ir, 1) / 2), 1);
%     
%     % zero-pad
%     assert(size(ir, 1) <= IR_LEN);
%     ir(end:IR_LEN, :) = 0;

    % interpolate to target length
    freqs = linspace(0, 1, size(spec, 1));
    freqs_ip = linspace(0, 1, IR_LEN/2 + 1);
    spec_ip = interp1(freqs, spec, freqs_ip, 'linear', 'extrap');
    if size(spec_ip, 1) == 1
        spec_ip = spec_ip.';
    end
    ir = ifft(AKsingle2bothSidedSpectrum(spec_ip, is_even));
end
