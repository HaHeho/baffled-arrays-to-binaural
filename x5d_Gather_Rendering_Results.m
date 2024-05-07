% Generate plots to compare specified combinations of rendered binaural 
% room impulse responses to a designated reference.
%
% The resulting plots are used in the publication manuscript.
% 
% -------------------------------------------------------------------------
%
% requires AKtools toolbox (run AKtoolsStart.m)
% $ git clone https://github.com/f-brinkmann/AKtools.git
%
% -------------------------------------------------------------------------
%
% Hannes Helmholz 27.03.2023
%
% -------------------------------------------------------------------------
clear; clc; close all;

addpath(genpath('dependencies'));

%% configuration
configs_base_dir = 'resources/BRIR_rendered/HRIR_L2702';
config_azims_deg = 0 : 359;
config_azim_front_id = find(config_azims_deg == 0);

% horizontal rotation
[excluded_configs, configs] = deal({'_vertical', '_MagLS', '_eMagLSinCH', '_eMagLSwDC'}, { ...
    {'SIM_VSA_SMA_LE1202_PW_struct_SubSampled_SH44_SSR.wav', ... 'HRIR_L2702_SSR360_SSR.wav'
    {'SIM_VSA_SMA_LE14_PW*_SSR.wav'}, {'SIM_VSA_SMA_TD14_PW*_SSR.wav'}, {'SIM_VSA_EMA5_PW*_SSR.wav'}, ...
    {'SIM_VSA_SMA_LE38_PW*_SSR.wav'}, {'SIM_VSA_SMA_TD42_PW*_SSR.wav'}, {'SIM_VSA_EMA9_PW*_SSR.wav'}, {'Simulation_SMA_EM32_SrcEar*_SSR.wav'}, ...
    {'SIM_VSA_SMA_LE110_PW*_SSR.wav'}, {'SIM_VSA_SMA_TD146_PW*_SSR.wav'}, {'SIM_VSA_EMA17_PW*_SSR.wav'}, ...
    {'SIM_VSA_SMA_LE230_PW*_SSR.wav'}, {'SIM_VSA_SMA_TD314_PW*_SSR.wav'}, {'SIM_VSA_EMA25_PW*_SSR.wav'}, ...
    }, {'HRIR_L2702_SSR360_SSR.wav', ...
    {'SIM_VSA_SMA_LE1202_PW*_SSR.wav', 'SIM_VSA_EMA59_PW*_SSR.wav'} ...
    }, {'CR1_VSA_1202RS_L_struct_SSR.wav', ... 'BRIR_CR1_KU_ROTM_L_SSR.wav', ...
    {'CR1_VSA_SMA_LE14_L*_SSR.wav'}, {'CR1_VSA_SMA_TD14_L*_SSR.wav'}, {'CR1_VSA_EMA5_L*_SSR.wav'}, ...
    {'CR1_VSA_SMA_LE38_L*_SSR.wav'}, {'CR1_VSA_SMA_TD42_L*_SSR.wav'}, {'CR1_VSA_EMA9_L*_SSR.wav'}, ...
    {'CR1_VSA_SMA_LE110_L*_SSR.wav'},{'CR1_VSA_SMA_TD146_L*_SSR.wav'}, {'CR1_VSA_EMA17_L*_SSR.wav'}, ...
    {'CR1_VSA_SMA_LE230_L*_SSR.wav'},  {'CR1_VSA_SMA_TD314_L*_SSR.wav'}, {'CR1_VSA_EMA25_L*_SSR.wav'}, ...
    }, {'BRIR_CR1_KU_ROTM_L_SSR.wav', ...
    {'CR1_VSA_1202RS_L*_SSR.wav', 'CR1_VSA_EMA59_L*_SSR.wav'}, ... 
    }, {'LBS_VSA_1202RS_PAC_struct_SSR.wav', ... 'BRIR_LBS_KU_ROTM_PAC_SSR.wav', ...
    {'LBS_VSA_SMA_LE14_PAC*_SSR.wav'}, {'LBS_VSA_SMA_TD14_PAC*_SSR.wav'}, {'LBS_VSA_EMA5_PAC*_SSR.wav'}, ...
    }, {'LBS_VSA_1202RS_PAC_struct_SSR.wav', ... % split configurations as they are too memory intensive
    {'LBS_VSA_SMA_LE38_PAC*_SSR.wav'}, {'LBS_VSA_SMA_TD42_PAC*_SSR.wav'}, {'LBS_VSA_EMA9_PAC*_SSR.wav'}, ...
    }, {'LBS_VSA_1202RS_PAC_struct_SSR.wav', ... % split configurations as they are too memory intensive
    {'LBS_VSA_SMA_LE110_PAC*_SSR.wav'}, {'LBS_VSA_SMA_TD146_PAC*_SSR.wav'}, {'LBS_VSA_EMA17_PAC*_SSR.wav'}, ...
    }, {'LBS_VSA_1202RS_PAC_struct_SSR.wav', ... % split configurations as they are too memory intensive
    {'LBS_VSA_SMA_LE230_PAC*_SSR.wav'}, {'LBS_VSA_SMA_TD314_PAC*_SSR.wav'}, {'LBS_VSA_EMA25_PAC*_SSR.wav'}, ...
    }, {'BRIR_LBS_KU_ROTM_PAC_SSR.wav', ...
    {'LBS_VSA_1202RS_PAC*_SSR.wav', 'LBS_VSA_EMA59_PAC*_SSR.wav'}, ...
}});

% vertical rotation
[excluded_configs, configs] = deal({}, { ...
    {'SIM_VSA_SMA_LE1202_PW_struct_SubSampled_SH44_vertical_SSR.wav', ... 'HRIR_L2702_vertical_SSR360_SSR.wav'
    {'SIM_VSA_SMA_LE14_PW*_vertical_SSR.wav'}, {'SIM_VSA_SMA_TD14_PW*_vertical_SSR.wav'}, {'SIM_VSA_EMA5_PW*_vertical_SSR.wav'}, ...
    {'SIM_VSA_SMA_LE38_PW*_vertical_SSR.wav'}, {'SIMffYea_VSA_SMA_TD42_PW*_vertical_SSR.wav'}, {'SIM_VSA_EMA9_PW*_vertical_SSR.wav'}, {'Simulation_SMA_EM32_SrcEar*_vertical_SSR.wav'}, ...
    {'SIM_VSA_SMA_LE110_PW*_vertical_SSR.wav'}, {'SIM_VSA_SMA_TD146_PW*_vertical_SSR.wav'}, {'SIM_VSA_EMA17_PW*_vertical_SSR.wav'}, ...
    {'SIM_VSA_SMA_LE230_PW*_vertical_SSR.wav'}, {'SIM_VSA_SMA_TD314_PW*_vertical_SSR.wav'}, {'SIM_VSA_EMA25_PW*_vertical_SSR.wav'}, ...
    }, {'HRIR_L2702_vertical_SSR360_SSR.wav', ...
    {'SIM_VSA_SMA_LE1202_PW*_vertical_SSR.wav', 'SIM_VSA_EMA59_PW*_vertical_SSR.wav'} ...
    }, {'CR1_VSA_1202RS_L_struct_vertical_SSR.wav', ...
    {'CR1_VSA_SMA_LE14_L*_vertical_SSR.wav'}, {'CR1_VSA_SMA_TD14_L*_vertical_SSR.wav'}, {'CR1_VSA_EMA5_L*_vertical_SSR.wav'}, ...
    {'CR1_VSA_SMA_LE38_L*_vertical_SSR.wav'}, {'CR1_VSA_SMA_TD42_L*_vertical_SSR.wav'}, {'CR1_VSA_EMA9_L*_vertical_SSR.wav'}, ...
    {'CR1_VSA_SMA_LE110_L*_vertical_SSR.wav'}, {'CR1_VSA_SMA_TD146_L*_vertical_SSR.wav'}, {'CR1_VSA_EMA17_L*_vertical_SSR.wav'}, ...
    {'CR1_VSA_SMA_LE230_L*_vertical_SSR.wav'}, {'CR1_VSA_SMA_TD314_L*_vertical_SSR.wav'}, {'CR1_VSA_EMA25_L*_vertical_SSR.wav'}, ...
    }, {'CR1_VSA_1202RS_L_struct_vertical_SSR.wav', ...
    {'CR1_VSA_1202RS_L*_vertical_SSR.wav', 'CR1_VSA_EMA59_L*_vertical_SSR.wav'}, ... 
    }, {'LBS_VSA_1202RS_PAC_struct_vertical_SSR.wav', ...
    {'LBS_VSA_SMA_LE14_PAC*_vertical_SSR.wav'}, {'LBS_VSA_SMA_TD14_PAC*_vertical_SSR.wav'}, {'LBS_VSA_EMA5_PAC*_vertical_SSR.wav'}, ...
    }, {'LBS_VSA_1202RS_PAC_struct_vertical_SSR.wav', ... % split configurations as they are too memory intensive
    {'LBS_VSA_SMA_LE38_PAC*_vertical_SSR.wav'}, {'LBS_VSA_SMA_TD42_PAC*_vertical_SSR.wav'}, {'LBS_VSA_EMA9_PAC*_vertical_SSR.wav'}, ...
    }, {'LBS_VSA_1202RS_PAC_struct_vertical_SSR.wav', ... % split configurations as they are too memory intensive
    {'LBS_VSA_SMA_LE110_PAC*_vertical_SSR.wav'}, {'LBS_VSA_SMA_TD146_PAC*_vertical_SSR.wav'}, {'LBS_VSA_EMA17_PAC*_vertical_SSR.wav'}, ...
    }, {'LBS_VSA_1202RS_PAC_struct_vertical_SSR.wav', ... % split configurations as they are too memory intensive
    {'LBS_VSA_SMA_LE230_PAC*_vertical_SSR.wav'}, {'LBS_VSA_SMA_TD314_PAC*_vertical_SSR.wav'}, {'LBS_VSA_EMA25_PAC*_vertical_SSR.wav'}, ...
    }, {'LBS_VSA_1202RS_PAC_struct_vertical_SSR.wav', ...
    {'LBS_VSA_1202RS_PAC*_vertical_SSR.wav', 'LBS_VSA_EMA59_PAC*_vertical_SSR.wav'}, ...
}});

norm_freq_rng = [40, 200]; % in Hz, 500 Hz is around the aliasing frequency of SH1
plot_freq_rng = [100, 20e3]; % in Hz
plot_mag_rng = [-19.99, 19.99]; % in dB, to not print [-20, 20] axis labels
plot_percentiles = [5, 95]; % in percent
plot_frac_sm = 24; % in fractional octaves
plot_fig_size = [31, 16]; % width to fit the figure title
plot_resolution = 300; % in DPI
plot_export_dir = 'plots/5_Rendering/HRIR_L2702/results_eMagLS';

lgd_str = {sprintf('Raw         median and %0.f^{th} to %0.f^{th} percentile', plot_percentiles), '', ...
    sprintf('eMagLS   median and %0.f^{th} to %0.f^{th} percentile', plot_percentiles), ''};

global DO_EXPORT_PLOT %#ok<*GVMIS>
DO_EXPORT_PLOT = true;
% DO_EXPORT_PLOT = false;

%% 
tic; % start measuring execution time
current_file = fileparts(mfilename('fullpath'));

if ~exist('excluded_configs', 'var'); excluded_configs = {}; end

for g = 1 : length(configs)
    %% load data
    fprintf('Parsing configuration files ... in "%s" ... ', configs_base_dir);
    files = [];
    config_group = configs{g};
    for c = 1 : length(config_group)
        config = config_group{c};
        if iscell(config)
            for ci = 1 : length(config)
                new_files = get_files(configs_base_dir, config{ci}, excluded_configs);
                [new_files(:).group] = deal(c);
                files = [files(:); new_files];
            end
        else
            new_files = get_files(configs_base_dir, config, excluded_configs);
            [new_files(:).group] = deal(c);
            files = [files(:); new_files];
        end
        if c == 1 && length(files) ~= 1
            error('Did not find a reference file.');
        end
    end; clear c ci config new_files;
    fprintf('found %d files ... ', length(files));
    fprintf('done.\n');
    
    for f = 1 : length(files)
        file = fullfile(files(f).folder, files(f).name);
        file = file(length(current_file)+2:end); % remove absolute path
        [~, file_name, ~] = fileparts(file); % remove file path and extension
        if endsWith(file_name, '_SSR', 'IgnoreCase', true)
            file_name = file_name(1:end-4); % remove "_SSR" ending
        end
        files(f).config = file_name; %#ok<*SAGROW> 
    
        fprintf('Loading file "%s" ... ', file);
        [sig, files(f).fs] = audioread(file);
        files(f).sig = unstack_SSR_BRIRs(sig);
        fprintf('done.\n');

        file = [file(1:end-7), 'meta.mat']; % remove "SSR.wav"
        fprintf('Loading file "%s" ... ', file);
        try
            data = load(file, 'freq_SMA_alias');
            files(f).freq_SMA_alias = data.freq_SMA_alias;
            fprintf('done.\n');
        catch
            fprintf('skipped (file not found).\n');
        end
    end; clear f file file_name sig data;
    fprintf('\n');
    
    %% compute spectra
    fprintf('Zero-padding IRs ... ');
    max_len = 0;
    for f = 1 : length(files)
        max_len = max([max_len, size(files(f).sig, 1)]);
    end; clear f;
    fprintf('to %d samples ... ', max_len);
    for f = 1 : length(files)
        files(f).sig(end+1:max_len, :, :) = 0;
    end; clear f;
    fprintf('done.\n');
    
    fprintf('Applying spectral smoothing ... ');
    fprintf('of 1/%d octaves ... ', plot_frac_sm);
    for f = 1 : length(files)
        files(f).spec_sm = AKboth2singleSidedSpectrum(fft(files(f).sig));
        % apply spectral smoothing
        for e = 1 : size(files(f).spec_sm, 3) % each ear
            files(f).spec_sm(:, :, e) = AKfractOctSmooth(files(f).spec_sm(:, :, e), ...
                'amp', files(f).fs, plot_frac_sm);
        end; clear e;
    end; clear f;
    fprintf('done.\n');

    fprintf('Computing spectral differences after smoothing ... ');
    for f = 2 : length(files)
        % calculate spectral difference
        files(f).spec_sm = compute_spectral_difference(files(f).spec_sm, files(1).spec_sm);
    end; clear f;
    files(1).spec_sm = ones(size(files(1).spec_sm));
    fprintf('done.\n');

    fprintf('Computing spectral difference percentiles ... ');
    for f = 1 : length(files)
        files(f).spec_sm_perc = compute_spectral_percentiles( ...
            files(f).spec_sm, 2, plot_percentiles);
        % clear variables for lower RAM usage
        files(f).sig = [];
        files(f).spec_sm = [];
    end; clear f;
    fprintf('done.\n');

    fprintf('Computing spectral difference ear average ... ');
    for f = 1 : length(files)
        files(f).spec_sm_perc_avg = compute_spectral_average(files(f).spec_sm_perc, 3);
    end; clear f;
    fprintf('done.\n');

    %% generate plot
    for c = 2 : length(config_group)
        config_files = files([files.group] == 1 | [files.group] == c);
        if size(config_files) <= 1; continue; end
        fig_name = [config_files(2).config, '_result'];
        size_ear = size(config_files(1).spec_sm_perc, 3);
        lbl_config = {'Interpreter', 'Latex', 'LineWidth', 1, ...
            'HandleVisibility', 'off', 'LabelOrientation', 'horizontal', ...
            'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'right'};

        fprintf('Generating plot "%s" ... ', fig_name);
        if ~isempty(plot_fig_size)
            fig = AKf(plot_fig_size(1), plot_fig_size(2));
        else
            fig = AKf();
        end
        % set(fig, 'NumberTitle', 'Off', 'Name', fig_name, 'Renderer', 'painters');
        set(fig, 'NumberTitle', 'Off', 'Name', fig_name);
        tl = tiledlayout(fig, 2, size_ear, 'TileSpacing', 'tight', 'Padding', 'tight');
        drawnow;

        title_str = 'Spectal difference to BRIR reference';
        if plot_frac_sm
            title_str = sprintf('%s   (1/%d oct. smoothing)', title_str, plot_frac_sm);
        end
        for e = 1 : size_ear % each ear
            ax = nexttile(tl);
            for f = 2 : length(config_files)
                if f == 2; color = 'k'; else; color = 'r'; end
                freqs = linspace(0, config_files(f).fs/2, size(config_files(f).spec_sm_perc, 1));
                plot_spec_azims_percentiles(freqs, config_files(f).spec_sm_perc, ...
                    plot_percentiles, title_str, e, plot_freq_rng, true, norm_freq_rng, color);
            end; clear f color freqs;
            ylabel('Magnitude in dB');
            % legend(lgd_str, 'Location', 'NorthWest');
            legend({config_files(2).config, '', config_files(3).config, ''}, ...
                'Interpreter', 'none', 'Location', 'NorthWest');
            if e == size_ear; set(ax, 'YAxisLocation', 'Right'); end
            xline(ax, config_files(end).freq_SMA_alias, ':', '$f_\textrm{A}$', lbl_config{:});
            ylim(plot_mag_rng);
            drawnow;
        end; clear e;

        title_str = sprintf('%s   (ears average)', title_str);
        ax = nexttile(tl);
        for f = 2 : length(config_files)
            if f == 2; color = 'k'; else; color = 'r'; end
            freqs = linspace(0, config_files(f).fs/2, size(config_files(f).spec_sm_perc_avg, 1));
            plot_spec_azims_percentiles(freqs, config_files(f).spec_sm_perc_avg, ...
                plot_percentiles, title_str, [], plot_freq_rng, true, norm_freq_rng, color);
        end; clear f color freqs;
        ylabel('Magnitude in dB');
        legend(lgd_str, 'Location', 'NorthWest');
        xline(ax, config_files(end).freq_SMA_alias, ':', '$f_\textrm{A}$', lbl_config{:});
        ylim(plot_mag_rng);
        drawnow;

        if DO_EXPORT_PLOT %#ok<*UNRCH>
            fprintf('exporting ... ');
            [~, ~] = mkdir(plot_export_dir); % ignore warning if directory already exists
            fig_name = fullfile(plot_export_dir, [fig_name, '.pdf']);
            exportgraphics(fig, fig_name, 'Resolution', plot_resolution);
        end
        fprintf('done.\n');
    
    end; clear c config_files fig_name size_ear lbl_config fig tl;
    fprintf('\n');

    if DO_EXPORT_PLOT; close all; end
end; clear g;

%%
fprintf(' ... finished in %.0fh %.0fm %.0fs.\n', ...
    toc/3600, mod(toc,3600)/60, mod(toc,60));
