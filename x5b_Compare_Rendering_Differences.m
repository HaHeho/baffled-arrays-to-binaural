% Compare arbitrary combinations of (rendered) binaural room impulse 
% responses by plotting various time domain and averaged frequency domain
% reprensetations. This is helpful for the direct comparison of similar
% rendering configurations and to validate the rendering method or detect
% potentially flawed data sets.
% 
% -------------------------------------------------------------------------
%
% requires AKtools toolbox (run AKtoolsStart.m)
% $ svn checkout https://svn.ak.tu-berlin.de/svn/AKtools --username aktools --password ak
%
% requires `natsortfiles()` from
% https://se.mathworks.com/matlabcentral/fileexchange/47434-natural-order-filename-sort
%
% -------------------------------------------------------------------------
%
% Hannes Helmholz 19.04.2023
%
% -------------------------------------------------------------------------
clear; clc; close all;

addpath(genpath('dependencies'));

%% configuration
config_azims_deg = 0 : 359; % in deg
config_azim_front_id = find(config_azims_deg == 0); % index

% horizontal rotation
[configs_base_dir, plot_export_dir, plot_mag_dr, excluded_configs, configs] = deal( ....
    'resources/BRIR_rendered/HRIR_L2702', 'plots/5_Rendering/HRIR_L2702/compare_all/', ...
    55, {'vertical'}, { ...
    {'SIM_VSA_SMA_LE1202_PW_struct_SubSampled_SH44_SSR.wav', ... 'HRIR_L2702_SSR360_SSR.wav', ...
    {'SIM_VSA_SMA_LE14_PW*_SSR.wav'}, {'SIM_VSA_SMA_TD14_PW*_SSR.wav'}, {'SIM_VSA_EMA5_PW*_SSR.wav'}, ...
    {'SIM_VSA_SMA_LE38_PW*_SSR.wav'}, {'SIM_VSA_SMA_TD42_PW*_SSR.wav'}, {'SIM_VSA_EMA9_PW*_SSR.wav'}, {'Simulation_SMA_EM32_SrcEar*_SSR.wav'}, ...
    {'SIM_VSA_SMA_LE110_PW*_SSR.wav'}, {'SIM_VSA_SMA_TD146_PW*_SSR.wav'}, {'SIM_VSA_EMA17_PW*_SSR.wav'}, ...
    {'SIM_VSA_SMA_LE230_PW*_SSR.wav'}, {'SIM_VSA_SMA_TD314_PW*_SSR.wav'}, {'SIM_VSA_EMA25_PW*_SSR.wav'}, ...
    {'SIM_VSA_SMA_LE1202_PW*_SSR.wav', 'SIM_VSA_EMA59_PW*_SSR.wav', 'HRIR_L2702_SSR360_SSR.wav'} ...
    }, {'CR1_VSA_1202RS_L_struct_SSR.wav', ... 'BRIR_CR1_KU_ROTM_L_SSR.wav', ...
    {'CR1_VSA_SMA_LE14_L*_SSR.wav'}, {'CR1_VSA_SMA_TD14_L*_SSR.wav'}, {'CR1_VSA_EMA5_L*_SSR.wav'}, ...
    {'CR1_VSA_SMA_LE38_L*_SSR.wav'}, {'CR1_VSA_SMA_TD42_L*_SSR.wav'}, {'CR1_VSA_EMA9_L*_SSR.wav'}, ...
    {'CR1_VSA_SMA_LE110_L*_SSR.wav'}, {'CR1_VSA_SMA_TD146_L*_SSR.wav'}, {'CR1_VSA_EMA17_L*_SSR.wav'}, ...
    {'CR1_VSA_SMA_LE230_L*_SSR.wav'}, {'CR1_VSA_SMA_TD314_L*_SSR.wav'}, {'CR1_VSA_EMA25_L*_SSR.wav'}, ...
    {'CR1_VSA_1202RS_L*_SSR.wav', 'CR1_VSA_EMA59_L*_SSR.wav', 'BRIR_CR 1_KU_ROTM_L_SSR.wav'} ...
    }, {'LBS_VSA_1202RS_PAC_struct_SSR.wav', ... 'BRIR_LBS_KU_ROTM_PAC_SSR.wav', ...
    {'LBS_VSA_SMA_LE14_PAC*_SSR.wav'}, {'LBS_VSA_SMA_TD14_PAC*_SSR.wav'}, {'LBS_VSA_EMA5_PAC*_SSR.wav'}, ...
    }, {'LBS_VSA_1202RS_PAC_struct_SSR.wav', ... % split configurations as they are too memory intensive
    {'LBS_VSA_SMA_LE38_PAC*_SSR.wav'}, {'LBS_VSA_SMA_TD42_PAC*_SSR.wav'}, {'LBS_VSA_EMA9_PAC*_SSR.wav'}, ...
    }, {'LBS_VSA_1202RS_PAC_struct_SSR.wav', ... % split configurations as they are too memory intensive
    {'LBS_VSA_SMA_LE110_PAC*_SSR.wav'}, {'LBS_VSA_SMA_TD146_PAC*_SSR.wav'}, {'LBS_VSA_EMA17_PAC*_SSR.wav'}, ...
    }, {'LBS_VSA_1202RS_PAC_struct_SSR.wav', ... % split configurations as they are too memory intensive
    {'LBS_VSA_SMA_LE230_PAC*_SSR.wav'}, {'LBS_VSA_SMA_TD314_PAC*_SSR.wav'}, {'LBS_VSA_EMA25_PAC*_SSR.wav'}, ...
    }, {'LBS_VSA_1202RS_PAC_struct_SSR.wav', ... % split configurations as they are too memory intensive
    {'LBS_VSA_1202RS_PAC*_SSR.wav', 'LBS_VSA_EMA59_PAC*_SSR.wav', 'BRIR_LBS_KU_ROTM_PAC_SSR.wav'} ...
}});

% vertical rotation
[configs_base_dir, plot_export_dir, plot_mag_dr, excluded_configs, configs] = deal( ....
    'resources/BRIR_rendered/HRIR_L2702', 'plots/5_Rendering/HRIR_L2702/compare_all/', ...
    55, {}, { ...
    {'SIM_VSA_SMA_LE1202_PW_struct_SubSampled_SH44_vertical_SSR.wav', ...
    {'SIM_VSA_SMA_LE14_PW_*_vertical_SSR.wav', 'SIM_VSA_SMA_TD14_PW_*_vertical_SSR.wav', 'SIM_VSA_EMA5_PW_*_vertical_SSR.wav'}, ...
    {'SIM_VSA_SMA_LE38_PW_*_vertical_SSR.wav', 'SIM_VSA_SMA_TD42_PW_*_vertical_SSR.wav', 'SIM_VSA_EMA9_PW_*_vertical_SSR.wav', 'Simulation_SMA_EM32_SrcEar_*_vertical_SSR.wav'}, ...
    {'SIM_VSA_SMA_LE110_PW_*_vertical_SSR.wav', 'SIM_VSA_SMA_TD146_PW_*_vertical_SSR.wav', 'SIM_VSA_EMA17_PW_*_vertical_SSR.wav'}, ...
    {'SIM_VSA_SMA_LE230_PW_*_vertical_SSR.wav', 'SIM_VSA_SMA_TD314_PW_*_vertical_SSR.wav', 'SIM_VSA_EMA25_PW_*_vertical_SSR.wav'}, ...
    {'SIM_VSA_SMA_LE1202_PW_*_vertical_SSR.wav', 'SIM_VSA_EMA59_PW_*_vertical_SSR.wav', 'HRIR_L2702_vertical_SSR360_SSR.wav'} ...
    }, {'CR1_VSA_1202RS_L_struct_vertical_SSR.wav', ...
    {'CR1_VSA_SMA_LE14_L_*_vertical_SSR.wav', 'CR1_VSA_SMA_TD14_L_*_vertical_SSR.wav', 'CR1_VSA_EMA5_L_*_vertical_SSR.wav'}, ...
    {'CR1_VSA_SM A_LE38_L_*_vertical_SSR.wav', 'CR1_VSA_SMA_TD42_L_*_vertical_SSR.wav', 'CR1_VSA_EMA9_L_*_vertical_SSR.wav'}, ...
    {'CR1_VSA_SMA_LE110_L_*_vertical_SSR.wav', 'CR1_VSA_SMA_TD146_L_*_vertical_SSR.wav', 'CR1_VSA_EMA17_L_*_vertical_SSR.wav'}, ...
    {'CR1_VSA_SMA_LE230_L_*_vertical_SSR.wav', 'CR1_VSA_SMA_TD314_L_*_vertical_SSR.wav', 'CR1_VSA_EMA25_L_*_vertical_SSR.wav'}, ...
    {'CR1_VSA_1202RS_L_*_vertical_SSR.wav', 'CR1_VSA_EMA59_L_*_vertical_SSR.wav'} ...
    }, {'LBS_VSA_1202RS_PAC_struct_vertical_SSR.wav', ... % split configurations as they are too memory intensive
    {'LBS_VSA_SMA_LE14_PAC_*_vertical_SSR.wav', 'LBS_VSA_SMA_TD14_PAC_*_vertical_SSR.wav', 'LBS_VSA_EMA5_PAC_*_vertical_SSR.wav'}, ...
    }, {'LBS_VSA_1202RS_PAC_struct_vertical_SSR.wav', ... % split configurations as they are too memory intensive
    {'LBS_VSA_SMA_LE38_PAC_*_vertical_SSR.wav', 'LBS_VSA_SMA_TD42_PAC_*_vertical_SSR.wav', 'LBS_VSA_EMA9_PAC_*_vertical_SSR.wav'}, ...
    }, {'LBS_VSA_1202RS_PAC_struct_vertical_SSR.wav', ... % split configurations as they are too memory intensive
    {'LBS_VSA_SMA_LE110_PAC_*_vertical_SSR.wav', 'LBS_VSA_SMA_TD146_PAC_*_vertical_SSR.wav', 'LBS_VSA_EMA17_PAC_*_vertical_SSR.wav'}, ...
    }, {'LBS_VSA_1202RS_PAC_struct_vertical_SSR.wav', ... % split configurations as they are too memory intensive
    {'LBS_VSA_SMA_LE230_PAC_*_vertical_SSR.wav', 'LBS_VSA_SMA_TD314_PAC_*_vertical_SSR.wav', 'LBS_VSA_EMA25_PAC_*_vertical_SSR.wav'}, ...
    }, {'LBS_VSA_1202RS_PAC_struct_vertical_SSR.wav', ... % split configurations as they are too memory intensive
    {'LBS_VSA_1202RS_PAC_*_vertical_SSR.wav', 'LBS_VSA_EMA59_PAC_*_vertical_SSR.wav'} ...
}});

norm_freq_rng = [40, 500]; % in Hz, 500 Hz is around the aliasing frequency of SH1
plot_freq_rng = [40, 24e3]; % in Hz
plot_percentiles = [25, 75]; % in percent
% plot_percentiles = [5, 95]; % in percent
plot_frac_sm = 24; % in fractional octaves
plot_fig_size = [40, 25];
plot_resolution = 300; % in DPI

global DO_PLOT_DIFF DO_EXPORT_PLOT %#ok<*GVMIS>
DO_PLOT_DIFF = true;
DO_EXPORT_PLOT = true;

% DO_PLOT_DIFF = false;
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
            break;
        end
    end; clear c ci new_files;
    fprintf('found %d files ... ', length(files));
    fprintf('done.\n');
    
    if isempty(files)
        fprintf('Configuration skipped (did not find reference file "%s").\n\n', config);
        continue;
    end; clear config;

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
    end; clear f file file_name sig;
    fprintf('\n');
    
    %% compute spectra
    if DO_PLOT_DIFF
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
    end

    fprintf('Gathering frontal IRs ... ');
    for f = 1 : length(files)
        files(f).sig_fro = files(f).sig(:, config_azim_front_id, 1); % left ear only
    end; clear f;
    fprintf('done.\n');
    
    fprintf('Applying spectral smoothing ... ');
    fprintf('of 1/%d octaves ... ', plot_frac_sm);
    for f = 1 : length(files)
        files(f).spec_sm = AKboth2singleSidedSpectrum(fft(files(f).sig));
        files(f).sig = []; % clear variables for lower RAM usage
        % apply spectral smoothing
        for e = 1 : size(files(f).spec_sm, 3)
            files(f).spec_sm(:, :, e) = AKfractOctSmooth(files(f).spec_sm(:, :, e), ...
                'amp', files(f).fs, plot_frac_sm);
        end; clear e;
    end; clear f;
    fprintf('done.\n');

    if DO_PLOT_DIFF
        fprintf('Computing spectral differences after smoothing ... ');
        for f = 2 : length(files)
            % calculate spectral difference
            files(f).spec_sm = compute_spectral_difference(files(f).spec_sm, files(1).spec_sm);
        end; clear f;
        files(1).spec_sm = ones(size(files(1).spec_sm));
        fprintf('done.\n');
    end

    fprintf('Computing spectral difference percentiles ... ');
    for f = 1 : length(files)
        files(f).spec_sm_perc = compute_spectral_percentiles( ...
            files(f).spec_sm, 2, plot_percentiles);
        files(f).spec_sm = []; % clear variables for lower RAM usage
    end; clear f;
    fprintf('done.\n');

    fprintf('Computing spectral difference ear average ... ');
    for f = 1 : length(files)
        files(f).spec_sm_perc = compute_spectral_average(files(f).spec_sm_perc, 3);
    end; clear f;
    fprintf('done.\n');

    %% generate plot
    for c = 2 : length(config_group)
        config_files = files([files.group] == 1 | [files.group] == c);
        if size(config_files) <= 1; continue; end
    
        if DO_PLOT_DIFF
            fig_name = [config_files(2).config, '_difference'];
        else
            fig_name = [config_files(2).config, '_compare'];
        end
        fprintf('Generating plot "%s" ... ', fig_name);
        if plot_fig_size %#ok<BDLGI,BDSCI>
            fig = AKf(plot_fig_size(1), plot_fig_size(2));
        else
            fig = AKf();
        end
        set(fig, 'NumberTitle', 'Off', 'Name', fig_name);
    
        plot_spec_compare_percentiles({config_files.config}, 'BRIRs', ...
            {config_files.sig_fro}, {config_files.spec_sm_perc}, ...
            plot_percentiles, [config_files.fs], plot_frac_sm, true, ...
            norm_freq_rng, plot_freq_rng, plot_mag_dr);

        if DO_EXPORT_PLOT %#ok<*UNRCH>
            fprintf('exporting ... ');
            [~, ~] = mkdir(plot_export_dir); % ignore warning if directory already exists
            fig_name = fullfile(plot_export_dir, [fig_name, '.pdf']);
            exportgraphics(fig, fig_name, 'Resolution', plot_resolution);
        end
        fprintf('done.\n');
    
    end; clear c config_files fig_name fig;
    fprintf('\n');

    if DO_EXPORT_PLOT; close all; end
end; clear g;

%%
fprintf(' ... finished in %.0fh %.0fm %.0fs.\n', ...
    toc/3600, mod(toc,3600)/60, mod(toc,60));
