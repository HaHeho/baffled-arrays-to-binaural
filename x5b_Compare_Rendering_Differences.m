% Compare arbitrary combinations of (rendered) binaural room impulse 
% responses by plotting various time domain and averaged frequency domain 
% representations. This is helpful for directly comparing similar rendering 
% configurations and validating the rendering method or detecting 
% potentially flawed data sets.
% 
% -------------------------------------------------------------------------
%
% requires AKtools toolbox (run AKtoolsStart.m)
% $ git clone https://github.com/f-brinkmann/AKtools.git
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
configs_base_dir = 'resources/BRIR_rendered/Kemar_HRTF_sofa_N44_adjusted';

% all configurations by order
configs = { ... 
    {'Kemar_HRTF_sofa_N44_adjusted_SSR_SSR.wav', ...
    {'Simulation_SMA_TD14_SrcEar*_SSR.wav', 'Simulation_EMA5_SrcEar*_SSR.wav'}, ...
    {'Simulation_SMA_TD42_SrcEar*_SSR.wav', 'Simulation_EMA9_SrcEar*_SSR.wav'}, ...
    {'Simulation_SMA_TD146_SrcEar*_SSR.wav', 'Simulation_EMA17_SrcEar*_SSR.wav'}, ...
    {'Simulation_SMA_TD314_SrcEar*_SSR.wav', 'Simulation_EMA25_SrcEar*_SSR.wav', 'Simulation_SMA_LE2702_SrcEar*_SSR.wav', 'Simulation_EMA89_SrcEar*_SSR.wav'} ...
    }, {'Anechoic_KEMAR_SrcEar_SSR.wav', ...
    {'Anechoic_SMA_TD14_SrcEar*_SSR.wav', 'Anechoic_EMA5_SrcEar*_SSR.wav', 'Anechoic_XMA6_SrcEar*_SSR.wav'}, ...
    {'Anechoic_SMA_TD42_SrcEar*_SSR.wav', 'Anechoic_EMA9_SrcEar*_SSR.wav', 'Anechoic_XMA9_SrcEar*_SSR.wav'}, ...
    {'Anechoic_SMA_TD146_SrcEar*_SSR.wav', 'Anechoic_EMA17_SrcEar*_SSR.wav', 'Anechoic_XMA18_SrcEar*_SSR.wav'}, ...
    {'Anechoic_SMA_TD314_SrcEar*_SSR.wav', 'Anechoic_EMA25_SrcEar*_SSR.wav', 'Anechoic_SMA_LE2702_SrcEar*_SSR.wav', 'Anechoic_EMA89_SrcEar*_SSR.wav'} ...
    }, {'LabWet_KEMAR_SrcEar_SSR.wav', ...
    {'LabWet_SMA_TD14_SrcEar*_SSR.wav', 'LabWet_EMA5_SrcEar*_SSR.wav', 'LabWet_XMA6_SrcEar*_SSR.wav'}, ...
    {'LabWet_SMA_TD42_SrcEar*_SSR.wav', 'LabWet_EMA9_SrcEar*_SSR.wav', 'LabWet_XMA9_SrcEar*_SSR.wav'}, ...
    {'LabWet_SMA_TD146_SrcEar*_SSR.wav', 'LabWet_EMA17_SrcEar*_SSR.wav', 'LabWet_XMA18_SrcEar*_SSR.wav'}, ...
    {'LabWet_SMA_TD314_SrcEar*_SSR.wav', 'LabWet_EMA25_SrcEar*_SSR.wav', 'LabWet_SMA_LE2702_SrcEar*_SSR.wav', 'LabWet_EMA89_SrcEar*_SSR.wav'} ...
}};
% configs = { ...
%     {'Kemar_HRTF_sofa_N44_adjusted_SSR_SSR.wav', ...
%     {'Simulation_SMA_TD14_SrcEar*_SSR.wav', 'Simulation_EMA5_SrcEar*_SSR.wav'}, ...
%     {'Simulation_SMA_TD42_SrcEar*_SSR.wav', 'Simulation_EMA9_SrcEar*_SSR.wav'}, ...
%     {'Simulation_SMA_TD146_SrcEar*_SSR.wav', 'Simulation_EMA17_SrcEar*_SSR.wav'}, ...
%     {'Simulation_SMA_TD314_SrcEar*_SSR.wav', 'Simulation_EMA25_SrcEar*_SSR.wav', 'Simulation_SMA_LE2702_SrcEar*_SSR.wav', 'Simulation_EMA89_SrcEar*_SSR.wav'} ...
% }};

norm_freq_rng = [40, 500]; % in Hz, 500 Hz is around the aliasing frequency of SH1
plot_freq_rng = [40, 24e3]; % in Hz
plot_percentiles = [25, 75]; % in percent
% plot_percentiles = [5, 95]; % in percent
plot_frac_sm = 24; % in fractional octaves
plot_fig_size = [40, 25];
plot_resolution = 300; % in DPI
plot_mag_dr = 45; % in dB
plot_export_dir = 'plots/5_Rendering';

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
