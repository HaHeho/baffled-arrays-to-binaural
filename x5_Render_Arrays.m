% Perform binaural rendering of SMA and EMA impulse response data sets
% at a specified spherical harmonics order and with a desired HRTF set.
%
% Thereby, a reference binaural room impulse response data set (either
% from a dummy head measurement or from a former high-resolution array
% rendering) can be specified. The reference data is compared against the
% currently rendered binaural ear signals by generating extensive plots to
% evaluate time domain and frequency domain differences individually at all 
% rendered head orientations.
%
% The script contains references and functionality for XMA rendering,
% which were left for consistency but should be obsolete for the DAGA
% publication.
%
% -------------------------------------------------------------------------
%
% requires Ambisonic Encoding toolbox
% $ git clone https://github.com/AppliedAcousticsChalmers/ambisonic-encoding.git
%
% requires eMagLS toolbox
% $ git clone https://github.com/thomasdeppisch/eMagLS.git
%
% requires AKtools toolbox (run AKtoolsStart.m)
% $ git clone https://github.com/f-brinkmann/AKtools.git
%
% -------------------------------------------------------------------------
%
% Hannes Helmholz 07.06.2023
%
% -------------------------------------------------------------------------
is_batch_mode = length(dbstack) > 1;
clearvars -except is_batch_mode;
if ~is_batch_mode; clear global; end
clc; close all;

addpath(genpath('dependencies'));

%% configuration
data_export_dir = 'resources/BRIR_rendered';
plot_export_dir = 'plots/5_Rendering';

global DO_PLOT DO_PLOT_PRESEN DO_EXPORT_PLOT %#ok<*GVMIS> 
global DO_EXPORT_META DO_EXPORT_FLAC DO_EXPORT_SSR 
if ~is_batch_mode
    DO_PLOT        = true;
    DO_PLOT_PRESEN = true;
    DO_EXPORT_PLOT = true;
    DO_EXPORT_META = true;
    DO_EXPORT_SSR  = true;
    DO_EXPORT_FLAC = false;

    % DO_PLOT        = false;
    % DO_PLOT_PRESEN = false;
    % DO_EXPORT_PLOT = false;
    % DO_EXPORT_META = false;
    % DO_EXPORT_SSR  = false;
    % DO_EXPORT_FLAC = true;
end

fprintf('Parsing configuration parameters ... ');
global params
if ~isfield(params, 'hrir_file')
    [params.hrir_max_N, params.hrir_file] = deal(44, 'resources/HRIR_KU100/HRIR_L2702.sofa'); % KU100 for DAGA publication for compatibility to reference BRIRs
end
if ~isfield(params, 'sh_type')
    params.sh_type = 'real';
    % params.sh_type = 'complex';
end
if ~isfield(params, 'sh_with_weights')
    params.sh_with_weights = true;
    % params.sh_with_weights = false; % do not use qudrature weights for SH transformations (pseudo-inverse will be used)
end
if ~isfield(params, 'half_block_length')
    params.half_block_length = 8192; % in samples
    % params.half_block_length = 4096; % in samples, required minimum for SMA and EMA
end
if ~isfield(params, 'rf_reg_type')
    params.rf_reg_type = 'tikhonov';
    % params.rf_reg_type = 'soft';
    % params.rf_reg_type = 'hard';
end
if ~isfield(params, 'rf_lim_dB')
    params.rf_lim_dB = 18; % in dB, according to Ambisonics convention
end
if ~isfield(params, 'rf_len')
    params.rf_len = 8192; % in samples, for DAGA publication
    % params.rf_len = 4096; % in samples, to be sure that zero-padded filters are well defined (e.g. at SH2)
    % params.rf_len = 2048; % in samples, recommended minimum for EMA
    % params.rf_len = 1024; % in samples, recommended minimum for SMA
end
if ~isfield(params, 'eq_type')
    params.eq_type = '';
    % params.eq_type = 'MagLS';
    % params.eq_type = 'MagLSwDC';
    % params.eq_type = 'MagLS+SHF';
    % params.eq_type = 'MagLS+ADF'
    % params.eq_type = 'eMagLS';
    % params.eq_type = 'eMagLSwDC';
    % params.eq_type = 'eMagLSinCH';
    % params.eq_type = 'eMagLSinCHwDC';
    % params.eq_type = 'eBreve';
end
if ~isfield(params, 'eq_magls_len')
    params.eq_magls_len = 512; % in samples
end
if ~isfield(params, 'eq_xma_len')
    params.eq_xma_len = 1024; % in samples, obsolete for the DAGA publication
end
if ~isfield(params, 'anec_simulation')
    params.anec_simulation = true;
    % params.anec_simulation = false;
end
if ~isfield(params, 'azim_align')
    params.azim_align = false;
    % params.azim_align = true;
end
if ~isfield(params, 'azim_align_frac_sm')
    params.azim_align_frac_sm = 24; % in fractional octaves
end
if ~isfield(params, 'circ_align')
    params.circ_align = true;
    % params.circ_align = false;
end
if ~isfield(params, 'trunc_threshold')
    params.trunc_threshold = -80; % in dB, relative to global peak
    % params.trunc_threshold = []; % no truncation
end
if ~isfield(params, 'trunc_win_len')
    params.trunc_win_len = 64; % in samples
end
if ~isfield(params, 'plot_fig_size')
    if DO_PLOT_PRESEN
        params.plot_fig_size = [34, 22];
    else
        params.plot_fig_size = [];
    end
end
if ~isfield(params, 'plot_percentiles')
    params.plot_percentiles = [5, 95]; % in percent
end
if ~isfield(params, 'plot_frac_sm')
    params.plot_frac_sm = 24; % in fractional octaves
end
if ~isfield(params, 'plot_freq_rng')
    params.plot_freq_rng = [40, 20e3]; % in Hz
end
if ~isfield(params, 'plot_norm_freq_rng')
    params.plot_norm_freq_rng = [40, 500]; % in Hz
end
if ~isfield(params, 'plot_spec_c_rng')
    if DO_PLOT_PRESEN
        params.plot_spec_c_rng = [-20, 20]; % in dB
    else
        params.plot_spec_c_rng = [-30, 30]; % in dB
    end
end
if ~isfield(params, 'plot_resolution')
    params.plot_resolution = 300; % in DPI
end

if ~isfield(params, 'array_file')
    %%
    % Vertical rotation
    params.reference_file = 'resources/HRIR_KU100/HRIR_L2702_vertical_SSR360.sofa'; % generate reference BRIRs from high-resolution SMA
    [params.N, params.array_file] = deal(29, 'resources/ARIR_WDR/SIM_VSA_SMA_LE1202_PW_struct_SubSampled_SH44.mat');
    % params.reference_file = 'resources/HRIR_KU100/HRIR_L2702_vertical_SSR360.sofa';
    params.reference_file = 'resources/BRIR_rendered/HRIR_L2702/SIM_VSA_SMA_LE1202_PW_struct_SubSampled_SH44_vertical_SSR.wav';
    [params.N, params.array_file] = deal(2, 'resources/ARIR_WDR/SIM_VSA_SMA_LE14_PW_struct_SubSampled_SH44.mat');
    % [params.N, params.array_file] = deal(4, 'resources/ARIR_processed/Simulation_SMA_EM32_SrcEar.sofa'); % 4.2 cm radius
    % [params.N, params.array_file] = deal(4, 'resources/ARIR_WDR/SIM_VSA_SMA_LE38_PW_struct_SubSampled_SH44.mat');
    % [params.N, params.array_file] = deal(8, 'resources/ARIR_WDR/SIM_VSA_SMA_LE110_PW_struct_SubSampled_SH44.mat');
    % [params.N, params.array_file] = deal(12, 'resources/ARIR_WDR/SIM_VSA_SMA_LE230_PW_struct_SubSampled_SH44.mat');
    % [params.N, params.array_file] = deal(2, 'resources/ARIR_WDR/SIM_VSA_EMA5_PW_struct_SubSampled_SH44.mat');
    % [params.N, params.array_file] = deal(4, 'resources/ARIR_WDR/SIM_VSA_EMA9_PW_struct_SubSampled_SH44.mat');
    % [params.N, params.array_file] = deal(8, 'resources/ARIR_WDR/SIM_VSA_EMA17_PW_struct_SubSampled_SH44.mat');
    % [params.N, params.array_file] = deal(12, 'resources/ARIR_WDR/SIM_VSA_EMA25_PW_struct_SubSampled_SH44.mat');
    % [params.N, params.array_file] = deal(29, 'resources/ARIR_WDR/SIM_VSA_EMA59_PW_struct_SubSampled_SH44.mat');

    % % Horizontal rotation
    % params.reference_file = 'resources/HRIR_KU100/HRIR_L2702_SSR360.sofa'; % generate reference BRIRs from high-resolution SMA
    % [params.N, params.array_file] = deal(29, 'resources/ARIR_WDR/SIM_VSA_SMA_LE1202_PW_struct_SubSampled_SH44.mat');
    % params.reference_file = 'resources/BRIR_rendered/HRIR_L2702/SIM_VSA_SMA_LE1202_PW_struct_SubSampled_SH44_SSR.wav';
    % [params.N, params.array_file] = deal(2, 'resources/ARIR_WDR/SIM_VSA_SMA_LE14_PW_struct_SubSampled_SH44.mat');
    % [params.N, params.array_file] = deal(4, 'resources/ARIR_processed/Simulation_SMA_EM32_SrcEar.sofa'); % 4.2 cm radius
    % [params.N, params.array_file] = deal(4, 'resources/ARIR_WDR/SIM_VSA_SMA_LE38_PW_struct_SubSampled_SH44.mat');
    % [params.N, params.array_file] = deal(8, 'resources/ARIR_WDR/SIM_VSA_SMA_LE110_PW_struct_SubSampled_SH44.mat');
    % [params.N, params.array_file] = deal(12, 'resources/ARIR_WDR/SIM_VSA_SMA_LE230_PW_struct_SubSampled_SH44.mat');
    % [params.N, params.array_file] = deal(2, 'resources/ARIR_WDR/SIM_VSA_EMA5_PW_struct_SubSampled_SH44.mat');
    % [params.N, params.array_file] = deal(4, 'resources/ARIR_WDR/SIM_VSA_EMA9_PW_struct_SubSampled_SH44.mat');
    % [params.N, params.array_file] = deal(8, 'resources/ARIR_WDR/SIM_VSA_EMA17_PW_struct_SubSampled_SH44.mat');
    % [params.N, params.array_file] = deal(12, 'resources/ARIR_WDR/SIM_VSA_EMA25_PW_struct_SubSampled_SH44.mat');
    % [params.N, params.array_file] = deal(29, 'resources/ARIR_WDR/SIM_VSA_EMA59_PW_struct_SubSampled_SH44.mat');

    % params.reference_file = 'resources/ARIR_WDR/BRIR_CR1_KU_ROTM_L.sofa'; % generate reference BRIRs from high-resolution SMA
    % [params.N, params.array_file] = deal(29, 'resources/ARIR_WDR/CR1_VSA_SMA_LE1202_L_struct.mat');
    % params.reference_file = 'resources/BRIR_rendered/HRIR_L2702/CR1_VSA_SMA_LE1202_L_struct_SSR.wav';
    % [params.N, params.array_file] = deal(2, 'resources/ARIR_WDR/CR1_VSA_SMA_LE14_L_struct_SubSampled_SH29.mat');
    % [params.N, params.array_file] = deal(4, 'resources/ARIR_WDR/CR1_VSA_SMA_LE38_L_struct_SubSampled_SH29.mat');
    % [params.N, params.array_file] = deal(8, 'resources/ARIR_WDR/CR1_VSA_SMA_LE110_L_struct_SubSampled_SH29.mat');
end

%% verify configuration
tic; % start measuring execution time
SOFAstart;

assert(~strcmpi(params.sh_type, 'complex'), ...
    'Implementation for the "%s" SH convention is not available anymore.', params.sh_type);
assert(params.N <= params.hrir_max_N, ...
    'Rendering for array at SH order %d not possible with HRIR at order %d.', ...
    params.N, params.hrir_max_N);

azEl2azCo = @(grid) [grid(1, :); pi/2 - grid(2, :)]; % in rad
azCo2azEl = azEl2azCo;
params.num_harmonics = (params.N+1)^2;
if contains(params.array_file, 'EMA', 'IgnoreCase', true)
    params.array_type = 'EMA';
    params.is_ch = contains(params.eq_type, 'eMagLSinCH', 'IgnoreCase', true);
    if params.is_ch
        params.num_harmonics = 2*params.N + 1;
    end
elseif contains(params.array_file, {'SMA', 'VSA'}, 'IgnoreCase', true)
    % This has to be done after parsing the EMA case for the exception
    % with the THK WDR "_VSA_" SMAs which have not been subsampled.
    params.array_type = 'SMA';
elseif contains(params.array_file, 'XMA', 'IgnoreCase', true)
    params.array_type = 'XMA';
    if ~strcmpi(params.sh_type, 'real')
        warning('XMA rendering with "%s" SH convention was not validated yet ...', params.sh_type);
    end
else
    error('Unknown array type for "%s".', params.array_file);
end
if (contains(params.eq_type, 'MagLS', 'IgnoreCase', true) && (2 > params.N || params.N > 12)) ...
        || ( strcmpi(params.array_type, 'XMA') && ~any(strcmpi(params.eq_type, {'', 'eBreve'}))) ...
        || (~strcmpi(params.array_type, 'XMA') && strcmpi(params.eq_type, 'eBreve')) ...
        || (~strcmpi(params.array_type, 'EMA') && contains(params.eq_type, 'eMagLSinCH', 'IgnoreCase', true))
    fprintf('skipped rendering for "%s" equalization.\n', params.eq_type);
    return;
end
is_anec_simulation = params.anec_simulation;
if is_anec_simulation
    % Disable anechoic simulation for most configurations
    is_anec_simulation = params.N < 44 && ~contains(params.eq_type, 'eMagLSinCH', 'IgnoreCase', true) ...
        && contains(params.array_file, {'SIM', 'ANE'}, 'IgnoreCase', true);
end
fprintf('found %s configuration ... ', params.array_type);

% parse file names
[~, array_file, ~] = fileparts(params.array_file);
if ~isempty(params.eq_type)
    array_file = sprintf('%s_%s', array_file, params.eq_type);
end
if ~isfield(params, 'reference_file')
    params.reference_file = [];
end
if isempty(params.reference_file)
    reference_file = [];
else
    [~, reference_file, ~] = fileparts(params.reference_file);
end
[hrir_path, hrir_file, ~] = fileparts(params.hrir_file);
params.hrir_eq_file = fullfile(hrir_path, [hrir_file, '_filter_lin.wav']);
clear hrir_path;
data_export_dir = fullfile(data_export_dir, hrir_file);
plot_export_dir = fullfile(plot_export_dir, hrir_file); 
is_first_plot = true;
fprintf('done.\n\n');

%% load array and reference data
fprintf('Loading array file "%s" ... ', params.array_file);
if endsWith(params.array_file, '.sofa')
    tmp = SOFAload(params.array_file);
    array_grid_rad = azEl2azCo(deg2rad(tmp.ReceiverPosition(:, 1:2)).'); % elevation to colatitude
    try
        assert(params.sh_with_weights);
        array_grid_weights = tmp.ReceiverQuadWeight.';
        fprintf('with quadrature weights ... ');
    catch
        fprintf('without quadrature weights ... ');
        array_grid_weights = [];
    end
    array_signals = get_SOFA_IRs(tmp);
    fs = tmp.Data.SamplingRate;
    array_radius = mean(tmp.ReceiverPosition(:, 3));
    try
        array_temp = mean(tmp.AirTemperature); % in degrees Celcius
    catch
        array_temp = [];
    end

    if ~ismatrix(array_signals) % ndims(array_signals) > 2
        % This may happen for the XMA datasets with multiple listener
        % directions
        fprintf('found %d directions ... ', size(array_signals, 2));
        ARRAY_DIR_ID = 1;
        fprintf('picking %.1f deg direction ... ', tmp.ListenerPosition(ARRAY_DIR_ID, 1));
        array_signals = squeeze(array_signals(:, ARRAY_DIR_ID, :));
    end

elseif endsWith(params.array_file, '.mat')
    tmp = load(params.array_file);
    % Make sure to cast data to double in order to prevent limitation of 
    % data size during calculations!!
    array_grid_rad = double([tmp.azimuth; tmp.colatitude]);
    if params.sh_with_weights
        fprintf('with quadrature weights ... ');
        array_grid_weights = double(tmp.quadWeight);
    else
        fprintf('without quadrature weights ... ');
        array_grid_weights = [];
    end
    array_signals = double(tmp.irChOne);
    fs = double(tmp.fs);
    array_radius = double(tmp.radius); % in m
    array_temp = double(tmp.avgAirTemp); % in degrees Celcius

elseif endsWith(params.array_file, '.wav')
    [array_signals, fs] = audioread(params.array_file);
    array_grid_rad = linspace(0, 2*pi, size(array_signals, 2) + 1); % azimuth in rad
    array_grid_rad = array_grid_rad(:, 1:end - 1);
    array_grid_rad(2, :) = pi/2; % colatitude in rad
    array_grid_rad2 = -get_array_grid('SSR', [], true, true).'; % azimuth and colatitude in rad
    array_grid_weights = []; % no qudrature weights available (pseudo-inverse will be used)
    array_radius = []; % in m
    array_temp = []; % in degrees Celcius
    error('The updated array grid needs to be verified!');
    
else
    [~, ~, tmp] = fileparts(params.array_file); 
    error('Data loading for data type "*%s" not implemented yet.', tmp);
end
fprintf('with %d samples ... ', size(array_signals, 1));
clear tmp;
fprintf('done.\n');

is_eMagLS = contains(params.eq_type, 'eMagLS', 'IgnoreCase', true);
if is_eMagLS
    fprintf('Discarding array quadrature weights in case of eMagLS processing ... ');
    if isempty(array_grid_weights)
        fprintf('skipped (no weights were available).\n');
    else
        array_grid_weights = [];
        fprintf('done.\n');
    end
end

fprintf('Zero-padding array signals ... ');
if params.circ_align
%     % This setting would work for all except XMA configurations
%     % (it leads to the circshifted beginnings of the rendered BRIRs 
%     % exhibiting a too high level in order to be cut-off in the end by
%     % `trunc_threshold`)
%     params.circ_len = -(params.rf_len/2 + params.eq_magls_len/2 + params.eq_xma_len); % in samples
%     % This setting is slightly too much for XMA+eBreve configurations
%     params.circ_len = -(params.rf_len + params.eq_magls_len + params.eq_xma_len)/2; % in samples
%     % This setting works well for most XMA configurations
%     params.circ_len = -(params.rf_len + params.eq_xma_len)/2; % in samples
    % This setting is safe for all XMA configurations
    params.circ_len = -params.rf_len/3; % in samples
else
    params.circ_len = 0; % no shifting
end
% longest resulting SMA/EMA BRIRs for 'MagLS' equalization
array_len = size(array_signals, 1) + (params.rf_len - 1) ...
    + (params.eq_magls_len - 1) + (2 * (params.eq_xma_len - 1));
fprintf('for minimum of %d samples ... ', array_len);
array_len = ceil(array_len / params.half_block_length); % to contain multiple of block sizes
% array_len = max(2, array_len) * params.half_block_length; % at least two blocks
array_len = array_len * params.half_block_length;
fprintf('to %d samples ... ', array_len);
array_signals(end+1:array_len, :) = 0;
clear array_len;
fprintf('done.\n');

fprintf('Loading reference BRIR file ');
if isempty(params.reference_file)
    fprintf('... skipped.\n');

    fprintf('Generating head orientation target directions ... ');
    % % Horizontal sweep from 0 deg to -359 deg
    % % (according to the file convention that the SSR requires)
    % array_file = [array_file, '_horizontal'];
    % reference_azel_deg = get_array_grid('SSR', [], false, false).';

    % Vertical sweep from 0 deg to -359 deg
    % (according to the file convention that the SSR requires, for testing
    % purposes, as the scene will not be usable with the regular azimuth 
    % head tracking)
    reference_azel_deg = get_array_grid('SSR_vertical', [], false, false).';
    fprintf('done.\n');
else
    fprintf('"%s" ... ', params.reference_file);
    if endsWith(params.reference_file, '.sofa')
        tmp = SOFAload(params.reference_file);
        reference_IRs = get_SOFA_IRs(tmp);
        reference_fs = tmp.Data.SamplingRate;
        if size(tmp.ListenerView, 1) == 360
            reference_azel_deg = tmp.ListenerView(:, 1:2).';
        else
            reference_azel_deg = tmp.SourcePosition(:, 1:2).';
        end
        clear tmp;
        fprintf('done.\n');

    elseif endsWith(params.reference_file, '.wav')
        [reference_IRs, reference_fs] = audioread(params.reference_file);
        assert(size(reference_IRs, 2) == 720, 'Unknown data format.');
        % Decompose SSR BRIRs
        reference_IRs = reshape(reference_IRs, ... 
            size(reference_IRs, 1), 2, size(reference_IRs, 2)/2);
        reference_IRs = permute(reference_IRs, [1, 3, 2]);
        if contains(params.reference_file, '_vertical', 'IgnoreCase', true)
            reference_azel_deg = get_array_grid('SSR_vertical', [], false, false).';
        else
            reference_azel_deg = -get_array_grid('SSR', [], false, false).';
        end
        fprintf('done.\n');
    else
        error('Unknown file type.');
    end
    assert(reference_fs == fs, 'Mismatch in sampling frequencies.');
    clear reference_fs;
end

fprintf('Transforming head orientations ... ');
fprintf('to %d target directions ... ', size(reference_azel_deg, 2));
% Correct small values to 0
reference_azel_deg(abs(reference_azel_deg) < 1e-5) = 0;

% Determine if vertical rotation grid is used
is_vertical = all(reference_azel_deg(1, :) == 0 ...
    | reference_azel_deg(1, :) == 180 | reference_azel_deg(1, :) == -180);

if is_vertical
    % Translate SH spherical coordinate convention into continous elevation
    % rotation (azimuth 0 deg, elevation 0deg to -359 deg).
    % For anechoic data this could be skipped here and instead done in the 
    % plotting functions for a depcition over a 360 degree rotation.
    % For room data skipping this will result in different ear signals 
    % i.e., with a sudden orientation swith at +-90 degrees.
    idx = reference_azel_deg(1, :) == 180 | reference_azel_deg(1, :) == -180;
    reference_azel_deg(2, idx) = -reference_azel_deg(2, idx) - 180;
    reference_azel_deg(1, idx) = 0;
    idx = reference_azel_deg(2, :) > 0;
    reference_azel_deg(2, idx) = reference_azel_deg(2, idx) - 360;
    clear idx;
end

% Change to "navigational" azimuth and elevation grid, since it looks nicer
% (done separately for elevation, since it does not work as intended
% otherwise)
reference_azel_deg(1, :) = sph2nav(reference_azel_deg(1, :));
reference_azel_deg(2, :) = sph2nav(reference_azel_deg(2, :));
% reference_azel_deg(reference_azel_deg == 180) = -180; % for nicer plotting

% Generate matrix of independently rendered head orientations
output_azel_rad = repmat(deg2rad(reference_azel_deg), [1, 1, size(array_signals, 1)]); % [2; 360; samples]

% Change to single precision since it is sufficient
reference_azel_deg = single(reference_azel_deg);
output_azel_rad = single(output_azel_rad);

if is_vertical
    params.azim_align = false;
    array_file = [array_file, '_vertical'];
end
fprintf('done.\n');
fprintf('\n');

%% precompute HRTFs with or without MagLS/eMagLS
if is_anec_simulation
    % Load high-resolution HRTFs for full-spherical diffuse-field evaluation
    [~, hrtfs_ref_nm] = render_sh_coefficients_binaurally(params.hrir_file, ...
        zeros(params.half_block_length+1, params.num_harmonics), ...
        0, params.hrir_max_N, params.sh_type, params.sh_with_weights);
end

fig_name = array_file;

if any(strcmpi(params.eq_type, {'', 'eBreve'}))
    % ---------------------- raw HRTF SH coefficients ---------------------
    [~, hrtfs_nm] = render_sh_coefficients_binaurally(params.hrir_file, ...
        zeros(params.half_block_length+1, params.num_harmonics), ...
        0, params.N, params.sh_type, params.sh_with_weights);
    % introduce delay to have all resulting BRIRs time aligned
    hrtfs_nm = apply_delay_fd(hrtfs_nm, params.eq_magls_len/2);

elseif contains(params.eq_type, 'MagLS', 'IgnoreCase', true)
    fprintf('Loading HRIR file "%s" ... ', params.hrir_file);
    hrir_sofa = SOFAload(params.hrir_file);
    % extract measurement directions
    hrir_grid_rad = azEl2azCo(deg2rad(hrir_sofa.SourcePosition(:, 1:2)).'); % elevation to colatitude
    % extract IRs
    hrir = get_SOFA_IRs(hrir_sofa);
    fprintf('with %d samples ... ', size(hrir, 1));
    clear hrir_sofa;
    fprintf('done.\n');

    assert(params.eq_magls_len <= params.half_block_length * 2, ...
        'Block length is shorter than recommended length of %d samples for MagLS filters.', ...
        params.eq_magls_len);
    
    fprintf('Computing %s filters ... at SH order %d ... with %d samples ... ', ...
        params.eq_type, params.N, params.eq_magls_len);
    is_diffusenessConst = contains(params.eq_type, 'wDC', 'IgnoreCase', true);
    if is_eMagLS
        if strcmpi(params.array_type, 'SMA')
            [hrirs_nm(:, :, 1), hrirs_nm(:, :, 2)] = getEMagLsFilters( ...
                hrir(:, :, 1), hrir(:, :, 2), hrir_grid_rad(1, :).', hrir_grid_rad(2, :).', ...
                array_radius, array_grid_rad(1, :).', array_grid_rad(2, :).', params.N, fs, ...
                params.eq_magls_len, is_diffusenessConst, params.sh_type, @getSH_AmbiEnc);
        else % EMA
            if contains(params.eq_type, 'eMagLSinCH', 'IgnoreCase', true)
                [hrirs_nm(:, :, 1), hrirs_nm(:, :, 2)] = getEMagLsFiltersEMAinCH( ...
                    hrir(:, :, 1), hrir(:, :, 2), hrir_grid_rad(1, :).', hrir_grid_rad(2, :).', ...
                    array_radius, array_grid_rad(1, :).', params.N, fs, ...
                    params.eq_magls_len, is_diffusenessConst, params.sh_type, @getSH_AmbiEnc, @getCH);
            else % in SH
                [hrirs_nm(:, :, 1), hrirs_nm(:, :, 2)] = getEMagLsFiltersEMAinSH( ...
                    hrir(:, :, 1), hrir(:, :, 2), hrir_grid_rad(1, :).', hrir_grid_rad(2, :).', ...
                    array_radius, array_grid_rad(1, :).', params.N, fs, ...
                    params.eq_magls_len, is_diffusenessConst, params.sh_type, @getSH_AmbiEnc);
            end
%             [hrirs_nm2(:, :, 1), hrirs_nm2(:, :, 2)] = getEMagLsFiltersEMAinCHtoSH( ...
%                 hrir(:, :, 1), hrir(:, :, 2), hrir_grid_rad(1, :).', hrir_grid_rad(2, :).', ...
%                 array_radius, array_grid_rad(1, :).', params.N, fs, ...
%                 params.eq_magls_len, false, params.sh_type, @getSH_AmbiEnc, @getCH);
%             [hrirs_nm3(:, :, 1), hrirs_nm3(:, :, 2)] = getEMagLsFiltersEMAinEHtoSH( ...
%                 hrir(:, :, 1), hrir(:, :, 2), hrir_grid_rad(1, :).', hrir_grid_rad(2, :).', ...
%                 array_radius, array_grid_rad(1, :).', params.N, fs, ...
%                 params.eq_magls_len, false, params.sh_type, @getSH_AmbiEnc, @getCH);
%             [hrirs_nm3(:, :, 1), hrirs_nm3(:, :, 2)] = getEMagLsFiltersEMA_old( ...
%                 hrir(:, :, 1), hrir(:, :, 2), hrir_grid_rad(1, :).', hrir_grid_rad(2, :).', ...
%                 array_radius, array_grid_rad(1, :).', params.N, fs, ...
%                 params.eq_magls_len, false, params.sh_type, @getSH_AmbiEnc, @getCH);
%             [hrirs_nmT(:, :, 1), hrirs_nmT(:, :, 2)] = getEMagLsFiltersEMA_Tommi( ...
%                 hrir(:, :, 1), hrir(:, :, 2), hrir_grid_rad(1, :).', hrir_grid_rad(2, :).', ...
%                 array_radius, array_grid_rad(1, :).', params.N, fs, ...
%                 params.eq_magls_len);
        end
    else % MagLS
        [hrirs_nm(:, :, 1), hrirs_nm(:, :, 2)] = getMagLsFilters( ...
            hrir(:, :, 1), hrir(:, :, 2), hrir_grid_rad(1, :).', hrir_grid_rad(2, :).', ...
            params.N, fs, params.eq_magls_len, is_diffusenessConst, params.sh_type, @getSH_AmbiEnc);
        if contains(params.eq_type, '+', 'IgnoreCase', true)
            fprintf('done.\n');
            fprintf('Computing array diffuse-field compensation filter ... ');
            fprintf('at SH order %d ... with %d samples ... ', params.N, params.eq_magls_len);
            if contains(params.eq_type, 'SHF', 'IgnoreCase', true)
                hrirs_adf = getMagLsSphericalHeadFilter(array_radius, params.N, fs, params.eq_magls_len);
            elseif contains(params.eq_type, 'ADF', 'IgnoreCase', true)
                hrirs_adf = getMagLsArrayDiffuseFilter( ...
                    array_radius, array_grid_rad(1, :).', array_grid_rad(2, :).', ...
                    params.N, fs, params.eq_magls_len, params.sh_type, @getSH_AmbiEnc);
            else
                error('Unkown equalization type "%s".', params.eq_type);
            end
        end
    end
    clear hrir_grid_rad is_diffusenessConst;

    if strcmpi(params.sh_type, 'complex')
        % Apply correction of wrong rotation. Not sure why this is 
        % neccesary, maybe its a mistake in the (e)MagLS implementation?
        hrirs_nm = conj(hrirs_nm);
    end
    fprintf('done.\n');

    if params.eq_magls_len < params.half_block_length * 2
        fprintf('Zero-padding %s filters ... to %d samples ... ', ...
            params.eq_type, params.half_block_length * 2);
        hrirs_nm(end+1:params.half_block_length * 2, :, :) = 0;
        fprintf('done.\n');
    end
    if contains(params.eq_type, 'MagLS+', 'IgnoreCase', true)
        fprintf('Applying %s filter ... ', params.eq_type);
        % Apply zero-padding before for convenience
        hrirs_adf(end+1:size(hrirs_nm, 1), :, :) = 0;
        for e = 1 : size(hrirs_nm, 3)
            hrirs_nm(:, :, e) = fftfilt(hrirs_adf, hrirs_nm(:, :, e));
        end; clear e;
        % Remove filter delay
        hrirs_nm = applySubsampleDelay(hrirs_nm, -params.eq_magls_len/2);
        fprintf('done.\n');

        if DO_PLOT
            fprintf('Generating plot "%s" ... ', fig_name);
            if params.plot_fig_size
                fig = AKf(params.plot_fig_size(1), params.plot_fig_size(2));
            else
                fig = AKf();
            end
            set(fig, 'NumberTitle', 'Off', 'Name', fig_name, 'Renderer', 'painters');
            tl = tiledlayout(fig, 2, 1, 'TileSpacing', 'tight', 'Padding', 'tight');
            drawnow;
            nexttile(tl);
            AKp(hrirs_adf, 'et2d', 'xu', 'n', 'dr', 110);
            legend(sprintf('%s filter', params.eq_type), 'Location', 'East');
            text(0.005, 0.99, {sprintf('Array type "%s"', params.array_type); ...
                sprintf('Array radius %.2f cm', array_radius * 100); ...
                sprintf('Array sensors %d', size(array_grid_rad, 2))}, ...
                'Color', 'b', 'FontSize', 14, 'VerticalAlignment', 'Top', 'Units', 'normalized');
            drawnow;
            nexttile(tl);
            AKp(hrirs_adf, 'm2d', 'fs', fs, 'lw', 2, 'x', [10, fs/2]);
            legend(sprintf('%s filter', params.eq_type), 'Location', 'South');
            drawnow;
            if DO_EXPORT_PLOT
                fprintf('exporting ... ');
                [~, ~] = mkdir(plot_export_dir); % ignore warning if directory already exists
                file_name = fullfile(plot_export_dir, [fig_name, '.pdf']);
                exportgraphics(fig, file_name, 'Append', ~is_first_plot);
                is_first_plot = false;
            end
            if is_batch_mode; close(fig); end
            clear irs fig tl;
            fprintf('done.\n');
        else
            disp('Generating plot ... skipped.');
        end
    end
    % fft and remove redundant information
    hrtfs_nm = AKboth2singleSidedSpectrum(fft(hrirs_nm));
    hrtfs_nm = apply_delay_fd(hrtfs_nm, 0.21 * params.eq_magls_len);
    if is_eMagLS
        hrtfs_nm = apply_delay_fd(hrtfs_nm, params.rf_len/2);
    end
    clear hrir hrirs_nm;
else
    error('Unkown equalization type "%s".', params.eq_type);
end

%% load XMA filters
fprintf('Parsing XMA filter names ... ');
if strcmpi(params.array_type, 'XMA')
    base_name = regexpi(params.array_file, '.+\/.+\_(XMA.+?)(\_.*)?\.', 'tokens');
    base_name = fullfile(fileparts(params.array_file), ...
        sprintf('Anechoic_%s_SrcEar', base_name{1}{1}));
    fprintf('done.\n');

    % file_name = sprintf('%s_x_nm_%s_v13.mat', base_name, params.sh_type);
    file_name = sprintf('%s_x_nm_%s_v14.mat', base_name, params.sh_type);
    fprintf('Loading XMA calibration file "%s" ... ', file_name);
    tmp = load(file_name);
    x_nm = tmp.x_nm;
%     if size(x_nm, 2) > size(arrayant _signals, 2)
%         warning('received wrong number of input channels, trying to fix ... ');
%         x_nm = x_nm(:, 1 : size(x_nm, 2) / size(array_signals, 2) : end, :);
%     end
    array_radius = tmp.R_projection; % This has to be set before generating the radial filters!
    assert(tmp.fs == fs, 'Mismatch in sampling frequencies.');
    clear tmp;
    fprintf('done.\n');

    assert(size(x_nm, 1) < params.half_block_length, ...
        ['For XMAs, half_block_length has to be larger than %d ', ...
        'because of the calibration filters.'], size(x_nm, 1));
    if size(x_nm, 1) ~= params.eq_xma_len
        warning(['The XMA calibration filters do not have the provided ', ...
            'length of %d samples. The provided circshift length of the ', ...
            'rendering result will be adjusted by the difference.'], params.eq_xma_len);
        params.circ_len = params.circ_len + (params.eq_xma_len - size(x_nm, 1))/2;
    end

    fprintf('Loading XMA equalization file ');
    if strcmpi(params.eq_type, 'eBreve')
        % file_name = sprintf('%s_e_nm_%s_v13.mat', base_name, params.sh_type);
        file_name = sprintf('%s_e_nm_%s_v14.mat', base_name, params.sh_type);
        fprintf('"%s" ... ', file_name)
        tmp = load(file_name);
        e_breve_nm = cat(3, tmp.e_breve_nm_l, tmp.e_breve_nm_r);
        clear tmp;
        clear file_name;
        fprintf('done.\n');
    
        assert(size(x_nm, 1) + size(e_breve_nm, 1) - 1 < params.half_block_length, ...
            ['For XMAs, half_block_length has to be larger than %d ', ...
            'because of the calibration and equalization filters.'], ...
            size(x_nm, 1) + size(e_breve_nm, 1) - 1);
        if size(e_breve_nm, 1) ~= params.eq_xma_len
            warning(['The XMA equalization filters do not have the provided ', ...
            'length of %d samples. The provided circshift length of the ', ...
            'rendering result will be adjusted by the difference.'], params.eq_xma_len);
            params.circ_len = params.circ_len + (params.eq_xma_len - size(e_breve_nm, 1))/2;
        end
    else
        e_breve_nm = zeros(2, params.num_harmonics, 2);
        e_breve_nm(1, :, :) = 1;
        % Introduce delay to have all resulting BRIRs time aligned
        hrtfs_nm = apply_delay_fd(hrtfs_nm, params.eq_xma_len/2);
        fprintf('... skipped.\n');
    end

    fprintf('Zero-padding XMA filters ... to %d samples ... ', params.half_block_length * 2);
    x_nm(end+1:params.half_block_length*2, :, :) = 0;
    e_breve_nm(end+1:params.half_block_length*2, :, :) = 0;
    X_nm = AKboth2singleSidedSpectrum(fft(x_nm));
    E_breve_nm = AKboth2singleSidedSpectrum(fft(e_breve_nm));
    clear x_nm e_breve_nm;
    fprintf('done.\n');

    fprintf('Applying XMA equalization filters ... ');
    hrtfs_nm = hrtfs_nm .* E_breve_nm;
    clear E_breve_nm;
    fprintf('done.\n');
else
    % introduce delay to have all resulting BRIRs time aligned
    hrtfs_nm = apply_delay_fd(hrtfs_nm, params.eq_xma_len);
    fprintf('... skipped.\n');
end
fprintf('\n');

%% precompute radial filters
fprintf('Computing radial filters ... ');
C = compute_speed_of_sound(array_temp);
% Alternatives to calculate the SMA aliasing frequency could be e.g.
% https://github.com/polarch/Spherical-Array-Processing/blob/master/sphArrayAliasLim.m
freq_SMA_alias = C * params.N / (2*pi * array_radius);
if DO_PLOT_PRESEN
    % Adjust plotting lower frequency limit based on alias frequency
    params.plot_freq_rng(1) = floor(freq_SMA_alias / 4 / 100) * 100; % in Hz
end

if is_eMagLS
    one_over_b_n_limited = [];
    fprintf('skipped.\n');
else
    fprintf('at SH order %d ... with %d samples ... ', params.N, params.rf_len);
    if isempty(array_temp)
        fprintf('with constant speed of sound of %.1f m/s ... ', C);
    else
        fprintf('with calculated speed of sound of %.1f m/s ... ', C);
    end
    freqs = linspace(0, fs/2, params.rf_len/2 + 1).';
    k = 2*pi * freqs / C;
    assert(size(k, 1) < params.half_block_length)

    if any(strcmpi(params.array_type, {'SMA', 'XMA'}))
        assert(params.half_block_length >= 1024, ...
            ['For SMAs, half_block_length has to be larger than %d ', ...
            'because of the radial filters.'], params.rf_len);
        [~, one_over_b_n_limited, ~] = get_sma_radial_filters(k, array_radius, ...
            params.N, params.rf_lim_dB, params.rf_reg_type, 2);
    else % EMA
        assert(params.half_block_length >= 2048, ...
            ['For EMAs, half_block_length has to be larger than %d ', ...
            'because of the radial filters.'], params.rf_len);
        % The EMA function uses Jens' normalization of the radial filters
        % which is different from the Ambisonics convention
        [one_over_b_n_limited, ~] = get_ema_radial_filters(k, array_radius, ...
            params.N, params.rf_lim_dB + db(4*pi), params.rf_reg_type, 2);
    end
    clear freqs k;
    fprintf('done.\n');
    
    if params.rf_len < params.half_block_length * 2
        fprintf('Zero-padding radial filters ... to %d samples ... ', params.half_block_length * 2);
        one_over_b_n_limited = ifft(AKsingle2bothSidedSpectrum(one_over_b_n_limited));
        one_over_b_n_limited(end+1:params.half_block_length * 2, :, :) = 0;
        one_over_b_n_limited = AKboth2singleSidedSpectrum(fft(one_over_b_n_limited));
        fprintf('done.\n');
    end
end

if DO_PLOT
    fprintf('Generating plot "%s" ... ', fig_name);
    text_str = {sprintf('Array type "%s"', params.array_type); ...
        sprintf('Array radius %.2f cm', array_radius * 100); ...
        sprintf('Environment speed of sound %.1f m/s', C); ...
        sprintf('Alias frequency %.0f Hz', freq_SMA_alias); ...
        sprintf('SH type "%s"', params.sh_type)};
    if is_eMagLS
        text_str = {text_str{1:end-1}, ...
            sprintf('Filter length %d samples', params.eq_magls_len), ...
            text_str{end}};
        if contains(params.eq_type, 'eMagLSinCH', 'IgnoreCase', true)
            lgd_str = get_CH_strings(params.N);
        else
            lgd_str = get_SH_strings(params.N);
        end
        lgd_title = sprintf('%s filters [N,M]', params.eq_type);
        irs = ifft(AKsingle2bothSidedSpectrum(hrtfs_nm));
        if params.N > 8 && ~contains(params.eq_type, 'eMagLSinCH', 'IgnoreCase', true)
            % Limit to coefficients of lower orders
            % (otherwise plot is illegible)
            irs = irs(:, 1:(8+1)^2, :);
        end
    else
        text_str = {text_str{1:end-1}, ...
            sprintf('Filter length %d samples', params.rf_len), ...
            sprintf('Filter limit %d dB', params.rf_lim_dB), ...
            sprintf('Filter regularization "%s"', params.rf_reg_type), ...
            text_str{end}};
        lgd_str = arrayfun(@num2str, 0:params.N, 'UniformOutput', false);
        lgd_title = 'Radial filters, N';
        irs = ifft(AKsingle2bothSidedSpectrum(one_over_b_n_limited));
        if strcmpi(params.array_type, 'EMA')
            % Do not plot symmetric part of EMA filters
            irs = irs(:, params.N+1:end);
        end
    end

    if params.plot_fig_size
        fig = AKf(params.plot_fig_size(1), params.plot_fig_size(2));
    else
        fig = AKf();
    end
    set(fig, 'NumberTitle', 'Off', 'Name', fig_name);
    tl = tiledlayout(fig, 2, is_eMagLS+1, 'TileIndexing', 'columnmajor', ...
        'TileSpacing', 'tight', 'Padding', 'tight');
    drawnow;
    for e = 1 : is_eMagLS+1 % each ear
        ax = nexttile(tl);
        AKp(irs(:, :, e), 'et2d', 'xu', 'n', 'dr', 110);
        if ~ismatrix(irs)
            if e == 1
                title(ax, [ax.Title.String, '   (left ear)']);
            else
                title(ax, [ax.Title.String, '   (right ear)']);
            end
        end
        lgd = legend(lgd_str, 'FontName', 'FixedWidth', ...
            'Location', 'East', 'NumColumns', ceil(size(irs, 2) / 21));
        lgd.ItemTokenSize(1) = 15; % decrease line length
        title(lgd, lgd_title);
        text(0.005, 0.99, text_str, ...
            'Color', 'b', 'FontSize', 14, 'VerticalAlignment', 'Top', 'Units', 'normalized');
        drawnow;

        ax = nexttile(tl);
        AKp(irs(:, :, e), 'm2d', 'fs', fs, 'lw', 2, 'x', [10, fs/2], 'dr', params.rf_lim_dB + 40);
        if ~ismatrix(irs)
            if e == 1
                title(ax, [ax.Title.String, '   (left ear)']);
            else
                title(ax, [ax.Title.String, '   (right ear)']);
            end
        end
        drawnow;
    end; clear e;
    if DO_EXPORT_PLOT
        fprintf('exporting ... ');
        [~, ~] = mkdir(plot_export_dir); % ignore warning if directory already exists
        file_name = fullfile(plot_export_dir, [fig_name, '.pdf']);
        exportgraphics(fig, file_name, 'Append', ~is_first_plot);
        is_first_plot = false;
    end
    if is_batch_mode; close(fig); end
    clear fig tl ax lgd irs text_str lgd_str lgd_title;
    fprintf('done.\n');
else
    disp('Generating plot ... skipped.');
end
fprintf('\n');

%% render
fprintf('Pre-computing SMA basis functions ... ');
if strcmpi(params.array_type, 'SMA')
    fprintf('at SH order %d ... ', params.N);
    S_Ynm_inv = get_basis_inv_sh_coeffs_sma(params.N, array_grid_rad, ...
        array_grid_weights, params.sh_type);
    fprintf('done.\n');
else
    fprintf('skipped.\n');
end

% loop over input signal
output_signals = zeros(size(array_signals, 1), size(output_azel_rad, 2), 2);
S_signals = zeros(size(array_signals, 1), size(hrtfs_nm, 2));
for s = 1 : params.half_block_length : size(array_signals, 1)
    fprintf('Computing block %d of %d ... ', ...
        (s+params.half_block_length-1) / params.half_block_length, ...
        ceil(size(array_signals, 1) / params.half_block_length));

    % get input block
    if s == 1
        input_block = [zeros(params.half_block_length, size(array_signals, 2)); ...
            array_signals(1 : params.half_block_length, :)];
    else
        input_block = array_signals(s - params.half_block_length : s + params.half_block_length - 1, :);
    end
    input_block_spec = AKboth2singleSidedSpectrum(fft(input_block));
    
    % get sound field coefficients
    if strcmpi(params.array_type, 'SMA')
        if exist('S_Ynm_inv', 'var')
            % With pre-computed inverted SH basis functions (this should be
            % faster, particularly for high SH orders)
            S_breve = get_sound_field_sh_coeffs_sma(input_block_spec, ...
                one_over_b_n_limited, S_Ynm_inv);
        else
            % Without pre-computed inverted SH basis functions
            S_breve = get_sound_field_sh_coeffs_sma(input_block_spec, ...
                one_over_b_n_limited, params.N, array_grid_rad, ...
                array_grid_weights, params.sh_type);
        end
    elseif strcmpi(params.array_type, 'EMA')
        S_breve = get_sound_field_sh_coeffs_ema(input_block_spec, ...
            one_over_b_n_limited, params.N, array_grid_rad(1, :), ...
            params.sh_type, params.is_ch);
    elseif strcmpi(params.array_type, 'XMA')
        S_breve = get_sound_field_sh_coeffs_xma(input_block_spec, ...
            one_over_b_n_limited, params.N, X_nm);
    else
        error('Unknown array type "%s".', params.array_type);
    end

    % Nyquist bin
    S_breve(end, :) = real(S_breve(end, :));

    tic; % start measuring execution time
    output_block = render_sh_coefficients_binaurally_blockwise( ...
        hrtfs_nm, S_breve, double(output_azel_rad(:, :, s)), params.N, params.sh_type);
    fprintf('in %.0fm %.1fs ... ', mod(toc,3600)/60, mod(toc,60));

    % overlap save
    output_signals(s:s+params.half_block_length-1, :, :) = ...
        output_block(params.half_block_length+1:end, :, :);
    if is_anec_simulation
        S_block = ifft(AKsingle2bothSidedSpectrum(S_breve));
        S_signals(s:s+params.half_block_length-1, :, :) = ...
            S_block(params.half_block_length+1:end, :, :);
    end
    
    fprintf('done.\n');
end
clear s input_block input_block_spec S_breve output_block S_block;

if is_batch_mode
    % clear variables for lower RAM usage
    % (done for consecutive steps as well)
    clear array_signals output_azel_rad;
end
fprintf('\n');

%% validate diffuse-field response
fprintf('Computing anechoic array simulation ... ');
if is_anec_simulation
    is_EMA_sigs = strcmpi(params.array_type, 'EMA') && ~is_eMagLS;

    % Simulate plane wave impinging on array
    % (identical to `getEMagLsFilters()`)
    sim.returnRawMicSigs = is_EMA_sigs; % raw mic signals, no SHs!
    sim.order = params.N;
    sim.fs = fs;
    sim.irLen = params.eq_magls_len;
    sim.oversamplingFactor = 1;
    sim.simulateAliasing = true;
    sim.radialFilter = 'none';
    sim.smaRadius = array_radius;
    sim.smaDesignAziZenRad = array_grid_rad.';
    sim.waveModel = 'planeWave';
    sim.arrayType = 'rigid';
    sim.shDefinition = params.sh_type;
    sim.shFunction = @getSH_AmbiEnc;
    fprintf('with @%s("%s") ... ', func2str(sim.shFunction), sim.shDefinition);
    smairMat = getSMAIRMatrix(sim);
    simulationOrder = sqrt(size(smairMat, 2)) - 1;
    fprintf('at SH order %d ... ', simulationOrder);
    fprintf('done.\n');
    
    fprintf('Computing diffuse-field responses ... ');
    hrtfs_ref_rms = compute_spectral_average_SH(hrtfs_ref_nm, 2, false);
    smair_df = compute_spectral_average_SH(smairMat, 2, true).';
    if strcmpi(params.array_type, 'EMA')
        % Compensate magnitude for EMA
        smair_df = smair_df * sqrt(2*simulationOrder + 1) / sqrt(size(smairMat, 2));
    end
    clear sim smairMat simulationOrder;
    % Zero-padding to target length
    smair_df = apply_delay_fd(smair_df, size(smair_df, 1) / 2); % make causal
    smair_df_ir = ifft(smair_df, 'symmetric'); % symmetric is required
    smair_df_ir(end+1:params.half_block_length*2, :, :) = 0;
    smair_df = AKboth2singleSidedSpectrum(fft(smair_df_ir));
    % Apply radial filters
    if is_EMA_sigs
        % For EMA without eMagLS
        smair_df_ir = get_sound_field_sh_coeffs_from_ema_t( ...
            smair_df_ir, ifft(AKsingle2bothSidedSpectrum(one_over_b_n_limited)), ...
            params.N, array_grid_rad(1, :));
        smair_df = AKboth2singleSidedSpectrum(fft(smair_df_ir));
    elseif ~is_eMagLS
        % For SMA without eMagLS
        smair_df = smair_df .* sh_repToOrder(one_over_b_n_limited.').';
    end
    clear is_EMA_sigs smair_df_ir;
    % Combine with rendering HRTFs
    ears_df = hrtfs_nm .* smair_df;
    ears_df_rms = compute_spectral_average_SH(ears_df, 2, true);
    if strcmpi(params.array_type, 'EMA')
        % Match HRTF magnitude for EMA
        ears_df_rms = ears_df_rms / 2;
        if is_eMagLS
            ears_df_rms = ears_df_rms / sqrt(2) / (4*pi) * (params.N+1);
        end
    else
        % Match HRTF magnitude for SMA
        ears_df_rms = ears_df_rms * sqrt(2) / (4*pi);
    end
    ears_df_rms = ears_df_rms / sqrt(size(hrtfs_ref_nm, 2));
    clear smair_df ears_df hrtfs_ref_nm;

    % Truncating to target length
    S_signals = S_signals(1:params.half_block_length*2, :);
    S_breve = AKboth2singleSidedSpectrum(fft(S_signals));
    % Zero-padding to rendering length
    hrtfs_nm = ifft(AKsingle2bothSidedSpectrum(hrtfs_nm));
    hrtfs_nm(end+1:size(S_signals, 1), :, :) = 0;
    hrtfs_nm = AKboth2singleSidedSpectrum(fft(hrtfs_nm));
    % Combine with rendering HRTFs
    ears_S = hrtfs_nm .* S_breve;
    % Calculate spherical average of anechoic array response
    ears_S_rms = compute_spectral_average_SH(ears_S, 2, false);
    ears_S_amp = rms(ears_df_rms, 'all') / rms(ears_S_rms, 'all');
    ears_S_rms = ears_S_rms * ears_S_amp;
    clear S_signals S_breve ears_S;
    fprintf('done.\n');

    if DO_PLOT
        fprintf('Generating plot "%s" ... ', fig_name);
        if params.plot_fig_size
            fig = AKf(params.plot_fig_size(1), params.plot_fig_size(2));
        else
            fig = AKf();
        end
        set(fig, 'NumberTitle', 'Off', 'Name', fig_name, 'Renderer', 'painters');
        tl = tiledlayout(fig, 2, size(hrtfs_ref_rms, 2), 'TileIndexing', 'columnmajor', ...
            'TileSpacing', 'tight', 'Padding', 'tight');
        drawnow;
        lgd_str = {sprintf('HRTF reference (N=%d)', params.hrir_max_N), ...
            sprintf('Rendered simulated array (N=%d)', params.N), ...
            sprintf('Rendered anechoic array (N=%d) (adjusted by %+.1f dB)', params.N, db(ears_S_amp))};

        for e = 1 : size(hrtfs_ref_rms, 2) % each ear
            nexttile(tl);
            axis off;

            ax = nexttile(tl);
            pd = AKp(ifft(AKsingle2bothSidedSpectrum( ...
                [hrtfs_ref_rms(:, e), ears_df_rms(:, e), ears_S_rms(:, e)])), ...
                'm2d', 'fs', fs, 'x', [10, fs/2], 'dr', 40, 'lw', 2, 'c', 'k');
            pd.h(1).LineWidth = 5;
            pd.h(1).Color = [.7, .7, .7];
            pd.h(3).Color = 'b';
            if e == 1
                title(ax, [ax.Title.String, '   (left ear)']);
            else
                title(ax, [ax.Title.String, '   (right ear)']);
            end
            lgd = legend(lgd_str, 'Location', 'South');
            title(lgd, 'Diffuse field response');
            drawnow;
        end; clear e;

        if DO_EXPORT_PLOT
            fprintf('exporting ... ');
            [~, ~] = mkdir(plot_export_dir); % ignore warning if directory already exists
            file_name = fullfile(plot_export_dir, [fig_name, '.pdf']);
            exportgraphics(fig, file_name, 'Append', ~is_first_plot);
            is_first_plot = false;
        end
        if is_batch_mode; close(fig); end
        clear fig tl ax lgd lgd_str pd;
        fprintf('done.\n');
    else
        disp('Generating plot ... skipped.');
    end
    if is_batch_mode
        clear hrtfs_ref_rms ears_df_rms ears_S_amp ears_S_rms;
    end
else
    fprintf('skipped.\n');
end
fprintf('\n');

%% align and validate data against reference
fprintf('Time-shifting BRIRs ... ');
if params.circ_len
    params.circ_len = round(params.circ_len);
    fprintf('by %d samples ... ', params.circ_len);
    output_signals = circshift(output_signals, params.circ_len, 1);
    if params.trunc_win_len
        % Apply fade-in window
        trunc_win = cos(linspace(pi/2, 0, params.trunc_win_len).').^2;
        output_signals(1 : params.trunc_win_len, :, :) = trunc_win ...
            .* output_signals(1 : params.trunc_win_len, :, :);
        clear trunc_win;
    end
    fprintf('done.\n');
else
    fprintf('skipped.\n');
end

if ~isempty(params.reference_file)
    fprintf('Zero-padding reference BRIRs ... ');
    ir_len = max([size(output_signals, 1), size(reference_IRs, 1)]);
    if ir_len > size(reference_IRs, 1)
        fprintf('to %d samples ... ', ir_len);
        reference_IRs(end+1:ir_len, :) = 0;
        fprintf('done.\n');
    else
        fprintf('skipped.\n');
    end
    fprintf('Zero-padding result BRIRs ... ');
    if ir_len > size(output_signals, 1)
        fprintf('to %d samples ... ', ir_len);
        output_signals(end+1:ir_len, :) = 0;
        fprintf('done.\n');
    else
        fprintf('skipped.\n');
    end
    clear ir_len;

    fprintf('Computing azimuth alignment ... ');
    if params.azim_align
        % calculate frequency range based on potentially relevant (>500 Hz) 
        % and valid (mirophone array aliasing frequency, <20 kHz) ILD ceues
        azim_align_freq_rng = [500, min([20e3, freq_SMA_alias])];
        [azim_min_shift, azim_min_shift_ind, azim_specs_diff, azim_specs_diff_mean] = ...
            compute_azim_alignment(double(reference_azel_deg), reference_IRs, ...
            output_signals, fs, params.azim_align_frac_sm, azim_align_freq_rng);
        fprintf(['found mismatch minima at [%.1f, %.1f] deg ', ...
            'for left / right ear individually ... '], azim_min_shift_ind);
        fprintf('done.\n');
        
        fprintf('Aligning results BRIRs azimuth ... ');
        fprintf('by %.1f deg for combined minimum at both ears ... ', azim_min_shift);
        output_signals = circshift(output_signals, -azim_min_shift, 2);
        fprintf('done.\n');

        if DO_PLOT
            fprintf('Generating plot "%s" ... ', fig_name);
            if params.plot_fig_size
                fig = AKf(params.plot_fig_size(1), params.plot_fig_size(2));
            else
                fig = AKf();
            end
            set(fig, 'NumberTitle', 'Off', 'Name', fig_name);
            title_str = {'Differences averaged over incidence directions'};
            if DO_PLOT_PRESEN; title_str = [title_str(:)', {''}]; end
            if params.azim_align_frac_sm
                title_str{end} = sprintf('%s   (1/%d oct. smoothing)', ...
                    title_str{end}, params.azim_align_frac_sm);
            end
            plot_azim_alignment(reference_azel_deg, ...
                azim_min_shift, azim_min_shift_ind, azim_specs_diff, azim_specs_diff_mean, ...
                fs, title_str, azim_align_freq_rng, true);
            if DO_EXPORT_PLOT
                fprintf('exporting ... ');
                [~, ~] = mkdir(plot_export_dir); % ignore warning if directory already exists
                file_name = fullfile(plot_export_dir, [fig_name, '.pdf']);
                exportgraphics(fig, file_name, 'Append', ~is_first_plot, 'Resolution', params.plot_resolution);
                is_first_plot = false;
            end
            if is_batch_mode; close(fig); end
            clear fig file_name title_str;
            fprintf('done.\n');
        else
            disp('Generating plot ... skipped.');
        end
        clear azim_specs_diff azim_specs_diff_mean;
    else
        disp('skipped.');
    end

    fprintf('Time-aligning reference BRIRs ... ');
    [~, reference_peaks] = max(abs(reference_IRs), [], 1);
    [~, output_peaks] = max(abs(output_signals), [], 1);
    reference_shift = round(median(output_peaks, 'all') - median(reference_peaks, 'all'));
    fprintf('by %d samples of median BRIR peak ... ', reference_shift);
    reference_IRs = circshift(reference_IRs, reference_shift);
    clear reference_peaks output_peaks reference_shift;
    fprintf('done.\n');

    fprintf('Adjusting result BRIR amplitudes ... ');
    % norm_fact = mean(rssq(reference_IRs) ./ rssq(output_signals), 'all');
    % fprintf('by %.1f dB of reference mean RSSQ amplitude ... ', db(norm_fact));
    % output_signals = output_signals * norm_fact;
    % clear norm_fact;
    if is_vertical
        [~, front_idx] = min(abs(reference_azel_deg(2, :)));
    else
        [~, front_idx] = min(abs(reference_azel_deg(1, :)));
    end
    norm_fact = mean(rssq(reference_IRs(:, front_idx, :)) ...
        ./ rssq(output_signals(:, front_idx, :)), 'all');
    fprintf('by %.1f dB of reference frontal RSSQ amplitude ... ', db(norm_fact));
    output_signals = output_signals * norm_fact;
    clear front_idx norm_fact;
    fprintf('done.\n');

    fprintf('Adjusting reference and result BRIR amplitudes ... ');
    norm_fact = max(abs([reference_IRs; output_signals]), [], 'all');
    if norm_fact <= 1
        fprintf('skipped (no clipping detected).\n');
    else
        norm_fact = norm_fact / .99;
        fprintf('by %.1f dB of peak amplitude ... ', -db(norm_fact));
        reference_IRs = reference_IRs / norm_fact;
        output_signals = output_signals / norm_fact;
        fprintf('done.\n');
    end
    clear norm_fact;
end

fprintf('Truncating BRIRs ... ');
if ~isempty(params.trunc_threshold)
    fprintf('at %-.1f dB to peak ... ', params.trunc_threshold);
    etcs_max = db(max(abs(output_signals), [], [2, 3])); % of all channels
    trunc_level = max(etcs_max) + params.trunc_threshold;
    ir_len = find(etcs_max >= trunc_level, 1, 'last');
    if ~isempty(params.reference_file)
        etcs_max = db(max(abs(reference_IRs), [], [2, 3])); % of all channels
        trunc_level = max(etcs_max) + params.trunc_threshold;
        ir_len = max([ir_len, find(etcs_max >= trunc_level, 1, 'last')]);
    end
    if params.trunc_win_len
        if ir_len + params.trunc_win_len > size(output_signals, 1)
            ir_len = size(output_signals, 1);
            params.trunc_win_len = 0;
        else
            ir_len = ir_len + params.trunc_win_len;
        end
    end
    fprintf('to %d samples ... ', ir_len);
    if ~isempty(params.reference_file)
        reference_IRs = reference_IRs(1:ir_len, :, :);
    end
    output_signals = output_signals(1:ir_len, :, :);
    if params.trunc_win_len
        % Apply fade-out window
        trunc_win = cos(linspace(0, pi/2, params.trunc_win_len).').^2;
        output_signals(end - params.trunc_win_len + 1 : end, :, :) = trunc_win ...
            .* output_signals(end - params.trunc_win_len + 1 : end, :, :);
    end
    clear etcs_max trunc_level ir_len trunc_win;
    fprintf('done.\n');
else
    fprintf('skipped.\n');
end

if DO_EXPORT_META
    file_name = fullfile(data_export_dir, [array_file, '_meta.mat']);
    fprintf('Exporting file "%s" ... ', file_name);
    [~, ~] = mkdir(data_export_dir); % ignore warning if directory already exists
    file_vars = whos();
    file_vars = file_vars(~strcmp({file_vars.class}, 'logical'));
    file_vars = file_vars(~strcmp({file_vars.name}, 'array_file'));
    file_vars = file_vars(~strcmp({file_vars.name}, 'array_signals'));
    file_vars = file_vars(~strcmp({file_vars.name}, 'base_name'));
    file_vars = file_vars(~strcmp({file_vars.name}, 'data_export_dir'));
    file_vars = file_vars(~strcmp({file_vars.name}, 'fig_name'));
    file_vars = file_vars(~strcmp({file_vars.name}, 'file_name'));
    file_vars = file_vars(~strcmp({file_vars.name}, 'hrir_file'));
    file_vars = file_vars(~strcmp({file_vars.name}, 'hrtfs_nm'));
    file_vars = file_vars(~strcmp({file_vars.name}, 'one_over_b_n_limited'));
    file_vars = file_vars(~strcmp({file_vars.name}, 'output_signals'));
    file_vars = file_vars(~strcmp({file_vars.name}, 'plot_export_dir'));
    file_vars = file_vars(~strcmp({file_vars.name}, 'reference_file'));
    file_vars = file_vars(~strcmp({file_vars.name}, 'reference_IRs'));
    file_vars = file_vars(~strcmp({file_vars.name}, 'S_signals'));
    file_vars = file_vars(~strcmp({file_vars.name}, 'S_Ynm_inv'));
    file_vars = file_vars(~strcmp({file_vars.name}, 'X_nm'));
    save(file_name, file_vars.name);
    clear file_name file_vars;
    fprintf('done.\n');
else
    disp('Exporting array rendering metadata ... skipped.');
end

if DO_EXPORT_FLAC
    file_name = fullfile(data_export_dir, [array_file, '_binaural.flac']);
    fprintf('Exporting file "%s" ... ', file_name);
    [~, ~] = mkdir(data_export_dir); % ignore warning if directory already exists
    file_data = squeeze(output_signals);
    file_data = file_data / max(abs(file_data), [], 'all'); % normalize to 1
    audiowrite(file_name, file_data, fs, 'BitsPerSample', 16);
    clear file_name file_data;
    fprintf('done.\n');
else
    disp('Exporting array rendering FLACs ... skipped.');
end

fprintf('Exporting array SSR BRIRs ');
if DO_EXPORT_SSR
    export_SSR_BRIRs(fullfile(data_export_dir, array_file), ...
        output_signals, fs, reference_azel_deg);
else
    fprintf('... skipped.\n');
end

if ~isempty(params.reference_file)
    fprintf('Exporting reference SSR BRIRs ');
    if DO_EXPORT_SSR
        file = fullfile(data_export_dir, reference_file);
        if strcmp(filesep(), '\')
            file = strrep(file, '\', '/'); % fix path string compare on Windows
        end
        if contains(params.reference_file, file)
            fprintf('... skipped (identical to reference source file).\n');
        else
            export_SSR_BRIRs(file, reference_IRs, fs, reference_azel_deg);
        end
        clear file;
    else
        fprintf('... skipped.\n');
    end
    fprintf('\n');
end

if DO_PLOT
    fprintf('Computing spectra ... ');
    output_spec = AKboth2singleSidedSpectrum(fft(output_signals));
    times = (0 : size(output_signals, 1) - 1) * 1000 / fs; % in ms
    freqs = linspace(0, fs/2, size(output_spec, 1));
    if isempty(params.reference_file)
        reference_IRs = nan(size(output_signals));
        reference_spec = nan(size(output_spec));
    else
        reference_spec = AKboth2singleSidedSpectrum(fft(reference_IRs));
    end
    fprintf('done.\n');

    fprintf('Applying spectral smoothing ... ');
    if params.plot_frac_sm
        fprintf('of 1/%d octaves ... ', params.plot_frac_sm);
        for e = 1 : size(output_spec, 3) % each ear
            output_spec(:, :, e) = AKfractOctSmooth( ...
                output_spec(:, :, e), 'amp', fs, params.plot_frac_sm);
            if ~isempty(params.reference_file)
                reference_spec(:, :, e) = AKfractOctSmooth( ...
                    reference_spec(:, :, e), 'amp', fs, params.plot_frac_sm);
            end
        end; clear e;
        fprintf('done.\n');
    else
        fprintf('skipped.\n');
    end

    if ~isempty(params.reference_file)
        fprintf('Computing spectral differences after smoothing ... ');
        output_diff_spec = compute_spectral_difference(output_spec, reference_spec);
        fprintf('done.\n');

        fprintf('Normalizing spectral differences ... ');
        fprintf('between %.1f Hz and %.1f Hz ... ', params.plot_norm_freq_rng);
        norm_freqs = get_freqs_inclusive(freqs, params.plot_norm_freq_rng);
        norm_lvl = mean(db(abs(output_diff_spec(norm_freqs, :, :))), [1, 2]); % based on median
        fprintf('by [%+.1f, %+.1f] dB for left / right ear individually ... ', -norm_lvl);
        output_diff_spec = output_diff_spec ./ db2mag(norm_lvl);
        clear norm_freqs norm_lvl;
        fprintf('done.\n');

        fprintf('Computing spectral difference percentiles ... ');
        output_diff_spec_perc = compute_spectral_percentiles( ...
            output_diff_spec, 2, params.plot_percentiles);
        fprintf('done.\n');

        fprintf('Computing spectral differences with weighting ... ');
        norm_spec = abs(reference_spec) ./ max(abs(reference_spec), [], [2, 3]);
        output_diff_spec_w = db2mag(db(abs(output_diff_spec)) .* norm_spec);
        output_diff_spec_perc_w = compute_spectral_percentiles( ...
            output_diff_spec_w, 2, params.plot_percentiles);
        clear norm_spec;
        fprintf('done.\n');

        fprintf('Computing mean absolute spectral error ... ');
        [~, output_diff_spec_log_err, freqs_log] = compute_spectral_error(output_diff_spec, 2, freqs);
        [~, output_diff_spec_log_err_w, ~] = compute_spectral_error(output_diff_spec_w, 2, freqs);
        fprintf('between %.1f Hz and %.1f Hz ... ', params.plot_freq_rng);
        freqs_bin = get_freqs_inclusive(freqs_log, params.plot_freq_rng);
        freqs_log = freqs_log(freqs_bin);
        output_diff_spec_log_err = output_diff_spec_log_err(freqs_bin, :);
        output_diff_spec_log_err_w = output_diff_spec_log_err_w(freqs_bin, :);
        clear freqs_bin;
        fprintf('done.\n');
    end

    % if ~isempty(params.reference_file)
    %     fprintf('Computing spectral differences before smoothing ... ');
    %     output_diff_spec = compute_spectral_difference(output_spec, reference_spec);
    %     fprintf('done.\n');
    % 
    %     fprintf('Normalizing spectral differences ... ');
    %     fprintf('between %.1f Hz and %.1f Hz ... ', params.plot_norm_freq_rng);
    %     norm_freqs = get_freqs_inclusive(freqs, params.plot_norm_freq_rng);
    %     norm_lvl = mean(db(abs(output_diff_spec(norm_freqs, :, :))), [1, 2]); % based on median
    %     fprintf('by [%+.1f, %+.1f] dB for left / right ear individually ... ', -norm_lvl);
    %     output_diff_spec = output_diff_spec ./ db2mag(norm_lvl);
    %     clear norm_freqs norm_lvl;
    %     fprintf('done.\n');
    % 
    %     fprintf('Computing spectral difference percentiles ... ');
    %     output_diff_spec_perc = compute_spectral_percentiles( ...
    %         output_diff_spec, 2, params.plot_percentiles);
    %     fprintf('done.\n');
    % 
    %     fprintf('Computing spectral differences with weighting ... ');
    %     norm_spec = abs(reference_spec) ./ max(abs(reference_spec), [], [2, 3]);
    %     output_diff_spec_w = db2mag(db(abs(output_diff_spec)) .* norm_spec);
    %     output_diff_spec_perc_w = compute_spectral_percentiles( ...
    %         output_diff_spec_w, 2, params.plot_percentiles);
    %     clear norm_spec;
    %     fprintf('done.\n');
    % 
    %     fprintf('Computing mean absolute spectral error ... ');
    %     [~, output_diff_spec_log_err, freqs_log] = compute_spectral_error(output_diff_spec, 2, freqs);
    %     [~, output_diff_spec_log_err_w, ~] = compute_spectral_error(output_diff_spec_w, 2, freqs);
    %     fprintf('between %.1f Hz and %.1f Hz ... ', params.plot_freq_rng);
    %     freqs_bin = get_freqs_inclusive(freqs_log, params.plot_freq_rng);
    %     freqs_log = freqs_log(freqs_bin);
    %     output_diff_spec_log_err = output_diff_spec_log_err(freqs_bin, :);
    %     output_diff_spec_log_err_w = output_diff_spec_log_err_w(freqs_bin, :);
    %     clear freqs_bin;
    %     fprintf('done.\n');
    % 
    %     if params.plot_frac_sm
    %         fprintf('Applying spectral smoothing ... of 1/%d octaves ... ', params.plot_frac_sm);
    %         for e = 1 : size(output_spec, 3) % each ear
    %             output_spec(:, :, e) = AKfractOctSmooth( ...
    %                 output_spec(:, :, e), 'amp', fs, params.plot_frac_sm);
    %             reference_spec(:, :, e) = AKfractOctSmooth( ...
    %                 reference_spec(:, :, e), 'amp', fs, params.plot_frac_sm);
    %             output_diff_spec(:, :, e) = AKfractOctSmooth( ...
    %                 output_diff_spec(:, :, e), 'amp', fs, params.plot_frac_sm);
    %             output_diff_spec_perc(:, :, e) = AKfractOctSmooth( ...
    %                 output_diff_spec_perc(:, :, e), 'amp', fs, params.plot_frac_sm);
    %             output_diff_spec_w(:, :, e) = AKfractOctSmooth( ...
    %                 output_diff_spec_w(:, :, e), 'amp', fs, params.plot_frac_sm);
    %             output_diff_spec_perc_w(:, :, e) = AKfractOctSmooth( ...
    %                 output_diff_spec_perc_w(:, :, e), 'amp', fs, params.plot_frac_sm);
    %         end
    %         fprintf('done.\n');
    %     end
    % end

    plot_idx = 1 : size(reference_azel_deg, 2);
    if size(output_signals, 1) / fs > 1.5
        % reduce spatial resolution to 3 degrees for long BRIRs
        % (otherwise the figure may occupy too many system resources,
        % causing the export to fail)
        if is_vertical
            plot_idx = plot_idx(1 : 2 : end);
        else
            plot_idx = plot_idx(1 : 3 : end);
        end
    end

    lbl_config = {'Interpreter', 'Latex', 'LineWidth', 1, ...
        'HandleVisibility', 'off', 'LabelOrientation', 'horizontal', ...
        'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'right'};

    %% generate plot
    fprintf('Generating plot "%s" ... ', fig_name);
    if params.plot_fig_size
        fig = AKf(params.plot_fig_size(1), params.plot_fig_size(2));
    else
        fig = AKf();
    end
    set(fig, 'NumberTitle', 'Off', 'Name', fig_name);
    tl = tiledlayout(fig, 2, size(reference_IRs, 3), 'TileIndexing', 'columnmajor', ...
        'TileSpacing', 'tight', 'Padding', 'tight');
    colormap('hot');
    drawnow;
    
    for e = 1 : size(output_signals, 3) % each ear
        title_str = {reference_file};
        if DO_PLOT_PRESEN; title_str = [title_str(:)', {''}]; end
        nexttile(tl);
        plot_ETC_azims_IRs(reference_azel_deg(:, plot_idx), times, ...
            reference_IRs(:, plot_idx, :), title_str, e);
        xticklabels(''); xlabel('');
        if e == 1; set(get(gca, 'ColorBar'), 'Location', 'WestOutside'); end
        if e == size(reference_IRs, 3); yticklabels(''); ylabel(''); end
        drawnow;

        title_str = {[array_file, '   ']};
        if DO_PLOT_PRESEN; title_str = [title_str(:)', {''}]; end
        title_str{end} = sprintf('%s->   %s', title_str{end}, hrir_file);
        nexttile(tl);
        plot_ETC_azims_IRs(reference_azel_deg(:, plot_idx), times, ...
            output_signals(:, plot_idx, :), title_str, e);
        if e == 1; set(get(gca, 'ColorBar'), 'Location', 'WestOutside'); end
        if e == size(reference_IRs, 3); yticklabels(''); ylabel(''); end
        drawnow;
    end; clear e title_str;
    if is_batch_mode; clear reference_IRs output_signals; end

    if DO_EXPORT_PLOT
        fprintf('exporting ... ');
        [~, ~] = mkdir(plot_export_dir); % ignore warning if directory already exists
        file_name = fullfile(plot_export_dir, [fig_name, '.pdf']);
        exportgraphics(fig, file_name, 'Append', ~is_first_plot, 'Resolution', params.plot_resolution);
        is_first_plot = false;
    end
    if is_batch_mode; close(fig); end
    clear fig tl file_name;
    fprintf('done.\n');

    %% generate plot
    fprintf('Generating plot "%s" ... ', fig_name);
    if params.plot_fig_size
        fig = AKf(params.plot_fig_size(1), params.plot_fig_size(2));
    else
        fig = AKf();
    end
    set(fig, 'NumberTitle', 'Off', 'Name', fig_name);
    tl = tiledlayout(fig, 2, size(reference_spec, 3), 'TileIndexing', 'columnmajor', ...
        'TileSpacing', 'tight', 'Padding', 'tight');
    colormap(AKcolormaps('RdBu', 128, true));
    drawnow;
    
    for e = 1 : size(reference_spec, 3) % each ear
        title_str = {reference_file};
        if DO_PLOT_PRESEN; title_str = [title_str(:)', {''}]; end
        if params.plot_frac_sm
            title_str{end} = sprintf('%s   (1/%d oct. smoothing)', ...
                title_str{end}, params.plot_frac_sm);
        end
        ax = nexttile(tl);
        plot_spec_azims(reference_azel_deg(:, plot_idx), freqs, reference_spec(:, plot_idx, :), ...
            title_str, e, params.plot_freq_rng, true, params.plot_norm_freq_rng);
        if ~isempty(params.reference_file)
            if is_vertical
                xline(ax, freq_SMA_alias, ':', '$f_\textrm{A}$', lbl_config{:});
                yline(ax, 0, ':', 'LineWidth', 1, 'HandleVisibility','off');
            else
                if e == 1
                    xline(ax, 0, ':', 'ipsilateral direction $\rightarrow$', lbl_config{:});
                    yline(ax, freq_SMA_alias, ':', '$f_\textrm{A}$', lbl_config{:}, ...
                        'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'top');
                else
                    xline(ax, 0, ':', '$\leftarrow$ ipsilateral direction', lbl_config{:}, ...
                        'LabelHorizontalAlignment', 'left');
                    yline(ax, freq_SMA_alias, ':', '$f_\textrm{A}$', lbl_config{:}, ...
                        'LabelVerticalAlignment', 'top');
                end
            end
            xticklabels(''); xlabel('');
            if e == 1; set(get(ax, 'ColorBar'), 'Location', 'WestOutside'); end
            if e == size(reference_spec, 3); yticklabels(''); ylabel(''); end
        end
        drawnow;
        
        title_str = {[array_file, '   ']};
        if DO_PLOT_PRESEN; title_str = [title_str(:)', {''}]; end
        title_str{end} = sprintf('%s->   %s', title_str{end}, hrir_file);
        if params.plot_frac_sm
            title_str{end} = sprintf('%s   (1/%d oct. smoothing)', ...
                title_str{end}, params.plot_frac_sm);
        end
        ax = nexttile(tl);
        plot_spec_azims(reference_azel_deg(:, plot_idx), freqs, output_spec(:, plot_idx, :), ...
            title_str, e, params.plot_freq_rng, true, params.plot_norm_freq_rng);
        if is_vertical
            xline(ax, freq_SMA_alias, ':', '$f_\textrm{A}$', lbl_config{:});
            yline(ax, 0, ':', 'LineWidth', 1, 'HandleVisibility','off');
        else
            if e == 1
                xline(ax, 0, ':', 'ipsilateral direction $\rightarrow$', lbl_config{:});
                yline(ax, freq_SMA_alias, ':', '$f_\textrm{A}$', lbl_config{:}, ...
                    'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'top');
            else
                xline(ax, 0, ':', '$\leftarrow$ ipsilateral direction', lbl_config{:}, ...
                    'LabelHorizontalAlignment', 'left');
                yline(ax, freq_SMA_alias, ':', '$f_\textrm{A}$', lbl_config{:}, ...
                    'LabelVerticalAlignment', 'top');
            end
        end
        if e == 1; set(get(ax, 'ColorBar'), 'Location', 'WestOutside'); end
        if e == size(reference_spec, 3); yticklabels(''); ylabel(''); end
        drawnow;
    end; clear e ax;
    if is_batch_mode; clear reference_spec output_spec; end
    
    if DO_EXPORT_PLOT
        fprintf('exporting ... ');
        [~, ~] = mkdir(plot_export_dir); % ignore warning if directory already exists
        file_name = fullfile(plot_export_dir, [fig_name, '.pdf']);
        exportgraphics(fig, file_name, 'Append', true);
    end
    if is_batch_mode; close(fig); end
    clear fig tl title_str file_name;
    fprintf('done.\n');

    if ~isempty(params.reference_file)
        %% generate plot
        fprintf('Generating plot "%s" ... ', fig_name);
        if params.plot_fig_size
            fig = AKf(params.plot_fig_size(1), params.plot_fig_size(2));
        else
            fig = AKf();
        end
        set(fig, 'NumberTitle', 'Off', 'Name', fig_name);
        tl = tiledlayout(fig, 2, size(output_diff_spec, 3), 'TileIndexing', 'columnmajor', ...
            'TileSpacing', 'tight', 'Padding', 'tight');
        colormap(AKcolormaps('RdBu', 128, true));
        drawnow;
    
        title_str = 'Spectal difference';
        if params.plot_frac_sm
            title_str = sprintf('%s   (1/%d oct. smoothing)', title_str, params.plot_frac_sm);
        end
        for e = 1 : size(output_diff_spec, 3) % each ear
            ax = nexttile(tl);
            plot_spec_azims(reference_azel_deg(:, plot_idx), freqs, output_diff_spec(:, plot_idx, :), ...
                title_str, e, params.plot_freq_rng, true, params.plot_norm_freq_rng, params.plot_spec_c_rng);
            if is_vertical
                xline(ax, freq_SMA_alias, ':', '$f_\textrm{A}$', lbl_config{:});
                yline(ax, 0, ':', 'LineWidth', 1, 'HandleVisibility','off');
            else
                if e == 1
                    xline(ax, 0, ':', 'ipsilateral direction $\rightarrow$', lbl_config{:});
                    yline(ax, freq_SMA_alias, ':', '$f_\textrm{A}$', lbl_config{:}, ...
                        'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'top');
                else
                    xline(ax, 0, ':', '$\leftarrow$ ipsilateral direction', lbl_config{:}, ...
                        'LabelHorizontalAlignment', 'left');
                    yline(ax, freq_SMA_alias, ':', '$f_\textrm{A}$', lbl_config{:}, ...
                        'LabelVerticalAlignment', 'top');
                end
            end
            if e == 1; set(get(ax, 'ColorBar'), 'Location', 'WestOutside'); end
            if e == size(output_diff_spec, 3); yticklabels(''); ylabel(''); end
            drawnow;
    
            ax = nexttile(tl);
            plot_spec_azims_percentiles(freqs, output_diff_spec_perc, ...
                params.plot_percentiles, title_str, e, params.plot_freq_rng, ...
                false); % without normalization as it has been done before
            if e == size(output_diff_spec, 3); set(ax, 'YAxisLocation', 'Right'); end
            if ~DO_PLOT_PRESEN
                semilogx(freqs_log, db(output_diff_spec_log_err(:, e)), ...
                    'LineWidth', 2, 'Color', 'green', 'DisplayName', ...
                    sprintf('Average absolute  (log. mean = %.2f dB)', ...
                    mean(db(output_diff_spec_log_err(:, e)))));
            end
            xline(ax, freq_SMA_alias, ':', '$f_\textrm{A}$', lbl_config{:});
            drawnow;
        end; clear e ax;
        if is_batch_mode; clear output_diff_spec output_diff_spec_perc output_diff_spec_log_err; end
    
        if DO_EXPORT_PLOT
            fprintf('exporting ... ');
            [~, ~] = mkdir(plot_export_dir); % ignore warning if directory already exists
            file_name = fullfile(plot_export_dir, [fig_name, '.pdf']);
            exportgraphics(fig, file_name, 'Append', true, 'Resolution', params.plot_resolution);
        end
        if is_batch_mode; close(fig); end
        clear fig tl title_str file_name;
        fprintf('done.\n');

        %% generate plot
        fprintf('Generating plot "%s" ... ', fig_name);
        if params.plot_fig_size
            fig = AKf(params.plot_fig_size(1), params.plot_fig_size(2));
        else
            fig = AKf();
        end
        set(fig, 'NumberTitle', 'Off', 'Name', fig_name);
        tl = tiledlayout(fig, 2, size(output_diff_spec_w, 3), 'TileIndexing', 'columnmajor', ...
            'TileSpacing', 'tight', 'Padding', 'tight');
        colormap(AKcolormaps('RdBu', 128, true));
        drawnow;
    
        title_str = 'Weighted spectal difference';
        if params.plot_frac_sm
            title_str = sprintf('%s   (1/%d oct. smoothing)', title_str, params.plot_frac_sm);
        end
        for e = 1 : size(output_diff_spec_w, 3) % each ear
            ax = nexttile(tl);
            plot_spec_azims(reference_azel_deg(:, plot_idx), freqs, output_diff_spec_w(:, plot_idx, :), ...
                title_str, e, params.plot_freq_rng, true, params.plot_norm_freq_rng, params.plot_spec_c_rng);
            if is_vertical
                xline(ax, freq_SMA_alias, ':', '$f_\textrm{A}$', lbl_config{:});
                yline(ax, 0, ':', 'LineWidth', 1, 'HandleVisibility','off');
            else
                if e == 1
                    xline(ax, 0, ':', 'ipsilateral direction $\rightarrow$', lbl_config{:});
                    yline(ax, freq_SMA_alias, ':', '$f_\textrm{A}$', lbl_config{:}, ...
                        'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'top');
                else
                    xline(ax, 0, ':', '$\leftarrow$ ipsilateral direction', lbl_config{:}, ...
                        'LabelHorizontalAlignment', 'left');
                    yline(ax, freq_SMA_alias, ':', '$f_\textrm{A}$', lbl_config{:}, ...
                        'LabelVerticalAlignment', 'top');
                end
            end
            if e == 1; set(get(ax, 'ColorBar'), 'Location', 'WestOutside'); end
            if e == size(output_diff_spec_w, 3); yticklabels(''); ylabel(''); end
            drawnow;
    
            ax = nexttile(tl);
            plot_spec_azims_percentiles(freqs, output_diff_spec_perc_w, ...
                params.plot_percentiles, title_str, e, params.plot_freq_rng, ...
                false); % without normalization as it has been done before
            if e == size(output_diff_spec_w, 3); set(ax, 'YAxisLocation', 'Right'); end
            if ~DO_PLOT_PRESEN
                semilogx(freqs_log, db(output_diff_spec_log_err_w(:, e)), ...
                    'LineWidth', 2, 'Color', 'green', 'DisplayName', ...
                    sprintf('Average absolute  (log. mean = %.2f dB)', ...
                    mean(db(output_diff_spec_log_err_w(:, e)))));
            end
            xline(ax, freq_SMA_alias, ':', '$f_\textrm{A}$', lbl_config{:});
            drawnow;
        end; clear e ax;
        if is_batch_mode; clear output_diff_spec_w output_diff_spec_perc_w output_diff_spec_log_err_w; end
    
        if DO_EXPORT_PLOT
            fprintf('exporting ... ');
            [~, ~] = mkdir(plot_export_dir); % ignore warning if directory already exists
            file_name = fullfile(plot_export_dir, [fig_name, '.pdf']);
            exportgraphics(fig, file_name, 'Append', true, 'Resolution', params.plot_resolution);
        end
        if is_batch_mode; close(fig); end
        clear fig tl title_str file_name;
        fprintf('done.\n');
    end

    clear fig_name times freqs freqs_log plot_idx lbl_config;
else
    disp('Generating plot ... skipped.');
end

%%
[~, this_file, ~] = fileparts(mfilename('fullpath'));
fprintf('\n"%s" ... finished in %.0fh %.0fm %.0fs.\n', ...
    this_file, toc/3600, mod(toc,3600)/60, mod(toc,60));
clear this_file;


%% helper functions
function Y = getSH_AmbiEnc(order, gridAziZenRad, shDefinition)
    % from Ambisonic Encoding toolbox
    % $ git clone https://github.com/AppliedAcousticsChalmers/ambisonic-encoding.git
    Y = zeros(size(gridAziZenRad, 1), (order+1)^2);
    for n = 0 : order
        for m = -n : n
            Y(:, n^2+n+m+1) = sphharm(n, m, ...
                gridAziZenRad(:, 2), gridAziZenRad(:, 1), shDefinition);
        end
    end
end

function assertAllClose(x1, x2, norm_tol, spec_tol_dB)
    if nargin < 4; spec_tol_dB = 1; end
    if nargin < 3; norm_tol = 1e-13; end

    norm_diff = max(abs(x1 - x2), [], 'all') / max(abs([x1, x2]), [], 'all');
    if norm_diff == 0
        fprintf('no sample difference ... ');
    elseif norm_diff < norm_tol
        fprintf('%.3g max. norm. abs. sample difference ... ', norm_diff);
    else
        % compute spectral difference
        spec_diff_dB = mag2db(abs(fft(x1) ./ fft(x2)));
        % ignore 0 Hz bin
        spec_diff_dB = max(spec_diff_dB(2:end, :), [], 'all');

        if spec_diff_dB < 1e-3
            fprintf('sample but no spec. difference ... ');
        elseif spec_diff_dB < spec_tol_dB
            fprintf('sample and %.3f dB max. spec. difference ... ', spec_diff_dB);
        else
            error(['Maximum spectral mismatch of %.3f dB (greater than %.3f dB tolarance) ', ...
                'and normalized absolute sample mismatch of %.3g (greater than %.3g tolarance).'], ...
                spec_diff_dB, spec_tol_dB, norm_diff, norm_tol);
        end
    end
end
