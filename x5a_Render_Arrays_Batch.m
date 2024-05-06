% Consecutively perform an arbitrary number and combination of binaural 
% renderings from microphone array impulse responses and varying parameter
% sets.
% 
% -------------------------------------------------------------------------
%
% Hannes Helmholz 16.03.2022
%
% -------------------------------------------------------------------------
clear all; clc; close all; %#ok<CLALL> 

%% configuration
global params eq_type tStart;
eq_type = {'', 'MagLS', 'MagLSwDC', 'MagLS+SHF', 'MagLS+ADF', 'eMagLS', 'eMagLSwDC', 'eMagLSinCH', 'eMagLSinCHwDC'};

%% run
tStart = tic; % start measuring execution time

global DO_PLOT DO_PLOT_PRESEN DO_EXPORT_PLOT %#ok<*GVMIS> 
global DO_EXPORT_META DO_EXPORT_FLAC DO_EXPORT_SSR 
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

for t = 1 : length(eq_type)
    global eq_type; % this has to be restated here for some reason
    params.eq_type = eq_type{t};

    % Use HRIR set and settings that are compatible with older rendering examples
    [params.hrir_max_N, params.hrir_file] = deal(44, 'resources/HRIR_KU100/HRIR_L2702.sofa'); % KU100 for compatibility to reference BRIRs
    params.rf_len = 8192; % in samples
    params.anec_simulation = true;
    params.azim_align = false;
    params.circ_align = true;

    if any(strcmpi(params.eq_type , {'', 'eMagLS'}))
        % Simulations - vertical rotation
        params.reference_file = 'resources/HRIR_KU100/HRIR_L2702_vertical_SSR360.sofa'; % generate reference BRIRs from high-resolution SMA
        [params.N, params.array_file] = deal(29, 'resources/ARIR_WDR/SIM_VSA_SMA_LE1202_PW_struct_SubSampled_SH44.mat'); x5_Render_Arrays;
        % params.reference_file = 'resources/HRIR_KU100/HRIR_L2702_vertical_SSR360.sofa';
        params.reference_file = 'resources/BRIR_rendered/HRIR_L2702/SIM_VSA_SMA_LE1202_PW_struct_SubSampled_SH44_vertical_SSR.wav';
        [params.N, params.array_file] = deal(2, 'resources/ARIR_WDR/SIM_VSA_SMA_LE14_PW_struct_SubSampled_SH44.mat'); x5_Render_Arrays;
        [params.N, params.array_file] = deal(2, 'resources/ARIR_WDR/SIM_VSA_SMA_TD14_PW_struct_SubSampled_SH44.mat'); x5_Render_Arrays;
        [params.N, params.array_file] = deal(4, 'resources/ARIR_processed/Simulation_SMA_EM32_SrcEar.sofa'); x5_Render_Arrays; % 4.2 cm radius
        [params.N, params.array_file] = deal(4, 'resources/ARIR_WDR/SIM_VSA_SMA_LE38_PW_struct_SubSampled_SH44.mat'); x5_Render_Arrays;
        [params.N, params.array_file] = deal(4, 'resources/ARIR_WDR/SIM_VSA_SMA_TD42_PW_struct_SubSampled_SH44.mat'); x5_Render_Arrays;
        [params.N, params.array_file] = deal(8, 'resources/ARIR_WDR/SIM_VSA_SMA_LE110_PW_struct_SubSampled_SH44.mat'); x5_Render_Arrays;
        [params.N, params.array_file] = deal(8, 'resources/ARIR_WDR/SIM_VSA_SMA_TD146_PW_struct_SubSampled_SH44.mat'); x5_Render_Arrays;
        [params.N, params.array_file] = deal(12, 'resources/ARIR_WDR/SIM_VSA_SMA_LE230_PW_struct_SubSampled_SH44.mat'); x5_Render_Arrays;
        [params.N, params.array_file] = deal(12, 'resources/ARIR_WDR/SIM_VSA_SMA_TD314_PW_struct_SubSampled_SH44.mat'); x5_Render_Arrays;
        [params.N, params.array_file] = deal(2, 'resources/ARIR_WDR/SIM_VSA_EMA5_PW_struct_SubSampled_SH44.mat'); x5_Render_Arrays;
        [params.N, params.array_file] = deal(4, 'resources/ARIR_WDR/SIM_VSA_EMA9_PW_struct_SubSampled_SH44.mat'); x5_Render_Arrays;
        [params.N, params.array_file] = deal(8, 'resources/ARIR_WDR/SIM_VSA_EMA17_PW_struct_SubSampled_SH44.mat'); x5_Render_Arrays;
        [params.N, params.array_file] = deal(12, 'resources/ARIR_WDR/SIM_VSA_EMA25_PW_struct_SubSampled_SH44.mat'); x5_Render_Arrays;
        [params.N, params.array_file] = deal(29, 'resources/ARIR_WDR/SIM_VSA_EMA59_PW_struct_SubSampled_SH44.mat'); x5_Render_Arrays;

        % Control Room 1 (dry, source not centered) - vertical rotation
        params.reference_file = []; % generate reference BRIRs from high-resolution SMA
        [params.N, params.array_file] = deal(29, 'resources/ARIR_WDR/CR1_VSA_SMA_LE1202_L_struct.mat'); x5_Render_Arrays;
        params.reference_file = 'resources/BRIR_rendered/HRIR_L2702/CR1_VSA_SMA_LE1202_L_struct_vertical_SSR.wav';
        [params.N, params.array_file] = deal(2, 'resources/ARIR_WDR/CR1_VSA_SMA_LE14_L_struct_SubSampled_SH29.mat'); x5_Render_Arrays;
        [params.N, params.array_file] = deal(2, 'resources/ARIR_WDR/CR1_VSA_SMA_TD14_L_struct_SubSampled_SH29.mat'); x5_Render_Arrays;
        [params.N, params.array_file] = deal(4, 'resources/ARIR_WDR/CR1_VSA_SMA_LE38_L_struct_SubSampled_SH29.mat'); x5_Render_Arrays;
        [params.N, params.array_file] = deal(4, 'resources/ARIR_WDR/CR1_VSA_SMA_TD42_L_struct_SubSampled_SH29.mat'); x5_Render_Arrays;
        [params.N, params.array_file] = deal(8, 'resources/ARIR_WDR/CR1_VSA_SMA_LE110_L_struct_SubSampled_SH29.mat'); x5_Render_Arrays;
        [params.N, params.array_file] = deal(8, 'resources/ARIR_WDR/CR1_VSA_SMA_TD146_L_struct_SubSampled_SH29.mat'); x5_Render_Arrays;
        [params.N, params.array_file] = deal(2, 'resources/ARIR_WDR/CR1_VSA_EMA5_L_struct_SubSampled_SH29.mat'); x5_Render_Arrays;
        [params.N, params.array_file] = deal(4, 'resources/ARIR_WDR/CR1_VSA_EMA9_L_struct_SubSampled_SH29.mat'); x5_Render_Arrays;
        [params.N, params.array_file] = deal(8, 'resources/ARIR_WDR/CR1_VSA_EMA17_L_struct_SubSampled_SH29.mat'); x5_Render_Arrays;
        [params.N, params.array_file] = deal(29, 'resources/ARIR_WDR/CR1_VSA_EMA59_L_struct_SubSampled_SH29.mat'); x5_Render_Arrays;
    
        %% Large Broadcasting Studio (reverberant) - vertical rotation
        params.reference_file = []; % generate reference BRIRs from high-resolution SMA
        [params.N, params.array_file] = deal(29, 'resources/ARIR_WDR/LBS_VSA_1202RS_PAC_struct.mat'); x5_Render_Arrays;
        params.reference_file = 'resources/BRIR_rendered/HRIR_L2702/LBS_VSA_1202RS_PAC_struct_vertical_SSR.wav';
        [params.N, params.array_file] = deal(2, 'resources/ARIR_WDR/LBS_VSA_SMA_LE14_PAC_struct_SubSampled_SH29.mat'); x5_Render_Arrays;
        [params.N, params.array_file] = deal(2, 'resources/ARIR_WDR/LBS_VSA_SMA_TD14_PAC_struct_SubSampled_SH29.mat'); x5_Render_Arrays;
        [params.N, params.array_file] = deal(4, 'resources/ARIR_WDR/LBS_VSA_SMA_LE38_PAC_struct_SubSampled_SH29.mat'); x5_Render_Arrays;
        [params.N, params.array_file] = deal(4, 'resources/ARIR_WDR/LBS_VSA_SMA_TD42_PAC_struct_SubSampled_SH29.mat'); x5_Render_Arrays;
        [params.N, params.array_file] = deal(8, 'resources/ARIR_WDR/LBS_VSA_SMA_LE110_PAC_struct_SubSampled_SH29.mat'); x5_Render_Arrays;
        [params.N, params.array_file] = deal(8, 'resources/ARIR_WDR/LBS_VSA_SMA_TD146_PAC_struct_SubSampled_SH29.mat'); x5_Render_Arrays;
        [params.N, params.array_file] = deal(2, 'resources/ARIR_WDR/LBS_VSA_EMA5_PAC_struct_SubSampled_SH29.mat'); x5_Render_Arrays;
        [params.N, params.array_file] = deal(4, 'resources/ARIR_WDR/LBS_VSA_EMA9_PAC_struct_SubSampled_SH29.mat'); x5_Render_Arrays;
        [params.N, params.array_file] = deal(8, 'resources/ARIR_WDR/LBS_VSA_EMA17_PAC_struct_SubSampled_SH29.mat'); x5_Render_Arrays;
        [params.N, params.array_file] = deal(29, 'resources/ARIR_WDR/LBS_VSA_EMA59_PAC_struct_SubSampled_SH29.mat'); x5_Render_Arrays;
    end

    %% Simulations - horizontal rotation
    % The rendering arrays have a slightly different radius (8.75 cm
    % instead of 8.5 cm) than the array in `ARIR_processed` with 
    % qualitatively similar results.
    params.reference_file = 'resources/HRIR_KU100/HRIR_L2702_SSR360.sofa'; % generate reference BRIRs from high-resolution SMA
    [params.N, params.array_file] = deal(29, 'resources/ARIR_WDR/SIM_VSA_SMA_LE1202_PW_struct_SubSampled_SH44.mat'); x5_Render_Arrays;
    params.reference_file = 'resources/BRIR_rendered/HRIR_L2702/SIM_VSA_SMA_LE1202_PW_struct_SubSampled_SH44_SSR.wav';
    [params.N, params.array_file] = deal(2, 'resources/ARIR_WDR/SIM_VSA_SMA_LE14_PW_struct_SubSampled_SH44.mat'); x5_Render_Arrays;
    [params.N, params.array_file] = deal(2, 'resources/ARIR_WDR/SIM_VSA_SMA_TD14_PW_struct_SubSampled_SH44.mat'); x5_Render_Arrays;
    [params.N, params.array_file] = deal(4, 'resources/ARIR_processed/Simulation_SMA_EM32_SrcEar.sofa'); x5_Render_Arrays; % 4.2 cm radius
    [params.N, params.array_file] = deal(4, 'resources/ARIR_WDR/SIM_VSA_SMA_LE38_PW_struct_SubSampled_SH44.mat'); x5_Render_Arrays;
    [params.N, params.array_file] = deal(4, 'resources/ARIR_WDR/SIM_VSA_SMA_TD42_PW_struct_SubSampled_SH44.mat'); x5_Render_Arrays;
    [params.N, params.array_file] = deal(8, 'resources/ARIR_WDR/SIM_VSA_SMA_LE110_PW_struct_SubSampled_SH44.mat'); x5_Render_Arrays;
    [params.N, params.array_file] = deal(8, 'resources/ARIR_WDR/SIM_VSA_SMA_TD146_PW_struct_SubSampled_SH44.mat'); x5_Render_Arrays;
    [params.N, params.array_file] = deal(12, 'resources/ARIR_WDR/SIM_VSA_SMA_LE230_PW_struct_SubSampled_SH44.mat'); x5_Render_Arrays;
    [params.N, params.array_file] = deal(12, 'resources/ARIR_WDR/SIM_VSA_SMA_TD314_PW_struct_SubSampled_SH44.mat'); x5_Render_Arrays;
    [params.N, params.array_file] = deal(2, 'resources/ARIR_WDR/SIM_VSA_EMA5_PW_struct_SubSampled_SH44.mat'); x5_Render_Arrays;
    [params.N, params.array_file] = deal(4, 'resources/ARIR_WDR/SIM_VSA_EMA9_PW_struct_SubSampled_SH44.mat'); x5_Render_Arrays;
    [params.N, params.array_file] = deal(8, 'resources/ARIR_WDR/SIM_VSA_EMA17_PW_struct_SubSampled_SH44.mat'); x5_Render_Arrays;
    [params.N, params.array_file] = deal(12, 'resources/ARIR_WDR/SIM_VSA_EMA25_PW_struct_SubSampled_SH44.mat'); x5_Render_Arrays;
    [params.N, params.array_file] = deal(29, 'resources/ARIR_WDR/SIM_VSA_EMA59_PW_struct_SubSampled_SH44.mat'); x5_Render_Arrays;

    %% Anechoic measurements - horizontal rotation
    % % The rendering results are not very accurate which is probably due
    % % to non-optimal postporcessing of the anchoic SMA measurement data.
    % params.reference_file = 'resources/HRIR_KU100/HRIR_L2702_SSR360.sofa'; % generate reference BRIRs from high-resolution SMA
    % [params.N, params.array_file] = deal(29, 'resources/ARIR_WDR/ANE_VSA_SMA_LE1202_struct_SubSampled_SH44.mat'); x5_Render_Arrays;
    % params.reference_file = 'resources/BRIR_rendered/HRIR_L2702/ANE_VSA_SMA_LE1202_struct_SubSampled_SH44_SSR.wav';
    % [params.N, params.array_file] = deal(2, 'resources/ARIR_WDR/ANE_VSA_SMA_LE14_struct_SubSampled_SH44.mat'); x5_Render_Arrays;
    % [params.N, params.array_file] = deal(2, 'resources/ARIR_WDR/ANE_VSA_SMA_TD14_struct_SubSampled_SH44.mat'); x5_Render_Arrays;
    % [params.N, params.array_file] = deal(4, 'resources/ARIR_WDR/ANE_VSA_SMA_LE38_struct_SubSampled_SH44.mat'); x5_Render_Arrays;
    % [params.N, params.array_file] = deal(4, 'resources/ARIR_WDR/ANE_VSA_SMA_TD42_struct_SubSampled_SH44.mat'); x5_Render_Arrays;
    % [params.N, params.array_file] = deal(8, 'resources/ARIR_WDR/ANE_VSA_SMA_LE110_struct_SubSampled_SH44.mat'); x5_Render_Arrays;
    % [params.N, params.array_file] = deal(8, 'resources/ARIR_WDR/ANE_VSA_SMA_TD146_struct_SubSampled_SH44.mat'); x5_Render_Arrays;
    % [params.N, params.array_file] = deal(12, 'resources/ARIR_WDR/ANE_VSA_SMA_LE230_struct_SubSampled_SH44.mat'); x5_Render_Arrays;
    % [params.N, params.array_file] = deal(12, 'resources/ARIR_WDR/ANE_VSA_SMA_TD314_struct_SubSampled_SH44.mat'); x5_Render_Arrays;
    % [params.N, params.array_file] = deal(2, 'resources/ARIR_WDR/ANE_VSA_EMA5_struct_SubSampled_SH44.mat'); x5_Render_Arrays;
    % [params.N, params.array_file] = deal(4, 'resources/ARIR_WDR/ANE_VSA_EMA9_struct_SubSampled_SH44.mat'); x5_Render_Arrays;
    % [params.N, params.array_file] = deal(8, 'resources/ARIR_WDR/ANE_VSA_EMA17_struct_SubSampled_SH44.mat'); x5_Render_Arrays;
    % [params.N, params.array_file] = deal(12, 'resources/ARIR_WDR/ANE_VSA_EMA25_struct_SubSampled_SH44.mat'); x5_Render_Arrays;
    % [params.N, params.array_file] = deal(29, 'resources/ARIR_WDR/ANE_VSA_EMA59_struct_SubSampled_SH44.mat'); x5_Render_Arrays;

    %% Control Room 1 (dry, source not centered) - horizontal rotation
    % The rendering results provide a very good match with the reference
    % BRIRs, execpt a broadband attenuation at high frequencies. Also 
    % direction-dependent much smaller deviations are seen. This is may be
    % due to a better SMA microphone and sphere hardware and / or due not 
    % using the KEMAR dummy head with a rather loose fitting of the ears.
    params.reference_file = 'resources/ARIR_WDR/BRIR_CR1_KU_ROTM_L.sofa'; % generate reference BRIRs from high-resolution SMA
    [params.N, params.array_file] = deal(29, 'resources/ARIR_WDR/CR1_VSA_SMA_LE1202_L_struct.mat'); x5_Render_Arrays;
    params.reference_file = 'resources/BRIR_rendered/HRIR_L2702/CR1_VSA_SMA_LE1202_L_struct_SSR.wav';
    [params.N, params.array_file] = deal(2, 'resources/ARIR_WDR/CR1_VSA_SMA_LE14_L_struct_SubSampled_SH29.mat'); x5_Render_Arrays;
    [params.N, params.array_file] = deal(2, 'resources/ARIR_WDR/CR1_VSA_SMA_TD14_L_struct_SubSampled_SH29.mat'); x5_Render_Arrays;
    [params.N, params.array_file] = deal(4, 'resources/ARIR_WDR/CR1_VSA_SMA_LE38_L_struct_SubSampled_SH29.mat'); x5_Render_Arrays;
    [params.N, params.array_file] = deal(4, 'resources/ARIR_WDR/CR1_VSA_SMA_TD42_L_struct_SubSampled_SH29.mat'); x5_Render_Arrays;
    [params.N, params.array_file] = deal(8, 'resources/ARIR_WDR/CR1_VSA_SMA_LE110_L_struct_SubSampled_SH29.mat'); x5_Render_Arrays;
    [params.N, params.array_file] = deal(8, 'resources/ARIR_WDR/CR1_VSA_SMA_TD146_L_struct_SubSampled_SH29.mat'); x5_Render_Arrays;
    [params.N, params.array_file] = deal(12, 'resources/ARIR_WDR/CR1_VSA_SMA_LE230_L_struct_SubSampled_SH29.mat'); x5_Render_Arrays;
    [params.N, params.array_file] = deal(12, 'resources/ARIR_WDR/CR1_VSA_SMA_TD314_L_struct_SubSampled_SH29.mat'); x5_Render_Arrays;
    [params.N, params.array_file] = deal(2, 'resources/ARIR_WDR/CR1_VSA_EMA5_L_struct_SubSampled_SH29.mat'); x5_Render_Arrays;
    [params.N, params.array_file] = deal(4, 'resources/ARIR_WDR/CR1_VSA_EMA9_L_struct_SubSampled_SH29.mat'); x5_Render_Arrays;
    [params.N, params.array_file] = deal(8, 'resources/ARIR_WDR/CR1_VSA_EMA17_L_struct_SubSampled_SH29.mat'); x5_Render_Arrays;
    [params.N, params.array_file] = deal(12, 'resources/ARIR_WDR/CR1_VSA_EMA25_L_struct_SubSampled_SH29.mat'); x5_Render_Arrays;
    [params.N, params.array_file] = deal(29, 'resources/ARIR_WDR/CR1_VSA_EMA59_L_struct_SubSampled_SH29.mat'); x5_Render_Arrays;

    %% Large Broadcasting Studio (reverberant) - horizontal rotation
    % Rendering results provide a very good match with the reference
    % similar to CR1.
    params.reference_file = 'resources/ARIR_WDR/BRIR_LBS_KU_ROTM_PAC.sofa'; % generate reference BRIRs from high-resolution SMA
    [params.N, params.array_file] = deal(29, 'resources/ARIR_WDR/LBS_VSA_1202RS_PAC_struct.mat'); x5_Render_Arrays;
    params.reference_file = 'resources/BRIR_rendered/HRIR_L2702/LBS_VSA_1202RS_PAC_struct_SSR.wav';
    [params.N, params.array_file] = deal(2, 'resources/ARIR_WDR/LBS_VSA_SMA_LE14_PAC_struct_SubSampled_SH29.mat'); x5_Render_Arrays;
    [params.N, params.array_file] = deal(2, 'resources/ARIR_WDR/LBS_VSA_SMA_TD14_PAC_struct_SubSampled_SH29.mat'); x5_Render_Arrays;
    [params.N, params.array_file] = deal(4, 'resources/ARIR_WDR/LBS_VSA_SMA_LE38_PAC_struct_SubSampled_SH29.mat'); x5_Render_Arrays;
    [params.N, params.array_file] = deal(4, 'resources/ARIR_WDR/LBS_VSA_SMA_TD42_PAC_struct_SubSampled_SH29.mat'); x5_Render_Arrays;
    [params.N, params.array_file] = deal(8, 'resources/ARIR_WDR/LBS_VSA_SMA_LE110_PAC_struct_SubSampled_SH29.mat'); x5_Render_Arrays;
    [params.N, params.array_file] = deal(8, 'resources/ARIR_WDR/LBS_VSA_SMA_TD146_PAC_struct_SubSampled_SH29.mat'); x5_Render_Arrays;
    [params.N, params.array_file] = deal(12, 'resources/ARIR_WDR/LBS_VSA_SMA_LE230_PAC_struct_SubSampled_SH29.mat'); x5_Render_Arrays;
    [params.N, params.array_file] = deal(12, 'resources/ARIR_WDR/LBS_VSA_SMA_TD314_PAC_struct_SubSampled_SH29.mat'); x5_Render_Arrays;
    [params.N, params.array_file] = deal(2, 'resources/ARIR_WDR/LBS_VSA_EMA5_PAC_struct_SubSampled_SH29.mat'); x5_Render_Arrays;
    [params.N, params.array_file] = deal(4, 'resources/ARIR_WDR/LBS_VSA_EMA9_PAC_struct_SubSampled_SH29.mat'); x5_Render_Arrays;
    [params.N, params.array_file] = deal(8, 'resources/ARIR_WDR/LBS_VSA_EMA17_PAC_struct_SubSampled_SH29.mat'); x5_Render_Arrays;
    [params.N, params.array_file] = deal(12, 'resources/ARIR_WDR/LBS_VSA_EMA25_PAC_struct_SubSampled_SH29.mat'); x5_Render_Arrays;
    [params.N, params.array_file] = deal(29, 'resources/ARIR_WDR/LBS_VSA_EMA59_PAC_struct_SubSampled_SH29.mat'); x5_Render_Arrays;
end

%%
clear;
global params tStart; % this has to be restated here for some reason
[~, this_file, ~] = fileparts(mfilename('fullpath'));
fprintf('\n"%s" ... finished in %.0fh %.0fm %.0fs.\n', ...
    this_file, toc(tStart)/3600, mod(toc(tStart),3600)/60, mod(toc(tStart),60));
clear this_file;
