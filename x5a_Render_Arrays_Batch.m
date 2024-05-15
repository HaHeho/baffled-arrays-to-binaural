% Consecutively perform an arbitrary combination of binaural renderings 
% from microphone array impulse responses and varying parameter sets.
% 
% -------------------------------------------------------------------------
%
% Hannes Helmholz 13.10.2023
%
% -------------------------------------------------------------------------
clear all; clc; close all; %#ok<CLALL> 

%% configuration
global params eq_type tStart;
eq_type = {'', 'eMagLS', 'eBreve'};

%% run
tStart = tic; % start measuring execution time

global DO_PLOT DO_PLOT_PRESEN DO_EXPORT_PLOT %#ok<*GVMIS> 
global DO_EXPORT_META DO_EXPORT_AURA DO_EXPORT_BRIR 
DO_PLOT        = true;
DO_PLOT_PRESEN = false;
DO_EXPORT_PLOT = true;
DO_EXPORT_META = true;
DO_EXPORT_BRIR = true;
DO_EXPORT_AURA = false;

% DO_PLOT        = false;
% DO_PLOT_PRESEN = true;
% DO_EXPORT_PLOT = false;
% DO_EXPORT_META = false;
% DO_EXPORT_BRIR = false;
% DO_EXPORT_AURA = true;

for t = 1 : length(eq_type)
    global eq_type; % this has to be restated here for some reason
    params.eq_type = eq_type{t};

    %% Simulations - useful for debugging
    % params.reference_file = 'resources/HRIR_KEMAR/Kemar_HRTF_sofa_N44_adjusted_SSR.wav';
    % [params.N, params.array_file] = deal(44, 'resources/ARIR_processed/Simulation_SMA_LE2702_SrcEar.sofa'); x5_Render_Arrays; % generate reference BRIRs from high-resolution SMA
    % params.reference_file = 'resources/BRIR_rendered/Kemar_HRTF_sofa_N44_adjusted/Simulation_SMA_LE2702_SrcEar_SSR.wav';
    % [params.N, params.array_file] = deal(1, 'resources/ARIR_processed/Simulation_SMA_TD6_SrcEar_SubSampled_SH44.sofa'); x5_Render_Arrays;
    % [params.N, params.array_file] = deal(2, 'resources/ARIR_processed/Simulation_SMA_TD14_SrcEar_SubSampled_SH44.sofa'); x5_Render_Arrays;
    % [params.N, params.array_file] = deal(4, 'resources/ARIR_processed/Simulation_SMA_EM32_SrcEar.sofa'); x5_Render_Arrays;
    % [params.N, params.array_file] = deal(4, 'resources/ARIR_processed/Simulation_SMA_TD42_SrcEar_SubSampled_SH44.sofa'); x5_Render_Arrays;
    % [params.N, params.array_file] = deal(8, 'resources/ARIR_processed/Simulation_SMA_TD146_SrcEar_SubSampled_SH44.sofa'); x5_Render_Arrays;
    % % [params.N, params.array_file] = deal(12, 'resources/ARIR_processed/Simulation_SMA_TD314_SrcEar_SubSampled_SH44.sofa'); x5_Render_Arrays;
    % % [params.N, params.array_file] = deal(29, 'resources/ARIR_processed/Simulation_SMA_TD1742_SrcEar_SubSampled_SH44.sofa'); x5_Render_Arrays;
    % [params.N, params.array_file] = deal(2, 'resources/ARIR_processed/Simulation_EMA5_SrcEar_SubSampled_SH44.sofa'); x5_Render_Arrays;
    % [params.N, params.array_file] = deal(4, 'resources/ARIR_processed/Simulation_EMA9_SrcEar_SubSampled_SH44.sofa'); x5_Render_Arrays;
    % [params.N, params.array_file] = deal(8, 'resources/ARIR_processed/Simulation_EMA17_SrcEar_SubSampled_SH44.sofa'); x5_Render_Arrays;
    % % [params.N, params.array_file] = deal(12, 'resources/ARIR_processed/Simulation_EMA25_SrcEar_SubSampled_SH44.sofa'); x5_Render_Arrays;
    % % [params.N, params.array_file] = deal(29, 'resources/ARIR_processed/Simulation_EMA59_SrcEar_SubSampled_SH44.sofa'); x5_Render_Arrays;
    % [params.N, params.array_file] = deal(44, 'resources/ARIR_processed/Simulation_EMA89_SrcEar_SubSampled_SH44.sofa'); x5_Render_Arrays;

    %% Anechoic measurements - useful for debugging
    % params.reference_file = 'resources/ARIR_processed/Anechoic_KEMAR_SrcEar.sofa';
    % [params.N, params.array_file] = deal(44, 'resources/ARIR_processed/Anechoic_SMA_LE2702_SrcEar.sofa'); x5_Render_Arrays; % generate reference BRIRs from high-resolution SMA
    % params.reference_file = 'resources/BRIR_rendered/Kemar_HRTF_sofa_N44_adjusted/Anechoic_SMA_LE2702_SrcEar_SSR.wav';
    % [params.N, params.array_file] = deal(1, 'resources/ARIR_processed/Anechoic_SMA_TD6_SrcEar_SubSampled_SH44.sofa'); x5_Render_Arrays;
    % [params.N, params.array_file] = deal(2, 'resources/ARIR_processed/Anechoic_SMA_TD14_SrcEar_SubSampled_SH44.sofa'); x5_Render_Arrays;
    % [params.N, params.array_file] = deal(4, 'resources/ARIR_processed/Anechoic_SMA_TD42_SrcEar_SubSampled_SH44.sofa'); x5_Render_Arrays;
    % [params.N, params.array_file] = deal(8, 'resources/ARIR_processed/Anechoic_SMA_TD146_SrcEar_SubSampled_SH44.sofa'); x5_Render_Arrays;
    % % [params.N, params.array_file] = deal(12, 'resources/ARIR_processed/Anechoic_SMA_TD314_SrcEar_SubSampled_SH44.sofa'); x5_Render_Arrays;
    % % [params.N, params.array_file] = deal(29, 'resources/ARIR_processed/Anechoic_SMA_TD1742_SrcEar_SubSampled_SH44.sofa'); x5_Render_Arrays;
    % [params.N, params.array_file] = deal(2, 'resources/ARIR_processed/Anechoic_EMA5_SrcEar_SubSampled_SH44.sofa'); x5_Render_Arrays;
    % [params.N, params.array_file] = deal(4, 'resources/ARIR_processed/Anechoic_EMA9_SrcEar_SubSampled_SH44.sofa'); x5_Render_Arrays;
    % [params.N, params.array_file] = deal(8, 'resources/ARIR_processed/Anechoic_EMA17_SrcEar_SubSampled_SH44.sofa'); x5_Render_Arrays;
    % % [params.N, params.array_file] = deal(12, 'resources/ARIR_processed/Anechoic_EMA25_SrcEar_SubSampled_SH44.sofa'); x5_Render_Arrays;
    % % [params.N, params.array_file] = deal(29, 'resources/ARIR_processed/Anechoic_EMA59_SrcEar_SubSampled_SH44.sofa'); x5_Render_Arrays;
    % [params.N, params.array_file] = deal(44, 'resources/ARIR_processed/Anechoic_EMA89_SrcEar_SubSampled_SH44.sofa'); x5_Render_Arrays;
    % [params.N, params.array_file] = deal(2, 'resources/ARIR_processed/Anechoic_XMA6_SrcEar.sofa'); x5_Render_Arrays;
    % [params.N, params.array_file] = deal(4, 'resources/ARIR_processed/Anechoic_XMA9_SrcEar.sofa'); x5_Render_Arrays;
    % [params.N, params.array_file] = deal(8, 'resources/ARIR_processed/Anechoic_XMA18_SrcEar.sofa'); x5_Render_Arrays;
    % [params.N, params.array_file] = deal(2, 'resources/ARIR_processed/Anechoic_XMA6v2_SrcEar.sofa'); x5_Render_Arrays;
    % [params.N, params.array_file] = deal(4, 'resources/ARIR_processed/Anechoic_XMA9v2_SrcEar.sofa'); x5_Render_Arrays;
    % [params.N, params.array_file] = deal(8, 'resources/ARIR_processed/Anechoic_XMA18v2_SrcEar.sofa'); x5_Render_Arrays;

    %% Lab (dry) with floor reflection
    params.reference_file = 'resources/ARIR_processed/LabWet_KEMAR_SrcEar.sofa';
    [params.N, params.array_file] = deal(44, 'resources/ARIR_processed/LabWet_SMA_LE2702_SrcEar.sofa'); x5_Render_Arrays; % generate reference BRIRs from high-resolution SMA
    params.reference_file = 'resources/BRIR_rendered/Kemar_HRTF_sofa_N44_adjusted/LabWet_SMA_LE2702_SrcEar_SSR.wav';
    [params.N, params.array_file] = deal(1, 'resources/ARIR_processed/LabWet_SMA_TD6_SrcEar_SubSampled_SH44.sofa'); x5_Render_Arrays;
    [params.N, params.array_file] = deal(2, 'resources/ARIR_processed/LabWet_SMA_TD14_SrcEar_SubSampled_SH44.sofa'); x5_Render_Arrays;
    [params.N, params.array_file] = deal(4, 'resources/ARIR_processed/LabWet_SMA_TD42_SrcEar_SubSampled_SH44.sofa'); x5_Render_Arrays;
    [params.N, params.array_file] = deal(8, 'resources/ARIR_processed/LabWet_SMA_TD146_SrcEar_SubSampled_SH44.sofa'); x5_Render_Arrays;
    % [params.N, params.array_file] = deal(12, 'resources/ARIR_processed/LabWet_SMA_TD314_SrcEar_SubSampled_SH44.sofa'); x5_Render_Arrays;
    % [params.N, params.array_file] = deal(29, 'resources/ARIR_processed/LabWet_SMA_TD1742_SrcEar_SubSampled_SH44.sofa'); x5_Render_Arrays;
    [params.N, params.array_file] = deal(2, 'resources/ARIR_processed/LabWet_EMA5_SrcEar_SubSampled_SH44.sofa'); x5_Render_Arrays;
    [params.N, params.array_file] = deal(4, 'resources/ARIR_processed/LabWet_EMA9_SrcEar_SubSampled_SH44.sofa'); x5_Render_Arrays;
    [params.N, params.array_file] = deal(8, 'resources/ARIR_processed/LabWet_EMA17_SrcEar_SubSampled_SH44.sofa'); x5_Render_Arrays;
    % [params.N, params.array_file] = deal(12, 'resources/ARIR_processed/LabWet_EMA25_SrcEar_SubSampled_SH44.sofa'); x5_Render_Arrays;
    % [params.N, params.array_file] = deal(29, 'resources/ARIR_processed/LabWet_EMA59_SrcEar_SubSampled_SH44.sofa'); x5_Render_Arrays;
    [params.N, params.array_file] = deal(44, 'resources/ARIR_processed/LabWet_EMA89_SrcEar_SubSampled_SH44.sofa'); x5_Render_Arrays;
    [params.N, params.array_file] = deal(2, 'resources/ARIR_processed/LabWet_XMA6_SrcEar.sofa'); x5_Render_Arrays;
    [params.N, params.array_file] = deal(4, 'resources/ARIR_processed/LabWet_XMA9_SrcEar.sofa'); x5_Render_Arrays;
    [params.N, params.array_file] = deal(8, 'resources/ARIR_processed/LabWet_XMA18_SrcEar.sofa'); x5_Render_Arrays;
end

%%
clear;
global params tStart; % this has to be restated here for some reason
[~, this_file, ~] = fileparts(mfilename('fullpath'));
fprintf('\n"%s" ... finished in %.0fh %.0fm %.0fs.\n', ...
    this_file, toc(tStart)/3600, mod(toc(tStart),3600)/60, mod(toc(tStart),60));
clear this_file;
