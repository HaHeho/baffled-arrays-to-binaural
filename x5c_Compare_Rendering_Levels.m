% Compare arbitrary combinations of (rendered) binaural room impulse 
% responses by plotting various resulting signal levels. This is helpful 
% for the preparation of the rendered BRIRs for the perceptual comparison 
% in a user study, where all stimuli should ideally be loudness normalized.
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
% Hannes Helmholz 14.05.2024
%
% -------------------------------------------------------------------------
clear; clc; close all;

addpath(genpath('dependencies'));

%% configuration
data_files = 'resources/BRIR_rendered/Kemar_HRTF_sofa_N44_adjusted/*SSR.wav';
% data_files = 'resources/BRIR_rendered/Kemar_HRTF_sofa_N44_adjusted/Simulation_*SSR.wav';
% data_files = 'resources/BRIR_rendered/Kemar_HRTF_sofa_N44_adjusted/Anechoic_*SSR.wav';
% data_files = 'resources/BRIR_rendered/Kemar_HRTF_sofa_N44_adjusted/LabWet_*SSR.wav';

level_plot_lim = [-5, 5]; % in dB

%%
tic; % start measuring execution time

fprintf('Parsing files ... in "%s" ... ', data_files);
SOFAstart;
files = dir(data_files);
files = files(~ismember({files.name}, {'.', '..', '.DS_Store'}));
files = natsortfiles(files);
fprintf('found %d files ... ', length(files));
current_file = fileparts(mfilename('fullpath'));
fprintf('done.\n');

assert(length(files), 'No files found.');
for f = 1 : length(files)
    file = fullfile(files(f).folder, files(f).name);
    file = file(length(current_file)+2:end); % remove absolute path

    fprintf('Loading file "%s" ... ', file);
    [~, files(f).config, ~] = fileparts(file);
    files(f).config = files(f).config(1:end-4);
    if endsWith(file, '.wav')
        [tmp, fs] = audioread(file);
        tmp = unstack_SSR_BRIRs(tmp);

        files(f).peak_front = db(max(abs(tmp(:, 1, :)), [], 'all'));
        files(f).peak_all = db(max(abs(tmp), [], 'all'));
        files(f).rssq_front = db(mean(rssq(tmp(:, 1, :)), 'all'));
        files(f).rssq_all = db(mean(rssq(tmp), 'all'));
    else
        error('Unknown data file format.');
    end
    clear tmp;
    fprintf('done.\n');
end; clear f file;

% normalize mean levels of all configurations by RSSQ
mean_lvl = mean([files.rssq_front]);
for f = 1 : length(files)
    files(f).rssq_front = files(f).rssq_front - mean_lvl;
end; clear f mean_lvl;
mean_lvl = mean([files.rssq_all]);
for f = 1 : length(files)
    files(f).rssq_all = files(f).rssq_all - mean_lvl;
end; clear f mean_lvl;

% % normalize mean levels of all configurations by peak
% mean_lvl = mean([files.peak_front]);
% for f = 1 : length(files)
%     files(f).peak_front = files(f).peak_front - mean_lvl;
% end; clear f mean_lvl;
% mean_lvl = mean([files.peak_all]);
% for f = 1 : length(files)
%     files(f).peak_all = files(f).peak_all - mean_lvl;
% end; clear f mean_lvl;

[files.offset] = deal(nan);
offset_lvl = 0;

%%
fig_name = sprintf('Compare levels peak "%s"', data_files); 
fprintf('Generating plot "%s" ... ', fig_name);
fig = AKf();
set(fig, 'NumberTitle', 'Off', 'Name', fig_name);
tl = tiledlayout(fig, 2, 1, 'TileSpacing', 'none', 'Padding', 'tight');
drawnow;

nexttile(tl);
bh = bar([files.peak_all], 'g', 'FaceColor', 'flat');
set(bh, 'ButtonDownFcn', @ButtonDownFcnBar);
set(bh, 'BaseValue', floor(min([files.peak_all])));
ylabel('Peak (all) level in dB');
grid on;
axis tight;
ylim([floor(min([files.peak_all])), ceil(max([files.peak_all]))]);
xticks(1 : length(files));
xticklabels(cellfun(@(a) strrep(a, '_', '\_'), {files.config}, 'UniformOutput', false));
set(gca, 'TickLabelInterpreter', 'Tex');
drawnow;

nexttile(tl);
bh = bar([files.peak_front], 'g', 'FaceColor', 'flat');
set(bh, 'ButtonDownFcn', @ButtonDownFcnBar);
set(bh, 'BaseValue', floor(min([files.peak_front])));
hold on;
plot(offset_lvl - [files.offset], 'LineStyle', 'none', 'Marker', '*', 'MarkerSize', 15);
ylabel('Peak (front) level in dB');
grid on;
axis tight;
ylim([floor(min([files.peak_front])), ceil(max([files.peak_front]))]);
xticks(1 : length(files));
xticklabels(cellfun(@(a) strrep(a, '_', '\_'), {files.config}, 'UniformOutput', false));
drawnow;

clear fig_name fig tl bh;
fprintf('done.\n');
fprintf('\n');

%%
fig_name = sprintf('Compare levels RSSQ "%s"', data_files); 
fprintf('Generating plot "%s" ... ', fig_name);
fig = AKf();
set(fig, 'NumberTitle', 'Off', 'Name', fig_name);
tl = tiledlayout(fig, 2, 1, 'TileSpacing', 'none', 'Padding', 'tight');
drawnow;

nexttile(tl);
bh = bar([files.rssq_all], 'g', 'FaceColor', 'flat');
set(bh, 'ButtonDownFcn', @ButtonDownFcnBar);
ylabel('Relative RSSQ (all) level in dB');
grid on;
axis tight;
ylim(level_plot_lim);
yticks(level_plot_lim(1) : level_plot_lim(2));
xticks(1 : length(files));
xticklabels(cellfun(@(a) strrep(a, '_', '\_'), {files.config}, 'UniformOutput', false));
set(gca, 'TickLabelInterpreter', 'Tex');
drawnow;

nexttile(tl);
bh = bar([files.rssq_front], 'g', 'FaceColor', 'flat');
set(bh, 'ButtonDownFcn', @ButtonDownFcnBar);
hold on;
plot(offset_lvl - [files.offset], 'LineStyle', 'none', 'Marker', '*', 'MarkerSize', 15);
ylabel('Relative RSSQ (front) level in dB');
grid on;
axis tight;
ylim(level_plot_lim);
yticks(level_plot_lim(1) : level_plot_lim(2));
xticks(1 : length(files));
xticklabels(cellfun(@(a) strrep(a, '_', '\_'), {files.config}, 'UniformOutput', false));
drawnow;

clear fig_name fig tl bh;
fprintf('done.\n');
fprintf('\n');

%%
fprintf(' ... finished in %.0fh %.0fm %.0fs.\n', ...
    toc/3600, mod(toc,3600)/60, mod(toc,60));
