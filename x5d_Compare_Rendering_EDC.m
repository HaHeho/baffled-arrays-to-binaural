% Generate an interactive plot to compare arbitrary combinations of 
% (rendered) binaural room impulse responses by their Energy Decay Curve. 
% This is helpful for directly comparing similar rendering configurations 
% and validating the rendering method or detecting potentially flawed 
% data sets.
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

edc_plot_dr = 40; % in dB

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

        files(f).etc_front = squeeze(db(abs(tmp(:, 1, :)))); % in dB
        files(f).edc_front = squeeze(AKedc(tmp(:, 1, :))); % in dB
        files(f).edc_all = squeeze(max(AKedc(tmp), [], 2)); % in dB

        % max of both ears
        files(f).etc_front = max(files(f).etc_front, [], 2);
        files(f).edc_front = max(files(f).edc_front, [], 2);
        files(f).edc_all = max(files(f).edc_all, [], 2);

        % normalize ETC
        files(f).etc_front = files(f).etc_front - max(files(f).etc_front);
    else
        error('Unknown data file format.');
    end
    clear tmp;
    fprintf('done.\n');
end; clear f file;

% pad all EDCs to the same length
max_len = max([cellfun(@length, {files.etc_front})]);
for f = 1 : length(files)
    files(f).etc_front(end+1:max_len, :) = -Inf;
    files(f).edc_front(end+1:max_len, :) = -Inf;
    files(f).edc_all(end+1:max_len, :) = -Inf;
end; clear f;

% get onset
edc_onset = Inf;
for f = 1 : length(files)
    edc_onset = min([edc_onset, find(files(f).edc_all >= -1, 1, 'last')]);
end; clear f;
edc_onset = floor(max([1, edc_onset - (.05 * max_len)]));

% get end
edc_end = -Inf;
for f = 1 : length(files)
    edc_end = max([edc_end, find(files(f).edc_all >= -edc_plot_dr, 1, 'last')]);
end; clear f;
edc_end = ceil(min([max_len, edc_end + (.05 * max_len)]));
clear max_len;

% define plot colors
colors = AKcolormaps('Dark2', length(files));

%%
fig_name = sprintf('Compare EDC "%s"', data_files); 
fprintf('Generating plot "%s" ... ', fig_name);
fig = AKf();
set(fig, 'NumberTitle', 'Off', 'Name', fig_name);
% tl = tiledlayout(fig, 3, 1, 'TileSpacing', 'none', 'Padding', 'tight');
tl = tiledlayout(fig, 2, 1, 'TileSpacing', 'tight', 'Padding', 'tight');
drawnow;

ax = nexttile(tl);
h_etc = plot(ax, [files.etc_front], 'LineWidth', 1);
colororder(ax, colors);
grid on;
xlim([edc_onset, edc_end]);
ylim([-edc_plot_dr, 2]);
ylabel('Normalized Energy Time Curve (front MAX) in dB');
xlabel('Time in samples');
drawnow;

ax = nexttile(tl);
h_front = plot(ax, [files.edc_front], 'LineWidth', 3);
colororder(ax, colors);
grid on;
xlim([edc_onset, edc_end]);
ylim([-edc_plot_dr, 2]);
ylabel('Energy Decay Curve (front MAX) in dB');
xlabel('Time in samples');
legend({files.config}, 'Interpreter', 'none', 'FontSize', 13, 'Location', 'NorthEast');
text(0.005, 0.02, {'Interact to select a configuration!'; ...
    'Subsequently, both the ETC and the EDC curves are highlighted for identification.'}, ...
    'Color', 'b', 'FontSize', 14, 'VerticalAlignment', 'Bottom', 'Units', 'normalized');
drawnow;

% ax = nexttile(tl);
% h_all = plot(ax, [files.edc_all], 'LineWidth', 3);
% colororder(ax, colors);
% grid on;
% xlim([edc_onset, edc_end]);
% ylim([-edc_plot_dr, 2]);
% ylabel('Energy Decay Curve (all MAX) in dB');
% xlabel('Time in samples');
% legend({files.config}, 'Interpreter', 'none', 'FontSize', 13, 'Location', 'NorthEast');
% drawnow;

% add triggers to highlight lines
for h = 1 : length(h_front)
    set(h_etc(h), 'ButtonDownFcn', {@mark_line_callback, h_etc, h_front, colors});
    set(h_front(h), 'ButtonDownFcn', {@mark_line_callback, h_front, h_etc, colors});
    % set(h_all(h), 'ButtonDownFcn', {@mark_line_callback, h_all, h_etc, colors});
end; clear h h_etc h_all h_front;

clear fig_name fig tl ax;
fprintf('done.\n');
fprintf('\n');

%%
fprintf(' ... finished in %.0fh %.0fm %.0fs.\n', ...
    toc/3600, mod(toc,3600)/60, mod(toc,60));


%% helper functions
function mark_line_callback(handle, ~, list, otherLists, orig_colors)
    MARK_COLOR = [0, 1, 0];

    list_idx = find(list == handle);
    fprintf('Selected position index %d.\n', list_idx);
    for i = 1 : length(otherLists)
        otherHandle(i) = otherLists(list_idx); %#ok<AGROW>
    end
    
    % check if line is already active
    if all(handle.Color == MARK_COLOR)
        return;
    end
    
    % reset inactive line color
    for i = 1 : length(list)
        set(list(i), 'Color', orig_colors(i, :));
        set(otherLists(i), 'Color', orig_colors(i, :));
    end

    % set active line color
    handle.Color = MARK_COLOR;
    for i = 1 : length(otherHandle)
        otherHandle(i).Color = MARK_COLOR;
    end

    % set active line before all others
    uistack(handle, 'top');
    for i = 1 : length(otherHandle)
        uistack(otherHandle(i), 'top');
    end
end
