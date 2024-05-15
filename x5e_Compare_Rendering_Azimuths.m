% Generate an interactive plot to compare arbitrary combinations of 
% rendered binaural room impulse responses by their azimuth alignment 
% angle applied during binaural rendering. This is helpful for directly 
% comparing similar rendering configurations and validating the 
% rendering method or detecting potentially flawed data sets.
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
data_files = 'resources/BRIR_rendered/Kemar_HRTF_sofa_N44_adjusted/*_meta.mat';
% data_files = 'resources/BRIR_rendered/Kemar_HRTF_sofa_N44_adjusted/Simulation_*_meta.mat';
% data_files = 'resources/BRIR_rendered/Kemar_HRTF_sofa_N44_adjusted/Anechoic_*_meta.mat';
% data_files = 'resources/BRIR_rendered/Kemar_HRTF_sofa_N44_adjusted/LabWet_*_meta.mat';

shift_plot_lim = [-50, 50]; % in deg

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
    files(f).config = files(f).config(1:end-5);
    if endsWith(file, '.mat')
        tmp = load(file);
        files(f).shift = tmp.azim_min_shift;
    else
        error('Unknown data file format.');
    end
    clear tmp;
    fprintf('done.\n');
end; clear f file;

%%
fig_name = sprintf('Compare azimuth shifts "%s"', data_files); 
fprintf('Generating plot "%s" ... ', fig_name);
fig = AKf();
set(fig, 'NumberTitle', 'Off', 'Name', fig_name);
tl = tiledlayout(fig, 1, 1, 'TileSpacing', 'tight', 'Padding', 'tight');
drawnow;
nexttile(tl);
bh = bar([files.shift], 'g', 'FaceColor', 'flat');
set(bh, 'ButtonDownFcn', @Button_Down_Fcn_Bar);
axis tight;
xticks(1 : length(files));
xticklabels(cellfun(@(a) strrep(a, '_', '\_'), {files.config}, 'UniformOutput', false));
yticks(-180 : 10 : 180);
ylim(shift_plot_lim);
ylabel('Azimuth in deg');
set(gca, 'TickLabelInterpreter', 'Tex');
grid on;
text(0.005, 0.99, {'Interact to select a configuration!'; ...
    'Subsequently, the configuration name and azimuth alignment angle are copied and may be pasted from the clipboard.'}, ...
    'Color', 'b', 'FontSize', 14, 'VerticalAlignment', 'Top', 'Units', 'normalized');
drawnow;
clear fig_name fig tl bh;
fprintf('done.\n');
fprintf('\n');

%%
fprintf(' ... finished in %.0fh %.0fm %.0fs.\n', ...
    toc/3600, mod(toc,3600)/60, mod(toc,60));
