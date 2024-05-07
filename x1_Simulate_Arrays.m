% Simulate a plane wave impinging from an arbitrary direction on SMAs and 
% EMAs with a desired sampling grid in an anechoic environment.
%
% Such simulations are helpful to evaluate the rendering method and to 
% investigate the influence of different sampling grids and equalization 
% methods on the rendered binaural signals.
%
% -------------------------------------------------------------------------
%
% requires Array Response Simulator toolbox
% $ git clone https://github.com/polarch/Array-Response-Simulator.git
% 
% requires AKtools toolbox (run AKtoolsStart.m)
% $ git clone https://github.com/f-brinkmann/AKtools.git
%
%  -------------------------------------------------------------------------
%
% Hannes Helmholz 13.03.2023
%
% -------------------------------------------------------------------------
clear; clc; close all;

addpath(genpath('dependencies'));

%% configuration
[array_order, array_type, array_r] = deal(4, 'Eigenmike', 0.042);
% % Radius matching the Chalmers VariSphear scattering body
% [array_order, array_type, array_r] = deal(8, 'tDesign', 0.085); % for subsampling validation
% [array_order, array_type, array_r] = deal(44, 'Lebedev', 0.085);

sim_DOA_deg = [0, 0]; % in deg, azimuth and elevation
sim_order = 44; % SH order for anechoic simulation
sim_array_type = 'rigid';
sim_len = 2048; % in samples
sim_fs = 48e3; % in Hz
sim_C = 343; % in m/s, speed of sound, in accordance with `simulateSphArray`

plot_export_dir = 'plots/1_Simulation';
data_export_dir = 'resources/ARIR_processed';

global DO_PLOT DO_PLOT_EXPORT DO_FILE_EXPORT %#ok<GVMIS> 
DO_PLOT        = true;
DO_PLOT_EXPORT = true;
DO_FILE_EXPORT = true;

% DO_PLOT        = false;
% DO_PLOT_EXPORT = false;
% DO_FILE_EXPORT = false;

%% prepare data
SOFAstart;
tic; % start measuring execution time

fprintf('Generating target %s grid at N=%d ... ', array_type, array_order);
[array_rad, array_weights] = get_array_grid( ...
    array_type, array_order, true, true); % in rad, azimuth and colatitude
array_ch = size(array_rad, 1);
if strcmpi(array_type, 'Fliege')
    array_name = sprintf('SMA_FM%d', array_ch);
elseif strcmpi(array_type, 'Lebedev')
    array_name = sprintf('SMA_LE%d', array_ch);
elseif strcmpi(array_type, 'tDesign')
    array_name = sprintf('SMA_TD%d', array_ch);
elseif any(strcmpi(array_type, {'Eigenmike', 'Eigenmike32'}))
    array_name = sprintf('SMA_EM%d', array_ch);
elseif strcmpi(array_type, 'Equatorial')
    array_name = sprintf('EMA%d', array_ch);
else
    error('Grid type "%s" not implemented yet.', array_type);
end
array_name = sprintf('Simulation_%s_', array_name);
if all(sim_DOA_deg == [0, 0])
    array_name = [array_name, 'SrcEar'];
else
    error('Parsing of the source position not implemented yet.');
end
fprintf('done.\n');

fprintf('Generating plot "%s" ... ', array_name);
if DO_PLOT
    fig = AKf;
    set(fig, 'NumberTitle', 'Off', 'Name', array_name);
    tl = tiledlayout(fig, 1, 1, 'TileSpacing', 'tight', 'Padding', 'tight');
    drawnow;
    nexttile(tl);
    plot_array_grid(array_rad, array_r);
    drawnow;
    if DO_PLOT_EXPORT
        fprintf('exporting ... ');
        [~, ~] = mkdir(plot_export_dir); % ignore warning if directory already exists
        file_name = fullfile(plot_export_dir, [array_name, '.pdf']);
        exportgraphics(fig, file_name, 'Append', false);
    end
    clear fig tl file_name;
    fprintf('done.\n');
else
    fprintf('skipped.\n'); %#ok<*UNRCH> 
end
fprintf('\n');

%%
% Helper function to convert directions from azimuth-elevation to azimuth-colatitude
azCo2azEl = @(dirs) [dirs(:, 1), pi/2 - dirs(:, 2)];

fprintf('Computing simulated array impulse responses ... ');
fprintf('for incidence direction [%.1f, %.1f] deg  ... ', sim_DOA_deg);
array_irs = simulateSphArray(sim_len, azCo2azEl(array_rad), deg2rad(sim_DOA_deg), ...
    sim_array_type, array_r, sim_order, sim_fs); % takes azimuth and elevation in rad
fprintf('done.\n');

fprintf('Generating plot "%s" ... ', array_name);
if DO_PLOT
    fig = AKf();
    set(fig, 'NumberTitle', 'Off', 'Name', array_name);
    tl = tiledlayout(fig, 2, 1, 'TileSpacing', 'tight', 'Padding', 'tight');
    drawnow;
    nexttile(tl);
    AKp(array_irs, 'et2d', 'du', 'n');
    drawnow;
    nexttile(tl);
    AKp(array_irs, 'm2d', 'fs', sim_fs, 'dr', 40);
    drawnow;
    if DO_PLOT_EXPORT
        fprintf('exporting ... ');
        [~, ~] = mkdir(plot_export_dir); % ignore warning if directory already exists
        file_name = fullfile(plot_export_dir, [array_name, '.pdf']);
        exportgraphics(fig, file_name, 'Append', true);
    end
    clear fig tl file_name;
    fprintf('done.\n');
else
    fprintf('skipped.\n');
end
fprintf('\n');

%%
fprintf('Determining plotting frequencies ... ');
array_alias_freq = sim_C * array_order / (2*pi * array_r);
fprintf('with aliasing above %.0f Hz ... ', array_alias_freq);
plot_freqs = ceil((array_alias_freq / 1e3) + 1) * 1e3; % above the aliasing frequency
plot_freqs = [plot_freqs/100, plot_freqs/10, plot_freqs/2, plot_freqs];
plot_freqs = arrayfun(@(freq) min([freq, sim_fs/2]), plot_freqs);
fprintf('at [%s] Hz ... ', AKlistIntegers(plot_freqs, false));
clear array_alias_freq;
fprintf('done.\n');

fprintf('Generating plot "%s" ... ', array_name);
if DO_PLOT
    fig = AKf();
    set(fig, 'NumberTitle', 'Off', 'Name', array_name);
    plot_sensor_signals(fig, array_irs, sim_fs, plot_freqs, array_rad, array_r);
    if DO_PLOT_EXPORT
        fprintf('exporting ... ');
        [~, ~] = mkdir(plot_export_dir); % ignore warning if directory already exists
        file_name = fullfile(plot_export_dir, [array_name, '.pdf']);
        exportgraphics(fig, file_name, 'Append', true);
    end
    clear fig file_name;
    fprintf('done.\n');
else
    fprintf('skipped.\n');
end
fprintf('\n');

%%

fprintf('Generating data file ')
if DO_FILE_EXPORT
    file_name = fullfile(data_export_dir, [array_name, '.sofa']);
    fprintf('"%s" ... ', file_name);

    % Generate SOFA struct
    sofa_struct = SOFAgetConventions('SingleRoomSRIR');
    sofa_struct.Data.IR = permute(array_irs, [3, 2, 1]);
    sofa_struct.API.N = sim_len;
    sofa_struct.API.R = array_ch;
    sofa_struct.Data.SamplingRate = double(sim_fs);
    sofa_struct.Data.Delay = zeros([1, array_ch]);

    sofa_struct.GLOBAL_ListenerShortName = sim_array_type;
    sofa_struct.ListenerPosition = [0, 0, 0];
    sofa_struct.ListenerPosition_Type = 'cartesian';
    sofa_struct.ListenerPosition_Units = 'metre';
    sofa_struct.ListenerView = [1, 0, 0];
    sofa_struct.ListenerView_Type = 'cartesian';
    sofa_struct.ListenerView_Units = 'metre';
    
    sofa_struct.GLOBAL_ReceiverShortName = 'Omnidirectional';
    sofa_struct.ReceiverPosition = [rad2deg(azCo2azEl(array_rad)), ...
        array_r * ones([array_ch, 1])]; % add radius in m
    sofa_struct.ReceiverPosition_Type = 'spherical';
    sofa_struct.ReceiverPosition_Units = 'degree, degree, metre';
    sofa_struct.ReceiverView = repmat([1, 0, 0], [array_ch, 1]);
    sofa_struct.ReceiverUp = repmat([0, 0, 1], [array_ch, 1]);

    sofa_struct.GLOBAL_SourceShortName = 'Plane wave';
    sofa_struct.SourcePosition = [1, 0, 0];
    sofa_struct.SourcePosition_Type = 'cartesian';
    sofa_struct.SourcePosition_Units = 'metre';
    sofa_struct.SourceView = [-1, 0, 0];
    sofa_struct.SourceView_Type = 'cartesian';
    sofa_struct.SourceView_Units = 'metre';

    % Remove unused non-mandatory SOFA data
    sofa_struct = SOFAremoveVariable(sofa_struct, 'RoomCornerA');
    sofa_struct = SOFAremoveVariable(sofa_struct, 'RoomCornerB');
    sofa_struct = SOFAremoveVariable(sofa_struct, 'RoomCorners_Type');
    sofa_struct = SOFAremoveVariable(sofa_struct, 'RoomCorners_Units');

    % Additional SOFA metadata
    sofa_struct.GLOBAL_Title = array_name;
    sofa_struct.GLOBAL_RoomType = 'anechoic';
    sofa_struct.GLOBAL_RoomDescription = 'Anechoic';
    sofa_struct.GLOBAL_Comment = 'Simulated with the Array Response Simulator toolbox'; 
    sofa_struct.GLOBAL_AuthorContact = 'hannes.helmholz@chalmers.se';
    sofa_struct.GLOBAL_Organization = 'Chalmers University of Technology';

    % Additional own metadata
    if ~any(isnan(array_weights)) && ~isempty(array_weights)
        sofa_struct = SOFAaddVariable(sofa_struct, ...
            'ReceiverQuadWeight', 'RI', array_weights);
    end

    fprintf('exporting ... ');
    [~, ~] = mkdir(data_export_dir); % ignore warning if directory already exists
    SOFAsave(file_name, sofa_struct);
    clear file_name;
    fprintf('done.\n');
else
    fprintf('... skipped.\n');
end
fprintf('\n');

%%
fprintf(' ... finished in %.0fh %.0fm %.0fs.\n', ...
    toc/3600, mod(toc,3600)/60, mod(toc,60));


%% helper functions
function plot_sensor_signals(fig, irs, fs, plot_freqs, grid_rad, array_r)
    tl = tiledlayout(fig, 'flow', 'TileSpacing', 'none', 'Padding', 'tight');
    colormap('hot');
    drawnow;
    
    grid_rad(:, 2) = pi/2 - grid_rad(:, 2); % colatitude to elevation
    [grid_cart(:, 1), grid_cart(:, 2), grid_cart(:, 3)] = ...
        sph2cart(grid_rad(:, 1), grid_rad(:, 2), array_r);

    specs_mag = db(abs(AKboth2singleSidedSpectrum(fft(irs))));
    freqs = linspace(0, fs/2, size(specs_mag, 1)).';
    [~, freq_ids] = min(abs(plot_freqs - freqs));
    clear freqs;
    
    for f = 1 : length(plot_freqs)
        spec = specs_mag(freq_ids(f), :).';
    
        ax(f) = nexttile(tl); %#ok<AGROW> 
        patch(ax(f), 'Faces', sphDelaunay(grid_rad), 'Vertices', grid_cart, ...
            'FaceVertexCData', spec, 'FaceColor', 'Interp', 'EdgeColor', 'none');
        if size(grid_cart, 1) < 1742
            hold on;
%             scatter3(ax(f), grid_cart(:, 1), grid_cart(:, 2), grid_cart(:, 3), ...
%                 (spec - min(spec) + 1) * 10, '*');
            scatter3(ax(f), grid_cart(:, 1), grid_cart(:, 2), grid_cart(:, 3), ...
                (spec - min(spec) + 5) * 5, 'ko', 'MarkerFaceColor', 'g');
        end
        title(ax(f), sprintf('f = %d Hz', plot_freqs(f)));
        view(3);
        axis equal;
        axis([-array_r, array_r, -array_r, array_r, -array_r, array_r]);
        axis off;
        drawnow;
    end; clear f spec;
    setappdata(fig, 'StoreTheLink', linkprop(ax, ...
        {'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim', 'CLim'}));
    caxis([min(specs_mag(freq_ids, :), [], 'all'), ...
        max(specs_mag(freq_ids, :), [], 'all')]);
    cb = colorbar;
    cb.Layout.Tile = 'East';
    cb.Label.String = 'Magnitude in dB';
    drawnow;
end
