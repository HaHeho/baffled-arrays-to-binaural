% Spatially subsample a high-resolution directional impulse response data 
% set into a different (lower-resolution) sampling grid in the spherical 
% harmonics domain. This is suitable for array and HRIR data sets.
%
% The script also performs a comparison of the subsampled data against a
% reference set if available.
% 
% -------------------------------------------------------------------------
%
% requires AKtools toolbox (run AKtoolsStart.m)
% $ svn checkout https://svn.ak.tu-berlin.de/svn/AKtools --username aktools --password ak
%
% requires Spherical-Harmonic-Transform toolbox
% $ git clone https://github.com/polarch/Spherical-Harmonic-Transform.git
% 
% -------------------------------------------------------------------------
%
% Hannes Helmholz 07.06.2023
%
% -------------------------------------------------------------------------
clear; clc; close all;

addpath(genpath('dependencies'));

%% configuration
% TH Cologne WDR data
[irs_order, irs_file] = deal(44, 'resources/ARIR_WDR/SIM_VSA_2702RS_PW_struct.mat');
[irs_order, irs_file] = deal(44, 'resources/ARIR_WDR/ANE_VSA_2702RS_struct.mat');
[irs_order, irs_file] = deal(29, 'resources/ARIR_WDR/CR1_VSA_SMA_LE1202_L_struct.mat');
[irs_order, irs_file] = deal(29, 'resources/ARIR_WDR/LBS_VSA_1202RS_PAC_struct.mat');
[target_orders, target_types] = deal([2, 4, 8, 12, 29], {'Lebedev', 'tDesign', 'Equatorial'});

% TH Cologne KU100 data
[irs_order, irs_file] = deal(44, 'resources/HRIR_KU100/HRIR_L2702.sofa');
[target_orders, target_types] = deal(44, {'SSR', 'vertical_SSR'});

sh_type = 'real'; % SH convention, or 'complex', results are identical
valid_frac_sm = 24; % in fractional octaves
plot_export_dir = 'plots/3_Subsampling'; % only for validation

global DO_FILE_EXPORT DO_PLOT_GRID DO_PLOT_EXPORT %#ok<*GVMIS> 
DO_FILE_EXPORT = true;
DO_PLOT_GRID   = true;
DO_PLOT_EXPORT = true;

% DO_FILE_EXPORT = false;
% DO_PLOT_GRID   = false;
% DO_PLOT_EXPORT = false;

%% load file and transform input data
tic; % start measuring execution time

fprintf('Loading file "%s" ... ', irs_file);
if endsWith(irs_file, '.mat')
    tmp = load(irs_file);
    irs = tmp.irChOne;
    irs_grid = [tmp.azimuth; tmp.colatitude].'; % azimut and colatitude in rad
    irs_weights = tmp.quadWeight.';
    irs_radius = tmp.radius;
elseif endsWith(irs_file, '.sofa')
    SOFAstart;
    tmp = SOFAload(irs_file);
    irs = get_SOFA_IRs(tmp);
    if contains(irs_file, {'SMA', 'XMA'}, 'IgnoreCase', true)
        irs_grid = deg2rad(tmp.ReceiverPosition(:, 1:2)); % in rad
        irs_radius = mean(tmp.ReceiverPosition(:, 3)); % in m
    else % HRIRs
        irs_grid = deg2rad(tmp.SourcePosition(:, 1:2)); % in rad
        irs_radius = mean(tmp.SourcePosition(:, 3)); % in m
    end
    irs_grid(:, 2) = pi/2 - irs_grid(:, 2); % elevation to colatitude
    if isfield(tmp, 'ReceiverQuadWeight')
        irs_weights = tmp.ReceiverQuadWeight;
    elseif isfield(tmp, 'SourceQuadWeight')
        irs_weights = tmp.SourceQuadWeight;
    else
        irs_weights = [];
    end
else
    error('Unknown data file format.');
end

% parse file names
[irs_path, irs_name, irs_end] = fileparts(irs_file);
irs_file = [irs_name, irs_end];
clear irs_name irs_end;
fprintf('done.\n');

if ~contains(irs_file, 'XMA', 'IgnoreCase', true)
    fprintf('Transforming into SH coefficients at N=%d ... ', irs_order);
    % normalize quadrature weights to 4*pi
    irs_weights = 4*pi / sum(irs_weights) * irs_weights;
    irs_SH = zeros(size(irs, 1), (irs_order + 1)^2, size(irs, 3));
    for e = 1 : size(irs, 3)
        if isempty(irs_weights)
            % transformation with Spherical-Harmonic-Transform and leasts-squares without quadrature weigths
            % -> identical to results with quadrature weights
            [irs_SH_SHT, ~] = leastSquaresSHT(irs_order, irs(:, :, e).', irs_grid, sh_type);
        else
            % transformation with Spherical-Harmonic-Transform and quadrature weigths
            [irs_SH_SHT, ~] = directSHT(irs_order, irs(:, :, e).', irs_grid, sh_type, irs_weights);
            % % transformation with Spherical-Harmonic-Transform and weighted least-squares
            % [irs_SH_SHT, ~] = leastSquaresSHT(irs_order, irs(:, :, e).', irs_grid, sh_type, irs_weights);
        end
        irs_SH(:, :, e) = irs_SH_SHT.';
        clear irs_SH_SHT;
    
        % % transformation with AKtools and with quadrature weigths
        % % -> different SH convention compared to Spherical-Harmonic-Transform
        % irs_SH_AK = AKsht(irs, false, [rad2deg(irs_grid), irs_weights], irs_order, 'complex', [], true, sh_type).';
    end; clear e;
    fprintf('done.\n');
end
fprintf('\n');

%% generete output data
for t = 1 : length(target_types)
    for o = 1 : length(target_orders)
        fprintf('Generating target %s grid ... ', target_types{t});
        clear target_grid target_weights;
        if ~contains(target_types{t}, 'SSR', 'IgnoreCase', true)
            fprintf('at N=%d ... ', target_orders(o));
        end
        if (any(strcmpi(target_types{t}, {'Fliege', 'Lebedev', 'tDesign'})) ...
                && target_orders(o) >= irs_order) || target_orders(o) > irs_order
            fprintf('skipped (given dataset SH order is equal or higher).\n\n');
            continue;
        end
        if strcmpi(target_types{t}, 'XMA')
            % TODO: This could surely be implemented in a smarter way
            IRS_XMA_CH = 18;
            assert(size(irs, 3) == IRS_XMA_CH, ...
                'This has only been implemented for datasets with %d channels.', IRS_XMA_CH);
            if target_orders(o) == 4
                target_ids = 1 : 2 : size(irs, 3);
            elseif target_orders(o) == 2
                target_ids = 1 : 3 : size(irs, 3);
            else
                error('This has not beem implemented for target order %d.', ...
                    target_orders(o));
            end
            target_irs = irs(:, :, target_ids);
            target_grid = irs_grid(target_ids, :);
            target_weights = [];
            clear target_ids;
        else
            [target_grid, target_weights] = get_array_grid( ...
                target_types{t}, target_orders(o), true, true);
        end
        target_ch = size(target_grid, 1);
        if strcmpi(target_types{t}, 'Fliege')
            target_ch_str = sprintf('SMA_FM%d', target_ch);
        elseif strcmpi(target_types{t}, 'Lebedev')
            target_ch_str = sprintf('SMA_LE%d', target_ch);
        elseif strcmpi(target_types{t}, 'tDesign')
            target_ch_str = sprintf('SMA_TD%d', target_ch);
        elseif any(strcmpi(target_types{t}, {'Eigenmike', 'Eigenmike32'}))
            target_ch_str = sprintf('SMA_EM%d', target_ch);
        elseif strcmpi(target_types{t}, 'Equatorial')
            target_ch_str = sprintf('EMA%d', target_ch);
        elseif strcmpi(target_types{t}, 'SSR')
            target_ch_str = sprintf('SSR%d', target_ch);
        elseif strcmpi(target_types{t}, 'vertical_SSR')
            target_ch_str = sprintf('vertical_SSR%d', target_ch);
        elseif ~strcmpi(target_types{t}, 'XMA')
            error('Grid type "%s" not implemented yet.', target_types{t});
        end
        fprintf('done.\n');
        
        if DO_PLOT_GRID && ~strcmpi(target_types{t}, 'XMA')
            fig_name = sprintf('Grid_%s', target_ch_str);
            file_name = fullfile(plot_export_dir, [fig_name, '.pdf']);
            
            fprintf('Generating plot "%s" ... ', fig_name);
            fig = AKf(20, 9);
            set(fig, 'NumberTitle', 'Off', 'Name', fig_name);
            tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'tight', 'Padding', 'tight');
            drawnow;
            nexttile(tl);
            plot_array_grid(target_grid, irs_radius);
            nexttile(tl);
            plot_array_grid(target_grid, irs_radius, false);
            if DO_PLOT_EXPORT
                fprintf('exporting ... ');
                [~, ~] = mkdir(plot_export_dir); % ignore warning if directory already exists
                exportgraphics(fig, file_name);
            end
            clear file_name file_name fig tl;
            fprintf('done.\n');
        end

        if ~strcmpi(target_types{t}, 'XMA')
            fprintf('Rendering to target grid ... ');
            target_irs = zeros(size(irs, 1), target_ch, size(irs, 3));
            for e = 1 : size(irs, 3)
                % transformation with Spherical-Harmonic-Transform
                target_irs(:, :, e) = inverseSHT(irs_SH(:, :, e).', target_grid, sh_type).';
            
                % % transformation with AKtools
                % % (result identical to Spherical-Harmonic-Transform)
                % target_irs(:, :, e) = AKisht(irs_SH_AK.', false, rad2deg(target_grid), 'complex', [], true, sh_type);
            end; clear e;
            fprintf('done.\n');
        end
        fprintf('\n');

        %%
        if contains(irs_file, {'SMA', 'VSA'}, 'IgnoreCase', true)
            target_config = strsplit(irs_file, '_SMA_LE2702_');
            if length(target_config) < 2
                target_config = strsplit(irs_file, '_SMA_LE1202_');
            end
            if length(target_config) < 2
                target_config = strsplit(irs_file, '_SMA_FM900_');
            end
            if length(target_config) < 2
                target_config = strsplit(irs_file, '_2702RS_'); % for THK WDR data
            end
            if length(target_config) < 2
                target_config = strsplit(irs_file, '_1202RS_'); % for THK WDR data
            end
            [~, target_file, target_end] = fileparts(target_config{2});
            reference_file = fullfile(irs_path, sprintf('%s_%s_%s%s', ... 
                target_config{1}, target_ch_str, target_file, target_end)); % before target_file!
            target_file = fullfile(irs_path, sprintf('%s_%s_%s_SubSampled_SH%d%s', ...
                target_config{1}, target_ch_str, target_file, irs_order, target_end));
        else
            [~, target_file, target_end] = fileparts(irs_file);
            if contains(irs_file, 'XMA', 'IgnoreCase', true)
                target_file = fullfile(irs_path, [regexprep(target_file, ...
                    'XMA\d+', sprintf('XMA%d', target_ch)), target_end]);
            else % HRIRs
                if contains(target_types{t}, 'SSR', 'IgnoreCase', true)
                    target_file = fullfile(irs_path, sprintf('%s_%s%s', ...
                        target_file, target_ch_str, target_end));
                else
                    target_file = fullfile(irs_path, sprintf('%s_N%d%s', ...
                        target_file, target_orders(o), target_end));
                end
            end
            reference_file = '';
        end
        clear target_ch_str target_config target_end;

        fprintf('Gathering target data ... ');
        target_data = tmp;
        if endsWith(target_file, '.mat')
            target_data.azimuth = target_grid(:, 1).';
            target_data.colatitude = target_grid(:, 2).';
            target_data.comments = sprintf(['Subsampled from "%s" with ', ...
                'Archontis Politis` Matlab toolbox'], irs_file); 
            target_data.nIr = target_ch;
            target_data.quadGrid = sprintf('%s %dSP', target_types{t}, target_ch);
            target_data.quadWeight = target_weights.';
            target_data.irChOne = target_irs;
            fprintf('done.\n');
            
            fprintf('Exporting file "%s" ... ', target_file);
            if DO_FILE_EXPORT
                save(target_file, '-struct', 'target_data', '-v7'); %#ok<*UNRCH> 
                fprintf('done.\n');
            else
                fprintf('skipped.\n');
            end
        
        elseif endsWith(target_file, '.sofa')
            if contains(irs_file, 'SMA', 'IgnoreCase', true)
                target_data.Data.IR = shiftdim(target_irs.', -1);
            else
                target_data.Data.IR = permute(target_irs, [2, 3, 1]);
            end
            if contains(irs_file, {'SMA', 'XMA'}, 'IgnoreCase', true)
                target_data.Data.Delay = zeros(1, target_ch);
                target_data.ReceiverPosition = rad2deg(target_grid); % in deg
                target_data.ReceiverPosition(:, 2) = ...
                    90 - target_data.ReceiverPosition(:, 2); % colatitude to elevation
                target_data.ReceiverPosition(:, 3) = irs_radius; % in m
                target_data.ReceiverView = target_data.ReceiverView(1:target_ch, :);
                target_data.ReceiverUp = target_data.ReceiverUp(1:target_ch, :);
            else % for HRIRs
                target_data.SourcePosition = rad2deg(target_grid); % in deg
                if contains(target_types{t}, 'SSR', 'IgnoreCase', true)
                    % revert incidence direction from source position
                    target_data.SourcePosition(:, 1) = -target_data.SourcePosition(:, 1);
                end
                target_data.SourcePosition(:, 2) = ...
                    90 - target_data.SourcePosition(:, 2); % colatitude to elevation
                target_data.SourcePosition(:, 3) = irs_radius; % in m
            end

            % Additional own metadata
            if ~isempty(target_weights) && ~any(isnan(target_weights))
                if contains(irs_file, 'SMA', 'IgnoreCase', true)
                    target_data = SOFAaddVariable(target_data, ...
                        'ReceiverQuadWeight', 'RI', target_weights);
                else % for HRIRs
                    target_data = SOFAaddVariable(target_data, ...
                        'SourceQuadWeight', 'MI', target_weights);
                end
            else
                target_data = SOFAremoveVariable(target_data, 'SourceQuadWeight');
                target_data = SOFAremoveVariable(target_data, 'ReceiverQuadWeight');
            end
            if isfield(target_data, 'AirTemperature') && ~strcmpi(target_types{t}, 'XMA')
                target_data.AirTemperature = repmat(...
                    mean(target_data.AirTemperature), [target_ch, 1]);
            end
            if isfield(target_data, 'AirHumidity') && ~strcmpi(target_types{t}, 'XMA')
                target_data.AirHumidity = repmat(...
                    mean(target_data.AirHumidity), [target_ch, 1]);
            end
            target_data = SOFAupdateDimensions(target_data);
        
            % Additional documentaiton
            target_data.GLOBAL_Comment = sprintf(['%s; Subsampled from "%s" with ', ...
                'Archontis Politis` Matlab toolbox'], target_data.GLOBAL_Comment, irs_file);
            fprintf('done.\n');
        
            fprintf('Exporting file "%s" ... ', target_file);
            if DO_FILE_EXPORT
                [~] = SOFAsave(target_file, target_data);
                fprintf('done.\n');
            else
                fprintf('skipped.\n');
            end
        else
            error('Unknown target file format for "%s".', target_file);
        end

        if DO_FILE_EXPORT && contains(target_types{t}, 'SSR', 'IgnoreCase', true)
            fprintf('Exporting file ');
            [ssr_dir, ssr_file, ~] = fileparts(target_file);
            ssr_file = fullfile(ssr_dir, [ssr_file, '.wav']);
            if endsWith(target_file, '.mat')
                fs = double(target_data.fs);
            else % if endsWith(target_file, '.sofa')
                fs = target_data.Data.SamplingRate;
            end
            export_SSR_BRIRs(ssr_file, target_irs, fs, [], false);
            clear ssr_dir ssr_file fs;
        end

        %% Validate data
        if ~exist('reference_file', 'var') || isempty(reference_file)
            fprintf('Skipped validation against reference.\n\n');
        elseif ~isfile(reference_file)
            fprintf('Skipped validation ... reference file "%s" not found.\n\n', reference_file);
        else
            fprintf('Performing validation against reference ... \n');
            close all;
            clear target_irs target_grid target_weights fs;
        
            if DO_FILE_EXPORT
                fprintf('Loading file "%s" ... ', target_file);
                clear target_data;
            else
                fprintf('Gathering target data ... ');
            end
            if endsWith(target_file, '.mat')
                if DO_FILE_EXPORT
                    target_data = load(target_file);
                end
                target_irs = target_data.irChOne;
                target_grid = rad2deg([target_data.azimuth; target_data.colatitude].'); % in deg
                target_grid(:, 2) = 90 - target_grid(:, 2); % colatitude to elevation
                target_weights = target_data.quadWeight.';
                fs = double(target_data.fs);
                
            elseif endsWith(target_file, '.sofa')
                if DO_FILE_EXPORT
                    target_data = SOFAload(target_file);
                end
                target_irs = get_SOFA_IRs(target_data);
                target_grid = target_data.ReceiverPosition(:, 1:2); % in deg
                target_weights = target_data.ReceiverQuadWeight;
                fs = target_data.Data.SamplingRate;
            else
                error('Unknown data file format.');
            end
            [~, target_file, ~] = fileparts(target_file);
            target_grid(target_grid(:, 1) < 0) = ...
                360 + target_grid(target_grid(:, 1) < 0); % wrap negative values
            fprintf('done.\n');
        
            fprintf('Loading file "%s" ... ', reference_file);
            if endsWith(reference_file, '.mat')
                reference_data = load(reference_file);
                reference_irs = reference_data.irChOne;
                reference_grid = rad2deg([reference_data.azimuth; reference_data.colatitude].'); % in deg
                reference_grid(:, 2) = 90 - reference_grid(:, 2); % colatitude to elevation
                reference_weights = reference_data.quadWeight.';
                if fs ~= double(reference_data.fs)
                    error('Mismatch in sampling frequncies.');
                end
            elseif endsWith(reference_file, '.sofa')
                reference_data = SOFAload(reference_file);
                reference_irs = get_SOFA_IRs(reference_data);
                reference_grid = reference_data.ReceiverPosition(:, 1:2); % in deg
                reference_weights = reference_data.ReceiverQuadWeight;
                if fs ~= reference_data.Data.SamplingRate
                    error('Mismatch in sampling frequncies.');
                end
            else
                error('Unknown data file format.');
            end
            reference_grid(reference_grid(:, 1) < 0) = ...
                360 + reference_grid(reference_grid(:, 1) < 0); % wrap negative values
            [~, reference_file, ~] = fileparts(reference_file);
            fprintf('done.\n');

            fprintf('Normalizing target quadrature weights ... ');
            if sum(target_weights) == sum(reference_weights)
                fprintf('skipped (already mathed).\n');
            else
                fprintf('from sum %.1f to sum %.1f ... ', sum(target_weights), sum(reference_weights));
                target_weights = target_weights / sum(target_weights) * sum(reference_weights);
                fprintf('done.\n');
            end

            fprintf('Aligning target grid ... ');
            [target_grid_cart(:, 1), target_grid_cart(:, 2), target_grid_cart(:, 3)] = ...
                sph2cart(target_grid(:, 1), target_grid(:, 2), 1);
            [reference_grid_cart(:, 1), reference_grid_cart(:, 2), reference_grid_cart(:, 3)] = ...
                sph2cart(reference_grid(:, 1), reference_grid(:, 2), 1);
            grid_ids = knnsearch(target_grid_cart, reference_grid_cart);
            target_irs = target_irs(:, grid_ids);
            target_grid = target_grid(grid_ids, :);
            target_weights = target_weights(grid_ids, :);
            clear target_grid_cart reference_grid_cart grid_ids;
            fprintf('done.\n');

            fprintf('Zero-padding IRs ... ');
            max_len = max(size(target_irs, 1), size(reference_irs, 1));
            fprintf('to %d samples ... ', max_len);
            target_irs(end:max_len, :) = 0;
            reference_irs(end:max_len, :) = 0;
            if size(target_irs, 2) ~= size(reference_irs, 2)
                error('Mismatch in data sizes.');
            end
            clear max_len;
            fprintf('done.\n');

            fprintf('Normalizing IR amplitude ... ');
            norm_fact = mean(rms(reference_irs) ./ rms(target_irs));
            fprintf('by %.1f dB of mean RMS amplitude ... ', db(norm_fact));
            target_irs = target_irs * norm_fact;
            clear norm_fact;
            fprintf('done.\n');
        
            fprintf('Calculating time difference ... ');
            diff_irs = abs(target_irs) ./ abs(reference_irs);
            fprintf('done.\n');
        
            fprintf('Calculating spectral difference ... ');
            target_spec = AKboth2singleSidedSpectrum(fft(target_irs));
            reference_spec = AKboth2singleSidedSpectrum(fft(reference_irs));
            % Apply smoothing before calculating difference
            target_spec = AKfractOctSmooth(target_spec, 'amp', fs, valid_frac_sm);
            reference_spec = AKfractOctSmooth(reference_spec, 'amp', fs, valid_frac_sm);
            diff_spec = compute_spectral_difference(target_spec, reference_spec);
%             % Apply smoothing after calculating difference
%             diff_spec = compute_spectral_difference(target_spec, reference_spec);
%             diff_spec = AKfractOctSmooth(diff_spec, 'amp', fs, valid_frac_sm);
            clear target_spec reference_spec;
            fprintf('done.\n');
            
            %% generate plot
            fig_name = target_file;
            file_name = fullfile(plot_export_dir, [fig_name, '.pdf']);
            
            plot_freq_rng = [40, 20e3];
            freqs = (0 : size(diff_spec, 1)-1)' * fs / size(diff_spec, 1) / 2;

            fprintf('Generating plot "%s" ... ', fig_name);
            fig = AKf();
            set(fig, 'NumberTitle', 'Off', 'Name', fig_name);
            tl = tiledlayout(fig, 2 + ~isempty(target_weights), 2, ...
                'TileSpacing', 'tight', 'Padding', 'tight');
            lgd_str = {reference_file, target_file};
            drawnow;
        
            nexttile(tl);
            plot(target_grid(:, 1), 'LineWidth', 5);
            hold on;
            plot(reference_grid(:, 1), 'LineWidth', 2);
            lgd = legend(lgd_str, 'Interpreter', 'none', 'Location', 'SouthEast');
            title(lgd, 'Azimuth');
            grid on; box on; axis tight;
            yticks(0:45:360); ylim([0, 360]);
            drawnow;

            nexttile(tl);
            plot(target_grid(:, 1) - reference_grid(:, 1), 'k', 'LineWidth', 2);
            lgd = legend('difference', 'Interpreter', 'none', 'Location', 'SouthWest');
            title(lgd, 'Azimuth');
            grid on; box on; axis tight;
            drawnow;

            nexttile(tl);
            plot(target_grid(:, 2), 'LineWidth', 5);
            hold on;
            plot(reference_grid(:, 2), 'LineWidth', 2);
            lgd = legend(lgd_str, 'Interpreter', 'none', 'Location', 'NorthEast');
            title(lgd, 'Elevation');
            grid on; box on; axis tight;
            yticks(-90:45:90); ylim([-90, 90]);
            drawnow;

            nexttile(tl);
            plot(target_grid(:, 2) - reference_grid(:, 2), 'k', 'LineWidth', 2);
            lgd = legend('difference', 'Interpreter', 'none', 'Location', 'NorthWest');
            title(lgd, 'Elevation');
            grid on; box on; axis tight;
            ylabel(tl, 'Angle in degree');
            drawnow;

            if ~isempty(target_weights)
                nexttile(tl);
                plot(target_weights, 'LineWidth', 5);
                hold on;
                plot(reference_weights, 'LineWidth', 2);
                lgd = legend(lgd_str, 'Interpreter', 'none', 'Location', 'SouthEast');
                title(lgd, 'Quadrature weight');
                grid on; box on; axis tight;
                drawnow;

                nexttile(tl);
                plot(target_weights - reference_weights, 'k', 'LineWidth', 2);
                lgd = legend('difference', 'Interpreter', 'none', 'Location', 'SouthWest');
                title(lgd, 'Quadrature weight');
                grid on; box on; axis tight;
            end
            xlabel(tl, 'Position');
            drawnow;

            if DO_PLOT_EXPORT
                fprintf('exporting ... ');
                [~, ~] = mkdir(plot_export_dir); % ignore warning if directory already exists
                exportgraphics(fig, file_name);
            end
        
            clear fig lgd_str lgd tl;
            fprintf('done.\n');
        
            %% generate plot
            fprintf('Generating plot "%s" ... ', fig_name);
            fig = AKf();
            set(fig, 'NumberTitle', 'Off', 'Name', fig_name);
            tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'tight', 'Padding', 'tight', ...
                'TileIndexing', 'columnmajor');
            lgd_str = {reference_file, target_file};
            drawnow;
        
            nexttile(tl);
            ch = 1;
            title_str = sprintf('Channel %d,  az = %.0f deg,  el = %.0f deg', ch, ...
                reference_grid(ch, 1), reference_grid(ch, 2));
            AKp([reference_irs(:, ch), target_irs(:, ch)], 'et2d', 'xu', 'n');
            lgd = legend(lgd_str, 'Interpreter', 'None', 'Location', 'NorthEast');
            title(lgd, title_str);
            drawnow;
            nexttile(tl);
            AKp([reference_irs(:, ch), target_irs(:, ch)], 'm2d', 'fs', fs);
            lgd = legend(lgd_str, 'Interpreter', 'None', 'Location', 'SouthWest');
            title(lgd, title_str);
            drawnow;
        
            nexttile(tl);
            ch = size(target_irs, 2);
            title_str = sprintf('Channel %d,  az = %.0f deg,  el = %.0f deg', ch, ...
                reference_grid(ch, 1), reference_grid(ch, 2));
            AKp([reference_irs(:, ch), target_irs(:, ch)], 'et2d', 'xu', 'n');
            lgd = legend(lgd_str, 'Interpreter', 'None', 'Location', 'NorthEast');
            title(lgd, title_str);
            drawnow;
            nexttile(tl);
            AKp([reference_irs(:, ch), target_irs(:, ch)], 'm2d', 'fs', fs);
            lgd = legend(lgd_str, 'Interpreter', 'None', 'Location', 'SouthWest');
            title(lgd, title_str);
            drawnow;
        
            if DO_PLOT_EXPORT
                fprintf('exporting ... ');
                exportgraphics(fig, file_name, 'Append', true);
            end
        
            clear lgd_str title_str fig tl ch lgd;
            fprintf('done.\n');
        
            %% generate plot
            fprintf('Generating plot "%s" ... ', fig_name);
            fig = AKf();
            set(fig, 'NumberTitle', 'Off', 'Name', fig_name);
            tl = tiledlayout(fig, 2, 1, 'TileSpacing', 'tight', 'Padding', 'tight');
            drawnow;
        
            sensors = 1 : size(diff_spec, 2);
            times = (0 : size(diff_irs, 1)-1)' * 1000 / fs;
        
            ax = nexttile(tl);
            title_str = sprintf('%s   vs.   %s', fig_name, reference_file);
            colormap(ax, AKcolormaps('Greys', 128, true));
            plot_ETC_azims_IRs(sensors, times, diff_irs, title_str);
            set(ax, 'XDir', 'Normal');
            xlabel(ax, 'Position');
            title(ax, title_str);
            drawnow;

            ax = nexttile(tl);
            if valid_frac_sm
                title_str = sprintf('%s   (1/%d oct. smoothing)', title_str, valid_frac_sm);
            end
            colormap(ax, AKcolormaps('RdBu', 128, true));
            plot_spec_azims(sensors, freqs, diff_spec, title_str, 1, plot_freq_rng);
            set(ax, 'XDir', 'Normal');
            set(ax.Children(1), 'FaceColor', 'Flat');
            caxis(ax, [-15, 15]);
            xlabel(ax, 'Position');
            title(ax, title_str);
            drawnow;

            if DO_PLOT_EXPORT
                fprintf('exporting ... ');
                exportgraphics(fig, file_name, 'Append', true);
            end
            
            clear fig tl sensors times ax;
            fprintf('done.\n');

            %% generate plot
            fprintf('Generating plot "%s" ... ', fig_name);
            fig = AKf();
            set(fig, 'NumberTitle', 'Off', 'Name', fig_name);
            tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'tight', 'Padding', 'tight');
            drawnow;

            r = 1;
            [pos(:, 1), pos(:, 2), pos(:, 3)] = sph2cart(...
                deg2rad(reference_grid(:, 1)), deg2rad(reference_grid(:, 2)), r);
            diff_freqs = freqs >= plot_freq_rng(1) & freqs <= plot_freq_rng(2);

            diff_spec_mean = mean(db(abs(diff_spec(diff_freqs, :))), 1);  % from smoothed spectra
            c_lim = max(abs(diff_spec_mean));

            nexttile(tl);
            quiver3(0, 0, 0, r * 0.75, 0, 0, 'off', 'filled', ...
                'LineWidth', 2, 'Color', 'k', 'MaxHeadSize', 3);
            hold on;
            scatter3(pos(:, 1), pos(:, 2), pos(:, 3), 3e3 / sqrt(size(pos, 1)), ...
                diff_spec_mean, 'LineWidth', 50 / sqrt(size(pos, 1)));
            grid on; box on;
            axis equal;
            axis([-r, r, -r, r, -r, r]);
            caxis([-c_lim, c_lim]);
            set(gca, 'Color', [0.6, 0.6, 0.6]);
            cb = colorbar('EastOutside');
            cb.Label.String = 'Mean spectral magnitude difference in dB';
            colormap(gca, AKcolormaps('RdBu', 128, true));
            drawnow;

            diff_spec_mean = mean(abs(db(abs(diff_spec(diff_freqs, :)))), 1);  % from smoothed spectra
            c_lim = max(abs(diff_spec_mean));

            nexttile(tl);
            quiver3(0, 0, 0, r * 0.75, 0, 0, 'off', 'filled', ...
                'LineWidth', 2, 'Color', 'k', 'MaxHeadSize', 3);
            hold on;
            scatter3(pos(:, 1), pos(:, 2), pos(:, 3), 3e3 / sqrt(size(pos, 1)), ...
                diff_spec_mean, 'LineWidth', 50 / sqrt(size(pos, 1)));
            grid on; box on;
            axis equal;
            axis([-r, r, -r, r, -r, r]);
            caxis([0, c_lim]);
            set(gca, 'Color', [0.6, 0.6, 0.6]);
            cb = colorbar('EastOutside');
            cb.Label.String = 'Mean absolute spectral magnitude difference in dB';
            colormap(gca, AKcolormaps('Greens', 128, false));
            title(tl, title_str, 'Interpreter', 'None');
            drawnow;

            if DO_PLOT_EXPORT
                fprintf('exporting ... ');
                exportgraphics(fig, file_name, 'Append', true);
            end

            clear fig tl cb c_lim diff_freqs diff_spec_mean r pos;
            fprintf('done.\n');

            %%
            clear fig_name file_name freqs plot_freq_rng title_str;
            fprintf('\n');
        end
    end
end; clear t o;

fprintf(' ... finished in %.0fh %.0fm %.0fs.\n', ...
    toc/3600, mod(toc,3600)/60, mod(toc,60));
