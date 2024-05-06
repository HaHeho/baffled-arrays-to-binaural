function [hrtfs, hrtfs_nm] = render_sh_coefficients_binaurally( ...
    hrir_file, S_breve, head_orientation, N, sphharm_type, is_with_weigths)

    if nargin < 6; is_with_weigths = false; end
    if nargin < 5; sphharm_type = 'complex'; end

    WEIGHTS_FILE = 'VSA_2702_leb_grid.mat';
%     WEIGHTS_FILE = ''; % do not use weights

    block_length = (size(S_breve, 1) - 1) * 2;
    
    % --------------------------- prepare HRTFs -------------------------------
    [file_path, file_name, ~] = fileparts(hrir_file);
    hrir_weights_file = fullfile(file_path, WEIGHTS_FILE);
    hrir_coeffs_file_str = fullfile(file_path, sprintf('%s_coeffs_%s_%%s_%d.mat', ...
        file_name, sphharm_type, block_length/2));
    if ~strcmp(filesep, '/') % required on Windows
        hrir_coeffs_file_str = strrep(hrir_coeffs_file_str, filesep, '/');
    end
    hrir_coeffs_withWeights_file = sprintf(hrir_coeffs_file_str, 'withWeights');
    hrir_coeffs_noWeights_file = sprintf(hrir_coeffs_file_str, sprintf('noWeights_sh%d', N));

    hrtfs_nm = [];
    if is_with_weigths && isfile(hrir_coeffs_withWeights_file)
        fprintf('Loading pre-comuputed SH HRIRs file "%s" ... ', hrir_coeffs_withWeights_file);
        hrtfs_data = load(hrir_coeffs_withWeights_file);
        hrtfs_nm = hrtfs_data.hrtfs_nm;
        if size(hrtfs_nm, 2) < (N+1)^2 ...
                || (size(hrtfs_nm, 2) > (N+1)^2 && ~hrtfs_data.used_grid_weights)
            fprintf('available at N=%d only, so they will be re-computed ... ', ...
                sqrt(size(hrtfs_nm, 2)) - 1);
            hrtfs_nm = [];
        else
            fprintf('truncated to N=%d ... ', N);
            hrtfs_nm = hrtfs_nm(:, 1:(N+1)^2, :);
        end
        fprintf('done.\n');
    end
    if isempty(hrtfs_nm) && isfile(hrir_coeffs_noWeights_file)
        fprintf('Loading pre-comuputed SH HRIRs file "%s" ... ', hrir_coeffs_noWeights_file);
        hrtfs_data = load(hrir_coeffs_noWeights_file);
        hrtfs_nm = hrtfs_data.hrtfs_nm;
        fprintf('done.\n');
    end
    clear hrtfs_data;

    if isempty(hrtfs_nm)
        hrtfs_nm = compute_hrtf_sh_coefficients(hrir_file, block_length, N, sphharm_type, ...
            hrir_coeffs_noWeights_file, hrir_coeffs_withWeights_file, hrir_weights_file, is_with_weigths);
    end

    if any(S_breve)
        % ------------------------------ rendering --------------------------------
        hrirs = render_sh_coefficients_binaurally_blockwise( ...
            hrtfs_nm, S_breve, head_orientation, N, sphharm_type);
        hrtfs = AKboth2singleSidedSpectrum(fft(hrirs));
    else
        hrtfs = NaN;
    end
end

%% helper functions
function hrtfs_nm = compute_hrtf_sh_coefficients(hrir_file, block_length, N, sphharm_type, ...
        hrir_coeffs_noWeights_file, hrir_coeffs_withWeights_file, hrir_weights_file, is_with_weigths)

    if endsWith(hrir_file, '.sofa')
        fprintf('Loading HRIR file "%s" ... ', hrir_file);
        hrir_sofa = SOFAload(hrir_file);
    
        % extract measurement directions
        grid_azi = deg2rad(     hrir_sofa.SourcePosition(:, 1)).';
        grid_col = deg2rad(90 - hrir_sofa.SourcePosition(:, 2)).'; % elevation to colatitude
            
        % extract IRs
        hrirs = get_SOFA_IRs(hrir_sofa);
        fprintf('done.\n');

        if is_with_weigths
            try
                fprintf('Loading quadrature weights from SOFA file ... ');
                grid_weights = hrir_sofa.SourceQuadWeight.';
                fprintf('done.\n');
            catch
                fprintf('skipped (field not found).\n');
                
                try
                    fprintf('Loading quadrature weights file "%s" ... ', hrir_weights_file);
                    tmp = load(hrir_weights_file);
                    grid_weights = tmp.grid_weights;
                    clear tmp;
                    fprintf('done.\n');
        
                    if length(grid_weights) ~= length(grid_azi)
                        error('Mismatch in HRTF sampling grid size and loaded number of quadrature weights.');
                    end
                catch
                    fprintf('skipped (file not found).\n');
                    grid_weights = [];
                end
            end
            clear hrir_sofa;
        else
            fprintf('Loading quadrature weights ... skipped (should not be used).\n');
            grid_weights = [];
        end

    elseif endsWith(hrir_file, '.mat')
        error('Not implemented yet.')
    else
        error('Unknown file type.');
    end

    % zero pad
    hrirs(end+1:block_length, :, :) = 0;

    % fft and remove redundant information
    hrtfs = AKboth2singleSidedSpectrum(fft(hrirs));

    fprintf('Transforming into SHs ... at N=%d ... ', N);
    hrtfs_nm = zeros(size(hrtfs, 1), (N+1)^2, size(hrtfs, 3), 'like', hrtfs);
    for e = 1 : size(hrtfs, 3) % each ear
        % Do not convert source position to propagation direction
        hrtfs_nm(:, :, e) = get_sound_field_sh_coeffs_sma( ...
            hrtfs(:, :, e), [], N, [grid_azi; grid_col], grid_weights, sphharm_type);
    end

    used_grid_weights = ~isempty(grid_weights);
    if used_grid_weights
        hrir_coeffs_file = hrir_coeffs_withWeights_file;
    else
        hrir_coeffs_file = hrir_coeffs_noWeights_file;
    end
    fprintf('exporting as "%s" ... ', hrir_coeffs_file);
    save(hrir_coeffs_file, 'hrtfs_nm', 'used_grid_weights');
    fprintf('done.\n');
end
