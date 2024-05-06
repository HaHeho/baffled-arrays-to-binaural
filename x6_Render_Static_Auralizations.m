% Perform convolution of (rendered) binaural room impulse responses with a
% source audio signal. This is done for a specified selection of static
% head orientations as well as continous rotation over all horizontal head
% orientations.
%
% -------------------------------------------------------------------------
%
% Hannes Helmholz 17.02.2023
%
% -------------------------------------------------------------------------
clear; clc; close all;

addpath(genpath('dependencies'));

%% configuration
[audio_file, audio_length, audio_fade] = deal('resources/Drums2_48.wav', 2.2, 0.4); % length and fade in s

data_exclude_configs = {'_MagLS', 'inCH', 'wDC'};
res_export_angles = [0, 90]; % in deg, exported head orientations
res_export_rotation = true;
res_export_dir = 'resources/BRIR_auralization/';
hpeq_file = '';

global DO_EXPORT_FILE %#ok<*GVMIS>
DO_EXPORT_FILE = true;
% DO_EXPORT_FILE = false;

%% 
tic; % start measuring execution time

render_static_auralizations({'resources/BRIR_rendered/HRIR_L2702/SIM*1202*', ...
    'resources/BRIR_rendered/HRIR_L2702/SIM*_EMA5_*', 'resources/BRIR_rendered/HRIR_L2702/SIM*_TD14_*', ...
    'resources/BRIR_rendered/HRIR_L2702/SIM*_EMA17_*', 'resources/BRIR_rendered/HRIR_L2702/SIM*_TD146_*'}, ...
    data_exclude_configs, audio_file, audio_length, audio_fade, ...
    res_export_angles, res_export_rotation, -6.0, res_export_dir, hpeq_file); % -6.4 dB matched to not yield clipping of exported signal
render_static_auralizations({'resources/BRIR_rendered/HRIR_L2702/SIM*_EM32_*'}, ...
    data_exclude_configs, audio_file, audio_length, audio_fade, ...
    res_export_angles, res_export_rotation, -8.6, res_export_dir, hpeq_file); % -6.4 dB matched to not yield clipping of exported signal
render_static_auralizations({'resources/BRIR_rendered/HRIR_L2702/LBS*1202*', ...
    'resources/BRIR_rendered/HRIR_L2702/LBS*_EMA5_*', 'resources/BRIR_rendered/HRIR_L2702/LBS*_TD14_*', ...
    'resources/BRIR_rendered/HRIR_L2702/LBS*_EMA17_*', 'resources/BRIR_rendered/HRIR_L2702/LBS*_TD146_*'}, ...
    data_exclude_configs, audio_file, audio_length, audio_fade, ...
    res_export_angles, res_export_rotation, -11.7, res_export_dir, hpeq_file); % -11.7 dB matched to not yield clipping of exported signal

fprintf(' ... finished in %.0fh %.0fm %.0fs.\n', ...
    toc/3600, mod(toc,3600)/60, mod(toc,60));


%% helper functions
function render_static_auralizations(data_base_dirs, data_exclude_configs, ...
    audio_file, audio_length, audio_fade, ...
    res_export_angles, res_export_rotation, res_export_amp, res_export_dir, hpeq_file)
    
    global DO_EXPORT_FILE
    DATA_END = '_SSR.wav';
    DO_INVERT_ANGLES = true; % 90 deg is left (horizontal) or above (vertical)
    % DO_INVERT_ANGLES = false; % 90 deg is right (horizontal) or below (vertical)
    ROTATION_BLOCK_SIZE = 128; % in samples

    fprintf('Loading audio file "%s" ... ', audio_file);
    [audio_data, audio_fs] = audioread(audio_file);
    if size(audio_data, 2) > 1
        fprintf('discarding %d channels ... ', size(audio_data, 2) - 1);
        audio_data = audio_data(:, 1);
    end
    fprintf('truncating to %.1f s length ... ', audio_length);
    audio_data = audio_data(1:round(audio_fs * audio_length));
    fprintf('appending %.1f s silence ... ', audio_fade);
    audio_data(end:end+round(audio_fs * audio_fade)) = 0;
    [~, audio_file, ~] = fileparts(audio_file);
    fprintf('done.\n');
    fprintf('\n'); 
    
    %%
    fprintf('Loading headphone equalization filter ');
    if ~isempty(hpeq_file)
        fprintf('"%s" ... ', hpeq_file);
        if endsWith(hpeq_file, '.mat')
            tmp = load(hpeq_file);
        else
            [~, ~, tmp] = fileparts(hpeq_file); 
            error('Data loading for data type "*%s" not implemented yet.', tmp);
        end
        hpcf_data = tmp.hpcf.linPhase; % linear-phase version
%         hpcf_data = tmp.hpcf.minPhase; % minimum-phase version
        assert(audio_fs == tmp.hpcf.fs, 'Mismatch in sampling frequencies.');
%         hpeq_name = strrep(tmp.hpcf.hpName, ' ', '_');
        [~, hpeq_name, ~] = fileparts(hpeq_file);
        clear tmp;
        fprintf('done.\n');
    
        fprintf('Applying headphone equalization filter ... ');
        % this makes the audio data two channels if the headphone compensation
        % filter is given for two ears
        audio_data = fftfilt(hpcf_data, audio_data);
        fprintf('done.\n');
    else
        fprintf('... skipped (no file provided).\n');
    end
    fprintf('\n'); 
    
    %%
    if ~iscell(data_base_dirs); data_base_dirs = {data_base_dirs}; end
    data_base_dirs = strip(data_base_dirs, 'right', '\');
    data_base_dirs = strip(data_base_dirs, 'right', '/');
    current_file = fileparts(mfilename('fullpath'));

    data_files = [];
    for d = 1 : length(data_base_dirs)
        fprintf('Parsing configuration files ... in "%s" ... ', data_base_dirs{d});
        new_files = dir([data_base_dirs{d}, DATA_END]);
        % exclude specified configurations
        new_files = new_files(~contains({new_files.name}, data_exclude_configs, 'IgnoreCase', true));
        data_files = [data_files; new_files]; %#ok<AGROW>
        fprintf('found %d files ... ', length(new_files));
        fprintf('done.\n');
    end; clear d new_files;
    fprintf('\n');
    
    %%
    export_level_max = -Inf;
    for f = 1 : length(data_files)
        data_file = fullfile(data_files(f).folder, data_files(f).name);
        data_file = data_file(length(current_file)+2:end); % remove absolute path

        % determine export file directory and name
        res_file = data_file(length(fileparts(data_file))+2:end);
        res_file = res_file(1:end-length(DATA_END)); % remove data ending

        [~, res_dir, ~] = fileparts(fileparts(data_base_dirs{1}));
        if ~isempty(hpeq_file)
            res_dir = sprintf('%s+%s', res_dir, hpeq_name); % add headphone name
        end
        res_file = fullfile(res_dir, res_file); % add subdirectory
        res_file = fullfile(res_export_dir, res_file);
        clear res_dir;
        
        fprintf('[%d/%d] ', f, length(data_files));
        fprintf('Loading file "%s" ... ', data_file);
        [data, fs] = audioread(data_file);
        assert(audio_fs == fs, 'Mismatch in sampling frequencies.');
        % change order of direction dimension to be last for convenience later
        data = permute(unstack_SSR_BRIRs(data), [1, 3, 2]);
        fprintf('done.\n');

        for a = 1 : length(res_export_angles) + res_export_rotation
            if a <= length(res_export_angles)
                % convolve the audio signal with the extracted head orientation
                % (this cannot handle negative angles at the moment)
                angle_id = round(res_export_angles(a)) + 1;
                if DO_INVERT_ANGLES; angle_id = size(data, 3) - angle_id; end
                res = fftfilt(data(:, :, angle_id), audio_data);
            else
                % ==============================================================
                % This version implements the rotation by partitioned
                % convolution with Matlab DSP objects.
                % ==============================================================
                % Set up DSP signal source
                audio_src = dsp.SignalSource(audio_data, ...
                    'SamplesPerFrame', ROTATION_BLOCK_SIZE, 'SignalEndAction', 'Cyclic repetition');

                % Initialize system objects
                data_FIR{1} = dsp.FIRFilter('NumeratorSource', 'Input port');
                data_FIR{2} = dsp.FIRFilter('NumeratorSource', 'Input port');
                data_FIR{3} = dsp.FIRFilter('NumeratorSource', 'Input port');
                data_FIR{4} = dsp.FIRFilter('NumeratorSource', 'Input port');

                % Fading window to blend BRIRs
                fade_win = hann(2 * ROTATION_BLOCK_SIZE);
                fade_win = reshape(fade_win, [ROTATION_BLOCK_SIZE, 2]);

                audio_len = size(audio_data, 1) * 2;
                data_angles = size(data, 3);
                angle_ids = round(linspace(1, data_angles, audio_len)).';
                if DO_INVERT_ANGLES
                    angle_ids = angle_ids - 1;
                    angle_ids = data_angles - angle_ids;
                end
                % Get BRIR for starting head orientation
                data_current = data(:, :, angle_ids(1)).';

                % Allocate result matrix
                blocks_num = ceil(audio_len / ROTATION_BLOCK_SIZE);
                res = zeros(ROTATION_BLOCK_SIZE * blocks_num, size(data, 2));
                angle_ids(end+1:size(res, 1)) = angle_ids(end);

                fprintf('Computing rotation ... ');
                for b = 0 : blocks_num - 1 % all blocks
                    b_ids = (1 : ROTATION_BLOCK_SIZE) + (b * ROTATION_BLOCK_SIZE);
                    angle_id = angle_ids(b_ids(1));
                    if mod(b, 20 * blocks_num / data_angles) < 1
                        fprintf('[%d]', angle_id);
                    end

                    % Get audio bloc from dsp source
                    audio_block = audio_src();
        
                    % Apply BRIR of previous head orientation
                    res_block(:, 1) = data_FIR{1}(audio_block, data_current(1, :)); % Left
                    res_block(:, 2) = data_FIR{2}(audio_block, data_current(2, :)); % Right
                    
                    % Get BRIR for current head orientation
                    data_current = data(:, :, angle_id).';
                
                    % Apply current BRIR and fade with previous
                    res_block(:, 1) = fade_win(:, 2) .* res_block(:, 1) ...
                        + fade_win(:, 1) .* data_FIR{3}(audio_block, data_current(1, :)); % Left
                    res_block(:, 2) = fade_win(:, 2) .* res_block(:, 2) ...
                        + fade_win(:, 1) .* data_FIR{4}(audio_block, data_current(2, :)); % Right

                    % Add result
                    res(b_ids, :) = res_block;
                end; clear b b_ids audio_block res_block;
                fprintf(' ... done.\n');

                clear fade_win audio_len data_angles angle_ids data_current blocks_num;
                release(audio_src);
                release(data_FIR{1});
                release(data_FIR{2});
                release(data_FIR{3});
                release(data_FIR{4});
                clear audio_src data_FIR;
                % ==============================================================

                % % ==============================================================
                % % This version implements the rotation by sample-wise weighted 
                % % addition of the BRIRs.
                % % ==============================================================
                % % Repeat the audio signal
                % audio_data = [audio_data; audio_data]; %#ok<AGROW>
                % 
                % audio_len = size(audio_data, 1);
                % [data_len, ~, data_angles] = size(data);
                % angle_ids = round(linspace(1, data_angles, audio_len)).';
                % if DO_INVERT_ANGLES
                %     angle_ids = angle_ids - 1;
                %     angle_ids = data_angles - angle_ids;
                % end
                % res = zeros(audio_len + data_len - 1, size(data, 2));
                % 
                % fprintf('Computing rotation ... ');
                % for s = 1 : audio_len
                %     angle_id = angle_ids(s);
                %     if mod(s, 10 * audio_len / data_angles) < 1
                %         fprintf('[%d]', angle_id);
                %     end
                %     % sample-wise weight with the audio signal
                %     s_ids = s : s+data_len-1;
                %     res(s_ids, :) = res(s_ids, :) + audio_data(s, :) .* data(:, :, angle_id);
                % end; clear s s_ids angle_id;
                % fprintf(' ... done.\n');
                % 
                % % restore the original audio signal length again
                % audio_data = audio_data(1 : ceil(size(audio_data, 1) / 2), :);
                % 
                % clear audio_len data_len data_angles angle_ids;
                % % ==============================================================
            end
            % apply amplification / attenuation to prevent clipping
            res = res * db2mag(res_export_amp);
    
            % determine export file name
            export_file = sprintf('%s_%s_%%sdeg.flac', res_file, audio_file);
            if a > length(res_export_angles)
                export_file = sprintf(export_file, '0to360');
            else
                export_file = sprintf(export_file, num2str(res_export_angles(a)));
            end
    
            fprintf('Exporting file "%s" ... ', export_file);
            res_level_max = max(db(abs(res)), [], 'all');
            export_level_max = max([export_level_max, res_level_max]);
            fprintf('with peak level of %+.1f dB ... ', res_level_max);
            if DO_EXPORT_FILE
                [~, ~] = mkdir(fileparts(export_file)); % ignore warning if directory already exists
                audiowrite(export_file, res, fs, 'BitsPerSample', 16);
                fprintf('done.\n');
            else
                fprintf('skipped.\n');
            end
        end; clear a angle_id res export_file res_level_max;
        fprintf('\n');
    end; clear f data_file res_file data fs;
    
    if ~isempty(data_files)
        fprintf(' ... export maxiumum level of %.1f dB.\n\n', export_level_max);
    end
end
