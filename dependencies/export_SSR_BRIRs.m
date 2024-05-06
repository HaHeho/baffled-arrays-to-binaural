function export_SSR_BRIRs(out_file, sig, fs, ...
    az, is_name_append, is_norm_excess_amp)

    if nargin < 6; is_norm_excess_amp = false; end
    if nargin < 5 || isempty(is_name_append); is_name_append = true; end
    if nargin < 4 || isempty(az); az = 0 : 359; end % in deg

    WAV_BIT_DEPTH = 32;
%     WAV_BIT_DEPTH = 16;

    % expects signals with the size [samples, 360, 2]
    assert(size(sig, 2) == 360, 'Invalid data size (%d) provided.', size(sig, 2));
    assert(length(az) == 360, 'Invalid azimuth size (%d) provided.', length(az));

    if ismatrix(az)
        % Split off azimuth if elevation data is given as well
        % (discard elevation, since the azimuth data is all 0 and will not
        % be sorted)
%         el = az(2, :); 
        az = az(1, :);
    end

    % rearange data to desired range
    az = nav2sph(az);
    [~, idx] = sort(az);
    sig = sig(:, idx, :);

    % stack ears consecutively per azimuth
    sig = reshape(permute(sig, [1, 3, 2]), ...
        size(sig, 1), size(sig, 2) * size(sig, 3));

    out_path = fileparts(out_file);
    [~, ~] = mkdir(out_path); % ignore warning if directory already exists

    if is_name_append
        out_file = [out_file, '_SSR.wav'];
    end
    fprintf('"%s" ... ', out_file);
    
    % WAV output requires signals in range [-1 .. +1]
    sig_max = max(abs(sig), [], 'all');
    if sig_max > .99
        if is_norm_excess_amp
            sig_max = sig_max / .99;
            warning('normalizing amplitude by %.1f dB ... ', -db(sig_max));
            sig = sig / sig_max;
        elseif sig_max > 1
            error('Exceeding signal range of [-1 .. +1] by %.1f dB.', db(sig_max));
        end
    end

    audiowrite(out_file, sig, fs, 'BitsPerSample', WAV_BIT_DEPTH);
    fprintf('done.\n');
end
