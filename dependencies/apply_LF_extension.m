function irs_ext = apply_LF_extension(irs, fs, ...
    lf_ref_freq_rng, target_len, target_start)
    % based on https://github.com/spatialaudio/lf-corrected-kemar-hrtfs

    if nargin < 5; target_start = 0; end
    if nargin < 4; target_len = 0; end
    if nargin < 3; lf_ref_freq_rng = [100, 300]; end % in Hz, from repository example
    assert(length(lf_ref_freq_rng) == 2);
    
    % IN_PAD_FACT = 10;
    IN_PAD_FACT = 1;

    if ~target_len
        target_len = size(irs, 1);
    elseif target_len > size(irs, 1)
        irs(end+1:target_len, :, :) = 0;
    end
    if target_start < 1
        target_start = 1;
    end

    % apply zero padding for better frequency resolution
    % (this may help with some direction dependent problems to achive
    % stable amplification at low frequencies)
    irs(end+1:size(irs, 1) * IN_PAD_FACT, :, :) = 0;

    freqs = linspace(0, fs/2, size(irs, 1)/2 + 1);
    % find indices for frequency range used to correct low frequency data
    freq_idx(1) = find(freqs >= lf_ref_freq_rng(1), 1, 'first');
    freq_idx(2) = find(freqs <= lf_ref_freq_rng(2), 1, 'last');

    irs_ext = zeros(target_len - target_start + 1, size(irs, 2), size(irs, 3));
    for d = 1 : size(irs_ext, 2) % each direction
        for e = 1 : size(irs_ext, 3) % each ear
            [spec, is_even] = AKboth2singleSidedSpectrum(fft(irs(:, d, e)));
            
            % correct magnitude
            % calculate magnitude RMS in given range
            mag_rms = rms(abs(spec(freq_idx(1):freq_idx(2))));
            % replace low frequency magnitude with mean value
            spec_mag_ext = abs(spec);
            spec_mag_ext(1:freq_idx(1)-1) = mag_rms;

            % correct phase
            spec_phase = unwrap(angle(spec));
            spec_phase_ext = spec_phase;
            % extrapolate phase from data in given range
            spec_phase_ext(1:freq_idx(1)-1) = ...
                interp1(freqs(freq_idx(1):freq_idx(2)), ...
                spec_phase(freq_idx(1):freq_idx(2)), freqs(1:freq_idx(1)-1), ...
                'linear', 'extrap');
            
            % compose complex spectrum
            spec_ext = spec_mag_ext .* exp(1i * spec_phase_ext);
            % calculate impulse response
            ir_ext = ifft(AKsingle2bothSidedSpectrum(spec_ext, is_even));

            % truncate corrected impulse response
            irs_ext(:, d, e) = ir_ext(target_start:target_len);
        end
    end

    % d = 143;
    % e = 1;
    % AKf();
    % AKp(irs(:, d, e), {'et2d', 'm2d'; '', ''}, 'fs', fs, 'xu', 'n');
    % AKp(irs_ext(:, d, e), {'et2d', 'm2d'; '', ''}, 'fs', fs, 'xu', 'n', 'N', target_len, 'c', 'r');
    % AKp(irs(:, d, e), {'', ''; 'et2d', 'm2d'}, 'fs', fs, 'xu', 'n', 'N', 4800);
    % AKp(irs_ext(:, d, e), {'', ''; 'et2d', 'm2d'}, 'fs', fs, 'xu', 'n', 'N', 4800, 'c', 'r');
end
