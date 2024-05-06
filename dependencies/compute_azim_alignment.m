function [min_shift, min_shift_ind, specs_diff, specs_diff_mean] = ...
    compute_azim_alignment(az, irs, irs_ref, fs, irs_frac_sm, diff_freq_rng)

    if nargin < 5; irs_frac_sm = 0; end % in fractional octaves
    if nargin < 6; diff_freq_rng = [500, 20e3]; end % in Hz
    
    if ~isvector(az)
        % Split off azimuth if elevation data is given as well
%         el = az(2, :);
        az = az(1, :);
%         az = -az(1, :); % transform azimuth into yaw, since it will be reversed in `plot_spec_azims()`
    end

    % zero-pad data so the spectral error can be computed
    ir_len = max([size(irs, 1), size(irs_ref, 1)]);
    irs(end+1:ir_len, :, :) = 0;
    irs_ref(end+1:ir_len, :, :) = 0;
    assert(all(size(irs) == size(irs_ref)));

%     % rearange data so there is no plotting artefacts
%     [az, az_idx] = sort(az);
%     irs = irs(:, az_idx, :);
%     irs_ref = irs_ref(:, az_idx, :);

%     % compute onsets in time domian (this can be used for comparison /
%     % validation of the aligment by spectral difference)
%     for e = 1 : size(irs, 3)
%         irs_onset(:, e) = AKonsetDetect(irs(:, :, e), 10, -6, 'rel', [2e3, fs]);
%         irs_ref_onset(:, e) = AKonsetDetect(irs_ref(:, :, e), 10, -6, 'rel', [2e3, fs]);
%     end
%     irs_onset = irs_onset - mean(irs_onset, 'all') + mean(irs_ref_onset, 'all');
%     AKf;
%     plot(irs_ref_onset, 'LineWidth', 2);
%     hold on;
%     plot(irs_onset, 'LineWidth', 2);
%     legend({'ref L', 'ref R', 'L', 'R'});
%     axis tight;
%     grid on;

    % calculate spectra
    specs = AKboth2singleSidedSpectrum(fft(irs));
    specs_ref = AKboth2singleSidedSpectrum(fft(irs_ref));
    if irs_frac_sm
        for e = 1 : size(specs, 3)
            specs(:, :, e) = AKfractOctSmooth(specs(:, :, e), ...
                'amp', fs, irs_frac_sm);
            specs_ref(:, :, e) = AKfractOctSmooth(specs_ref(:, :, e), ...
                'amp', fs, irs_frac_sm);
        end
    end
    freqs = linspace(0, fs/2, size(specs, 1)).';
    freqs_bin = get_freqs_inclusive(freqs, diff_freq_rng);

    % logarithmic weights from 1 to 0 in entire frequency range
    freqs_weight = [0; log10(freqs(end))-log10(freqs(2:end))];
    freqs_weight = freqs_weight / max(freqs_weight);

%     % logarithmic weights from 1 to 0 in specified frequency range
%     freqs_bins = [find(freqs_bin, 1, 'first'), find(freqs_bin, 1, 'last')];
%     freqs_weight = [zeros(freqs_bins(1), 1); ...
%         log10(freqs(freqs_bins(2)))-log10(freqs(freqs_bins(1):freqs_bins(2))); ...
%         zeros(length(freqs) - freqs_bins(2) - 1, 1)];
%     freqs_weight = freqs_weight / max(freqs_weight);

    % cut frequency bins outside of specified range
    % (for performance reasons)
    specs = specs(freqs_bin, :, :);
    specs_ref = specs_ref(freqs_bin, :, :);
    freqs = freqs(freqs_bin);
    freqs_weight = freqs_weight(freqs_bin);

    specs_diff_bin = zeros(length(freqs), size(irs, 2), size(irs, 3));
    for a = 1 : length(az) % each direction
        if a > 1
            % rotate dataset
            specs = circshift(specs, 1, 2);
        end
        % compute differences between datasets for all incidence directions
        spec_diff_bin = compute_spectral_difference(specs, specs_ref);

%         % this would be easier to understand but is less performant
%         specs_rot = circshift(specs, round(az(a)), 2);
%         spec_diff_bin = compute_spectral_difference(specs_rot, specs_ref);

        % compute absolute differences
        % (equivalent to: spec_diff = db2mag(abs(db(spec_diff)));)
        spec_diff_bin(spec_diff_bin < 1) = 1 ./ spec_diff_bin(spec_diff_bin < 1);

        % mean differences over all incidence directions
        specs_diff_bin(:, a, :) = mean(spec_diff_bin, 2);
    end
    
    % mean weighted differences over all frequency bins
    specs_diff_mean = squeeze(mean(specs_diff_bin .* freqs_weight, 1, 'omitnan'));

    % determine difference minima for individual ears
    min_idx_ind = zeros(size(specs_diff_mean, 2), 1);
    for e = 1 : size(specs_diff_mean, 2)
        [~, min_idx_ind(e)] = min(specs_diff_mean(:, e));
    end

    % determine difference minima for both ears combined
    [~, min_idx] = min(sum(specs_diff_mean, 2));
    
    % translate shift index into angle
    min_shift_ind = az(min_idx_ind);
    min_shift = az(min_idx);

    % refill frequency bins outside of specified range (for plotting)
    specs_diff = nan(length(freqs_bin), size(irs, 2), size(irs, 3));
    specs_diff(freqs_bin, :, :) = specs_diff_bin;
end
