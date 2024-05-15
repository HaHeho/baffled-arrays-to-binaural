function irs_diff = compute_spectral_difference_IRs(irs, irs_ref)

    % zero-pad data so the spectral error can be computed
    ir_len = max([size(irs, 1), size(irs_ref, 1)]);
    irs(end+1:ir_len, :, :) = 0;
    irs_ref(end+1:ir_len, :, :) = 0;

    [specs, is_even] = AKboth2singleSidedSpectrum(fft(irs));
    irs_diff = ifft(AKsingle2bothSidedSpectrum(compute_spectral_difference( ...
        specs, AKboth2singleSidedSpectrum(fft(irs_ref))), is_even));
end
