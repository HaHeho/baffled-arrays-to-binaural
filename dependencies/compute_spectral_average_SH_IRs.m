function irs_avg = compute_spectral_average_SH_IRs(irs_sh, dims)

    if nargin < 2; dims = 2; end

    [specs_sh, is_even] = AKboth2singleSidedSpectrum(fft(irs_sh));
    irs_avg = ifft(AKsingle2bothSidedSpectrum( ...
        compute_spectral_average_SH(specs_sh, dims), is_even));
end
