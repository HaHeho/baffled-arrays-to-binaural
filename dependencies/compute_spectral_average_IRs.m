function irs_avg = compute_spectral_average_IRs(irs, ...
    dims) %, before_frac_sm, fs)

    if nargin < 2; dims = 2; end
%     if nargin < 3; before_frac_sm = 0; end % in fractional octaves

    [specs, is_even] = AKboth2singleSidedSpectrum(fft(irs));
%     if before_frac_sm
%         for e = 1:size(specs, 3)
%             specs(:, :, e) = AKfractOctSmooth(specs(:, :, e), ...
%                 'amp', fs, before_frac_sm);
%         end; clear e;
%     end
    irs_avg = ifft(AKsingle2bothSidedSpectrum( ...
        compute_spectral_average(specs, dims), is_even));
end
