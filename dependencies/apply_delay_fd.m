function tfs = apply_delay_fd(tfs, delay_samples)
%     % instead of
%     irs = ifft(AKsingle2bothSidedSpectrum(tfs));
%     irs = circshift(irs, round(delay_samples));
%     tfs = AKboth2singleSidedSpectrum(fft(irs));

    omega = linspace(0, 0.5, size(tfs, 1)).';
    tfs = tfs .* exp(-1j * 2 * pi * omega * delay_samples);
end
