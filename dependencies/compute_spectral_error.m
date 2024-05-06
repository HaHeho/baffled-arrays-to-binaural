function [spec_err, spec_log_err, freqs_log] = compute_spectral_error(specs_diff, dims, freqs)
    
    if nargin < 2; dims = ndims(specs_diff); end
    if nargin < 3; freqs = []; end
    
    if ~isreal(specs_diff); specs_diff = abs(specs_diff); end

%     % absolute mean with logarithmic scaling
%     spec_err = db2mag(squeeze(mean(abs(db(specs_diff)), dims)));
%     % absolute RMS with logarithmic scaling
%     spec_err = db2mag(squeeze(rms(db(specs_diff), dims)));
%     % absolute mean without logarithmic scaling
%     spec_err = squeeze(mean(db2mag(abs(db(specs_diff))), dims));
    % absolute RMS without logarithmic scaling
    spec_err = squeeze(rms(db2mag(abs(db(specs_diff))), dims));
%     % equivalent
%     spec_err = compute_spectral_average(db2mag(abs(db(specs_diff))), dims);

    if nargout > 1
        assert(~isempty(freqs), 'Given frequency vector is required for interpolation.');

%         % generate logarithmic frequency vector with identical length
%         freqs_log = [0, logspace(log10(freqs(2)), log10(freqs(end)), length(freqs) - 1)];
        % generate logarithmic frequency vector with different length
        % (helpful to prevent accidentally plotting with wrong frequencies)
        freqs_log = logspace(log10(freqs(2)), log10(freqs(end)), length(freqs) + 1);
    
        % interpolate to target frequencies
        % (to calculate perceptually-spaced mean over all frequencies later)
        spec_log_err = interp1(freqs, spec_err, freqs_log, 'linear', 'extrap');
    end
end
