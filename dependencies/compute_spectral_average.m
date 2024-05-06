function spec_avg = compute_spectral_average(specs, dims)

    if nargin < 2; dims = ndims(specs); end
    
    if ~isreal(specs); specs = abs(specs); end
    
%     % mean with logarithmic scaling
%     specs_db = mag2db(abs(specs));
%     % calculate average
%     spec_avg_db = squeeze(mean(specs_db, dims));
%     % revert logarithmic scaling
%     spec_avg = db2mag(spec_avg_db);
    
%     % mean without logarithmic scaling
%     spec_avg = squeeze(mean(specs, dims));

    % RMS without logarithmic scaling
    spec_avg = squeeze(rms(specs, dims));
end
