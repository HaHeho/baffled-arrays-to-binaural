function spec_perc = compute_spectral_percentiles(specs, dims, percentiles)
    
    if nargin < 2; dims = ndims(specs); end
    if nargin < 3; percentiles = [10, 90]; end
    
    if ~isreal(specs); specs = abs(specs); end
    
%     % mean without logarithmic scaling
%     spec_mean = mean(specs, dims);
%     % percentiles without logarithmic scaling
%     spec_perc = prctile(specs, percentiles, dims);
%     % combine
%     spec_perc = cat(dims, spec_mean, spec_perc);

    % percentiles without logarithmic scaling
    spec_perc = prctile(specs, [50; percentiles(:)], dims);
end
