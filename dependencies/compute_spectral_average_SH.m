function spec_avg = compute_spectral_average_SH(specs_sh, dims, norm_dim_size)

    if nargin < 2; dims = 2; end
    if nargin < 3; norm_dim_size = false; end
        
    if ~isreal(specs_sh); specs_sh = abs(specs_sh); end
    
    spec_avg = squeeze(rms(specs_sh, dims)) * 4*pi / sqrt(2);

    if norm_dim_size
        spec_avg = spec_avg * sqrt(size(specs_sh, dims));
    end
end
