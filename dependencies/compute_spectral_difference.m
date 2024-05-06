function spec_diff = compute_spectral_difference(spec, spec_ref)
    
    if ~isreal(spec); spec = abs(spec); end
    if ~isreal(spec_ref); spec_ref = abs(spec_ref); end

    spec_diff = spec ./ spec_ref;
end
