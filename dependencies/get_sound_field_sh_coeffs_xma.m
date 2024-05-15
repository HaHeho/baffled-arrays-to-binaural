function S_breve = get_sound_field_sh_coeffs_xma( ...
    array_tfs, one_over_b_n, N, X_nm)
    
    % SH decomposition and radial filtering
    S_breve = zeros(size(array_tfs, 1), (N+1)^2, 'like', array_tfs);
    for n = 0 : N
        for m = -n : n
            if ~rem(n+m, 2) % non-zero coefficient
                acn = n^2+n+m+1;
                S_breve(:, acn) = sum(array_tfs .* X_nm(:, :, acn), 2) ...
                    .* one_over_b_n(:, n+1);

%                 % Use 1i^(-n) to be compatible to SMA binaural rendering 
%                 % equation i.e., compensate for using a different SMA
%                 % radial filter definition compared to the XMA publication
%                 S_breve(:, acn) = 1i^(-n) * sum(array_tfs .* X_nm(:, :, acn), 2) ...
%                     .* one_over_d_n(:, n+1);
            end
        end
    end
end
