function S_breve = get_sound_field_sh_coeffs_sma( ...
    array_tfs, one_over_d_n, N, alpha_beta, grid_weights, sh_type)
    
    if nargin < 6; sh_type = 'complex'; end
    if nargin < 4
        % The inverted basis functions have been pre-computed and are given
        % instead of the `N` parameter
        Y_nm_inv = N;
        N = sqrt(size(Y_nm_inv, 1)) - 1;
    else
        Y_nm_inv = get_basis_inv_sh_coeffs_sma(N, alpha_beta, grid_weights, sh_type);
    end
    
    % SH decomposition
    S_breve = zeros(size(array_tfs, 1), size(Y_nm_inv, 1), 'like', array_tfs);
    for n = 0 : N
        for m = -n : n
            acn = n^2+n+m+1;
            % 4*pi is not included here since the quadrature weights are
            % normalized to 4*pi
            S_breve(:, acn) = sum(array_tfs .* Y_nm_inv(acn, :), 2);

            % Apply radial filters
            if ~isempty(one_over_d_n)
                S_breve(:, acn) = S_breve(:, acn) .* one_over_d_n(:, n+1);
            end
        end
    end
end
