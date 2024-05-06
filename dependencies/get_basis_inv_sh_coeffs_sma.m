function Ynm = get_basis_inv_sh_coeffs_sma(N, alpha_beta, grid_weights, sh_type)
    
    if nargin < 4; sh_type = 'complex'; end

    % SH basis functions
    Ynm = zeros(size(alpha_beta, 2), (N+1)^2);
    for n = 0 : N
        for m = -n : n
            Ynm(:, n^2+n+m+1) = sphharm(n, m, ...
                alpha_beta(2, :), alpha_beta(1, :), sh_type);
        end
    end
    if isempty(grid_weights)
        % Perform pseudo-inverse (based on entire set of SH coefficients!)
        Ynm = pinv(Ynm);
    else
        % Apply quadrature weights
        Ynm = 4*pi / sum(grid_weights) * grid_weights .* Ynm';
    end
end
