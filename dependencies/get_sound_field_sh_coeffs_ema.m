function S_breve = get_sound_field_sh_coeffs_ema( ...
    array_tfs, one_over_d_n, N, alpha_ema, sphharm_type, ischarm)
    
    if nargin < 5; sphharm_type = 'complex'; end
    if nargin < 6; ischarm = false; end
    
    %  ------------- PWD of the sound pressure along the equator --------------
    S_ring_m_surf = zeros(size(array_tfs, 1), 2*N+1);
    
    if strcmpi(sphharm_type, 'complex')
        
        for m = -N : N
            % discretized transformation integral
            S_ring_m_surf(:, m + N + 1) = ...
                sum(array_tfs .* exp(-1i .* m .* alpha_ema), 2) / size(alpha_ema, 2);
        end
    
    elseif strcmpi(sphharm_type, 'real')
        
        for m = -N : -1
            % discretized transformation integral
            S_ring_m_surf(:, m + N + 1) = ...
                sum(array_tfs .* sqrt(2) .* sin(abs(m) .* alpha_ema), 2) / size(alpha_ema, 2);
        end
        
        S_ring_m_surf(:, N + 1) = sum(array_tfs, 2) / size(alpha_ema, 2);
        
        for m = 1 : N
            % discretized transformation integral
            S_ring_m_surf(:, m + N + 1) = ...
                sum(array_tfs .* sqrt(2) .* cos(m .* alpha_ema), 2) / size(alpha_ema, 2);
        end

    else
        error('Unknown sphharm_type.');
    end
    
    % -------------------- compute A_m --------------------
    if isempty(one_over_d_n)
        S_ring_bar_m = S_ring_m_surf;
    else
        S_ring_bar_m = S_ring_m_surf .* one_over_d_n;
    end

    % ------- get the coefficients --------
    if ischarm
        S_breve = zeros(size(array_tfs, 1), 2*N + 1, 'like', S_ring_bar_m);
        ms = ch_stackOrder(N);
        for m = -N : N
            S_breve(:, ms == m) = S_ring_bar_m(:, m + N + 1);
        end
    else
        S_breve = zeros(size(array_tfs, 1), (N+1)^2, 'like', S_ring_bar_m);
        for n = 0 : N
            for m = -n : n
                S_breve(:, n^2+n+m+1) = S_ring_bar_m(:, m + N + 1) .* N_nm(n, m, pi/2);
            end
        end
    end

end

function ms = ch_stackOrder(order)
    ms = 0;
    for m = 1 : order
        ms = [ms, -m, m]; %#ok<*AGROW> 
    end
end
