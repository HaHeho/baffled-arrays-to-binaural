function output_block = render_sh_coefficients_binaurally_blockwise( ...
    hrtfs_nm, S_breve, head_dir_azel_rad, N, sphharm_type)

    if nargin < 5; sphharm_type = 'complex'; end
    
    n_dirs = size(head_dir_azel_rad, 2);
    is_ch = size(S_breve, 2) == 2*N + 1;
    
    ear_specs = zeros(size(hrtfs_nm, 1), n_dirs, size(hrtfs_nm, 3), 'like', hrtfs_nm);
    if all(head_dir_azel_rad(2, :) == 0)
        for n = 0 : N
            for m = -n : n
                if is_ch
                    acn = [2*n, 2*n+1];
                    if n == 0
                        acn(1) = acn(2);
                    elseif m == n
                        acn = flip(acn);
                    elseif m ~= -n
                        continue;
                    end
                else
                    acn = [n^2+n+m+1, n^2+n-m+1];
                end
                if strcmpi(sphharm_type, 'complex')
                    ear_specs = ear_specs + S_breve(:, acn(1)) ...
                        .* hrtfs_nm(:, acn(2), :) .* exp(-1i .* m .* head_dir_azel_rad(1, :));
                elseif strcmpi(sphharm_type, 'real')
                    ear_specs = ear_specs + S_breve(:, acn(1)) ...
                         .* (hrtfs_nm(:, acn(1), :) .* cos(m * head_dir_azel_rad(1, :)) ...
                           - hrtfs_nm(:, acn(2), :) .* sin(m * head_dir_azel_rad(1, :)));
                else
                    error('SH type "%s" not implemented yet.', sphharm_type);
                end
            end
        end
    else
        assert(~is_ch, 'CHs for 3-DoF rotations not implemented yet.');

        % Rotate sound field in opposite direction
        % (this should be slightly more efficient, since the sound field is
        % only 2-dimensional, compared to the 3-dimensional HRTF)
        ear_specs = zeros([size(hrtfs_nm, 1), n_dirs, size(hrtfs_nm, 3)], 'like', S_breve);
        for d = 1 : n_dirs
            % Get rotation matrix
            M_zyx = euler2rotationMatrix(-head_dir_azel_rad(1, d), -head_dir_azel_rad(2, d), 0, 'zyx');
            % Compute rotation matrix in the SH domain
            S_breve_rot = S_breve * getSHrotMtx(M_zyx, N, sphharm_type);
            % Combine sound field and HRTFs
            ear_specs(:, d, :) = sum(S_breve_rot .* hrtfs_nm, 2);
        end
    end
    
    % Complement spectrum and iDFT
    output_block = ifft(AKsingle2bothSidedSpectrum(ear_specs));
end
