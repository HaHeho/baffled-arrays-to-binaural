function [grid, grid_weigths] = get_array_grid(type, order, is_rad, is_col)

% requires Spherical-Harmonic-Transform toolbox
% $ git clone https://github.com/polarch/Spherical-Harmonic-Transform.git
% 
% requires `getLebedevSphere()` from
% https://se.mathworks.com/matlabcentral/fileexchange/27097-getlebedevsphere
%
% requires `getEigenmikeNodes()` from
% https://www.mathworks.com/matlabcentral/fileexchange/73721-geteigenmikenodes

    if nargin < 3; is_rad = true; end
    if nargin < 4; is_col = true; end

    if strcmpi(type, 'Fliege')
        % from Spherical-Harmonic-Transform toolbox
        [~, grid, ~] = getFliegeNodes(order + 1); % in rad, elevation
        grid_weigths = ones(size(grid, 1), 1);

    elseif strcmpi(type, 'Lebedev')
        % from https://se.mathworks.com/matlabcentral/fileexchange/27097-getlebedevsphere
        grid_cart = getLebedevSphere(lebedevOrder2Degree(2*order + 1));
        [grid(:, 1), grid(:, 2), ~] = cart2sph( ...
            grid_cart.x, grid_cart.y, grid_cart.z); % in rad, elevation
        grid_weigths = grid_cart.w;

    elseif strcmpi(type, 'tDesign')
%         % from Spherical-Harmonic-Transform toolbox
%         % (not identical to pyfar)
%         [~, grid] = getTdesign(order * 2); % in rad

        % identical to pyfar as used for the SMA measurements
        % https://github.com/pyfar/pyfar/blob/41478c9192333d2febb33ab548b43494dbd2a887/pyfar/samplings/samplings.py#L1032
        target_degree = order * 2;
        switch target_degree
            case 3
                n_ch = 8;
            case 5
                n_ch = 18;
            case 7
                n_ch = 32;
            case 9
                n_ch = 50;
            case 11
                n_ch = 72;
            case 13
                n_ch = 98;
            case 15
                n_ch = 128;
            otherwise
                n_ch = ceil((target_degree + 1)^2 / 2) + 1;
        end
        grid_cart = readmatrix(sprintf( ...
            'https://web.maths.unsw.edu.au/~rsw/Sphere/Points/SF/SF29-Nov-2012/sf%03d.%05d', ...
            target_degree, n_ch), 'FileType', 'text');
        [grid(:, 1), grid(:, 2)] = cart2sph( ...
            grid_cart(:, 1), grid_cart(:, 2), grid_cart(:, 3)); % in rad, elevation
        grid_weigths = ones(size(grid, 1), 1);

    elseif any(strcmpi(type, {'Eigenmike', 'Eigenmike32'}))
        % from https://www.mathworks.com/matlabcentral/fileexchange/73721-geteigenmikenodes
        assert(order == 4, 'The Eigenmike grid is only defined for order 4.')
        grid = getEigenmikeNodes('rad', false); % in rad, elevation
        grid = grid(:, 1:2); % remove radius

    elseif strcmpi(type, 'Equatorial')
        n_ch = 2*order + 1;
        grid = linspace(0, 2*pi, n_ch + 1).'; % in rad
        grid = grid(1:end-1);
        grid = [grid, zeros(n_ch, 1)]; % add elevation

    elseif any(strcmpi(type, {'SSR', 'SSR_vertical', 'vertical_SSR'}))
        n_ch = 360;
        % use source position instead of incidence direction
        grid = -linspace(0, 2*pi, n_ch + 1).'; % in rad
        grid = grid(1:end-1);
        grid = [grid, zeros(n_ch, 1)]; % add elevation
        if any(strcmpi(type, {'SSR_vertical', 'vertical_SSR'}))
            grid = flip(grid, 2); % switch azimuth and elevation
        end

    else
        error('Grid type "%s" not implemented yet.', type);
    end

    % limit coordinate range to convention required for SH processing
    % (by converting to cartesian coordinates and back)
    % azimuth:     -pi   .. +pi     -180 deg .. +180 deg
    % elevation:   -pi/2 .. +pi/2    -90 deg ..  +90 deg
    [grid_x, grid_y, grid_z] = sph2cart(grid(:, 1), grid(:, 2), 1);
    [grid(:, 1), grid(:, 2), ~] = cart2sph(grid_x, grid_y, grid_z);

    if is_col
        % elevation to colatitude
        grid(:, 2) = pi/2 - grid(:, 2);
    end

    if ~is_rad
        % radians to degrees
        grid = rad2deg(grid);
    end

    if ~exist('grid_weigths', 'var')
        grid_weigths = NaN;
    else
        % normalize to 4*pi
        grid_weigths = 4*pi / sum(grid_weigths) * grid_weigths;
    end
end

function degree = lebedevOrder2Degree(order)
    % this combination is found in the header of getLebedevSphere()
    degrees = [6,14,26,38,50,74,86,110,146,170,194,230,266,302,350,434,590,...
        770,974,1202,1454,1730,2030,2354,2702,3074,3470,3890,4334,4802,5294,5810];
    orders = [3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,35,41,47,53,59,65,...
        71,77,83,89,95,101,107,113,119,125,131];
    
    i = find(orders==order, 1);
    if isempty(i)
        error('Invalid order %d.', order);
    end
    degree = degrees(i);
end
