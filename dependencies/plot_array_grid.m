function plot_array_grid(azi_col_rad, r, do_plot_labels)

    OUTER_SCALE = 1.4;
    DOT_SIZE = 30;
    SPHERE_RES = 100;
    
    if nargin < 3; do_plot_labels = true; end
    if nargin < 2; r = 1; end % in m
    if length(r) > 1 && size(azi_col_rad, 1) ~= length(r)
        error('Mismatch in parameter size.');
    end
    
    % get cartesian microphone coordinates
    [m_x, m_y, m_z] = sph2cart(azi_col_rad(:, 1), ...
        pi/2 - azi_col_rad(:, 2), r * 1.01); % colatitude to elevation
    
    % get cartesian body coordinates
    b_r          = mean(r); 
    b_az         = linspace(0, 2*pi, 2*SPHERE_RES);
    b_el         = pi/2 - linspace(0, pi, SPHERE_RES + 1);
    [b_az, b_el] = meshgrid(b_az, b_el);
    [b_x, b_y, b_z] = sph2cart(b_az, b_el, b_r * ones(size(b_az)));
    
    % plot body
    hold on;
    surf(b_x, b_y, b_z);

    % plot coordinate axes
    line([-1 1] * b_r * OUTER_SCALE, [.0 .0], [0 0], 'Marker', '.', ...
        'LineStyle', '-', 'Color', [.5 .5 .5], 'LineWidth', 1);
    line([ 0 0], [-1  1] * b_r * OUTER_SCALE, [0 0], 'Marker', '.', ...
        'LineStyle', '-', 'Color', [.5 .5 .5], 'LineWidth', 1);
    line([ 0 0], [ 0  0], [-1 1] * b_r * OUTER_SCALE, 'Marker', '.', ...
        'LineStyle', '-', 'Color', [.5 .5 .5], 'LineWidth', 1);

    % plot mics
    scatter3(m_x, m_y, m_z, DOT_SIZE, 'ko', 'MarkerFaceColor', 'k');

    colormap gray;
    caxis([-5*r, 1.2*r]);

    shading interp;
    axis equal;
    box on;
    grid on;

    axis([-1, 1, -1, 1, -1, 1] * b_r * OUTER_SCALE);
    view(30, 15);

    ax = gca;
    if do_plot_labels
        drawnow; % this is required before requesting the ticks
        yticks(xticks(ax));
        zticks(xticks(ax));

        xlabel('$x$ in m', 'Interpreter', 'Latex');
        ylabel('$y$ in m', 'Interpreter', 'Latex');
        zlabel('$z$ in m', 'Interpreter', 'Latex');
    else
        xlabel('');
        ylabel('');
        zlabel('');
        xticklabels([]);
        yticklabels([]);
        zticklabels([]);

        set(ax, 'TickLength', [0, 0]);
    end
    set(ax, 'Projection', 'perspective');

    drawnow;
end
