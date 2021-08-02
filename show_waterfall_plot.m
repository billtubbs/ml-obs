function p = show_waterfall_plot(t, Z, z_minmax, ax_labels, view_angles, filename)
% Make a Waterfall plot in current axis showing the
% time series in Z (y and z axes) over time t (x-axis).
%
% Arguments:
%   t : Time vector (1 column)
%   Z : Time series data in columns
%   ax_labels : cell array of labels for x, y, z axes
%   view : (optional) array [az el] containing azimuth and
%       elevation angles of plot view (degrees). If not
%       [0 84] is used.
%   filename : (optional) save figure to file.
%

    if nargin == 3
        view_angles = [0 84];
    end
    n_series = size(Z, 2);
    [X, Y] = meshgrid(t', 1:n_series);
    p = waterfall(X, Y, Z', Y);
    ylim([1 n_series])
    xlabel(ax_labels{1}, 'Interpreter', 'Latex')
    ylabel(ax_labels{2}, 'Interpreter', 'Latex')
    zlabel(ax_labels{3}, 'Interpreter', 'Latex')
    zlim(z_minmax);
    set(gca, 'TickLabelInterpreter', 'latex')
    set(p, 'FaceColor', 'flat');
    set(p, 'EdgeColor', 'k');
    set(p, 'FaceVertexCData', sqrt(parula(n_series)));
    set(p, 'FaceAlpha', 0.75)
    %position = get(gcf, 'Position');
    %set(gcf, 'Position', [position(1:2) 560 420]);

    % Default viewing angles are view(-37.5, 30)
    % From front at high angle from above: view(0, 85)
    view(view_angles)

    if nargin > 5
        saveas(gcf, filename)
    end
end