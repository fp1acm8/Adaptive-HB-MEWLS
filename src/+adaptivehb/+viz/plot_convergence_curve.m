function fig = plot_convergence_curve(solverResults)
%PLOT_CONVERGENCE_CURVE Plot RMSE vs polynomial degree for each solver.
%   FIG = PLOT_CONVERGENCE_CURVE(SOLVERRESULTS) overlays the RMSE
%   convergence curves of all solvers on a single axes.  Each solver's
%   convergence table (fields: degree, rmse) provides the data points,
%   plotted as lines with circle markers.
%
%   Input:
%       SOLVERRESULTS - 1-by-K struct array; each element must have fields
%                       convergence (table with degree and rmse columns)
%                       and name (string used in the legend).
%
%   Output:
%       FIG - figure handle.
%
%   See also adaptivehb.viz.plot_solution_surface,
%            adaptivehb.viz.plot_error_distribution.

% --- Input validation -------------------------------------------------
arguments
    solverResults (1, :) struct
end

% Create a new figure with white background.
fig = figure('Name', 'Convergence Curves', 'Color', 'w');

% Enable hold so that multiple plot() calls overlay on the same axes.
hold on;

for i = 1:numel(solverResults)
    % Extract the convergence table for this solver.
    % Table columns: degree (1..maxDegree), rmse, maxAbsError, mae.
    conv = solverResults(i).convergence;

    % Plot RMSE vs polynomial degree as a line with circle markers.
    % 'DisplayName' feeds the legend with the solver's name.
    plot(conv.degree, conv.rmse, '-o', 'DisplayName', solverResults(i).name);
end

hold off;

% Label the axes clearly.
xlabel('Polynomial degree');  % x-axis: total degree d of the polynomial
ylabel('RMSE');               % y-axis: root mean squared error

% Show a legend identifying each solver's curve.
legend('Location', 'northeast');

% Add grid lines for easier visual comparison of RMSE values.
grid on;
end
