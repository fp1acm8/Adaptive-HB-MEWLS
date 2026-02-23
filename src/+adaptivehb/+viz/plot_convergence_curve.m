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

arguments
    solverResults (1, :) struct
end

fig = figure('Name', 'Convergence Curves', 'Color', 'w');
hold on;
for i = 1:numel(solverResults)
    conv = solverResults(i).convergence;
    plot(conv.degree, conv.rmse, '-o', 'DisplayName', solverResults(i).name);
end
hold off;
xlabel('Polynomial degree');
ylabel('RMSE');
legend('Location', 'northeast');
grid on;
end
