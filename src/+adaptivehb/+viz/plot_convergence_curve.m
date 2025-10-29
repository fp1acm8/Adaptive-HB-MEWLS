function fig = plot_convergence_curve(solverResults)
%PLOT_CONVERGENCE_CURVE Plot RMSE decay for each solver.

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
