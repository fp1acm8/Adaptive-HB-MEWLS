function fig = plot_solution_surface(dataset, solverResults)
%PLOT_SOLUTION_SURFACE Compare solver predictions against the reference.
%   FIG = PLOT_SOLUTION_SURFACE(DATASET, SOLVERRESULTS) returns a figure
%   handle with subplots showing the ground-truth field and the predicted
%   surfaces produced by each solver.

arguments
    dataset (1, 1) struct
    solverResults (1, :) struct
end

nSolvers = numel(solverResults);
fig = figure('Name', 'Solution Surfaces', 'Color', 'w');

% Truth subplot
subplot(1, nSolvers + 1, 1);
scatter3(dataset.xy(:, 1), dataset.xy(:, 2), dataset.fTrue, 15, dataset.fTrue, 'filled');
title('Ground truth');
xlabel('x'); ylabel('y'); zlabel('f');
colormap turbo;
view(45, 30);

displayRange = [min(dataset.fTrue), max(dataset.fTrue)];

for i = 1:nSolvers
    subplot(1, nSolvers + 1, i + 1);
    pred = solverResults(i).prediction;
    scatter3(dataset.xy(:, 1), dataset.xy(:, 2), pred, 15, pred, 'filled');
    title(sprintf('%s prediction', solverResults(i).name));
    xlabel('x'); ylabel('y'); zlabel('f');
    colormap turbo;
    view(45, 30);
    caxis(displayRange);
end

end
