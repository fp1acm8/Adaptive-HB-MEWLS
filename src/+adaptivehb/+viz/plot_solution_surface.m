function fig = plot_solution_surface(dataset, solverResults)
%PLOT_SOLUTION_SURFACE Compare solver predictions against ground truth.
%   FIG = PLOT_SOLUTION_SURFACE(DATASET, SOLVERRESULTS) creates a figure
%   with (nSolvers + 1) side-by-side 3-D scatter subplots: the first shows
%   the ground-truth field (dataset.fTrue), and each subsequent subplot
%   shows one solver's prediction.  All subplots share the same colour axis
%   (clamped to the ground-truth range) and use the "turbo" colourmap with
%   a 45/30-degree view angle.
%
%   Inputs:
%       DATASET       - struct with xy (N-by-2) and fTrue (N-by-1).
%       SOLVERRESULTS - 1-by-K struct array; each element must have fields
%                       prediction (N-by-1) and name (string).
%
%   Output:
%       FIG - figure handle (caller is responsible for saving/closing).
%
%   See also adaptivehb.viz.plot_error_distribution,
%            adaptivehb.viz.plot_convergence_curve.

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
