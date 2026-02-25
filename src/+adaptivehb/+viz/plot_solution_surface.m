function fig = plot_solution_surface(dataset, solverResults)
%PLOT_SOLUTION_SURFACE Compare solver predictions against ground truth.
%   FIG = PLOT_SOLUTION_SURFACE(DATASET, SOLVERRESULTS) creates a figure
%   with (nSolvers + 1) side-by-side 3-D scatter subplots: the first shows
%   the ground-truth field (dataset.f_true), and each subsequent subplot
%   shows one solver's QI approximation.  All subplots share the same colour
%   axis (clamped to the ground-truth range) and use the "turbo" colourmap
%   with a 45/30-degree view angle.
%
%   Inputs:
%       DATASET       - struct with data (N-by-2) and f_true (N-by-1).
%       SOLVERRESULTS - 1-by-K struct array; each element must have fields
%                       QI (N-by-1) and name (string).
%
%   Output:
%       FIG - figure handle (caller is responsible for saving/closing).
%
%   See also adaptivehb.viz.plot_error_distribution,
%            adaptivehb.viz.plot_convergence_curve.

% --- Input validation -------------------------------------------------
arguments
    dataset (1, 1) struct
    solverResults (1, :) struct
end

% Count how many solvers we need to display.
nSolvers = numel(solverResults);

% Create a new figure with white background.
fig = figure('Name', 'Solution Surfaces', 'Color', 'w');

% --- First subplot: ground truth --------------------------------------
% This provides the visual reference that all solver panels are compared to.
subplot(1, nSolvers + 1, 1);
% scatter3 draws coloured 3-D points: x, y, z=f_true, markerSize=15, colour=f_true.
scatter3(dataset.data(:, 1), dataset.data(:, 2), dataset.f_true, 15, dataset.f_true, 'filled');
title('Ground truth');
xlabel('x'); ylabel('y'); zlabel('f');
colormap turbo;       % perceptually uniform colormap (red-yellow-green-blue)
view(45, 30);         % azimuth=45°, elevation=30° for a clear 3-D perspective

% Compute the colour axis range from the ground truth so that all subplots
% share the same scale, making visual comparison meaningful.
displayRange = [min(dataset.f_true), max(dataset.f_true)];

% --- Solver subplots --------------------------------------------------
for i = 1:nSolvers
    % Place each solver in the next subplot position.
    subplot(1, nSolvers + 1, i + 1);

    % Extract this solver's QI approximation values.
    QI = solverResults(i).QI;

    % Plot QI surface with the same scatter3 style as ground truth.
    scatter3(dataset.data(:, 1), dataset.data(:, 2), QI, 15, QI, 'filled');
    title(sprintf('%s (QI)', solverResults(i).name));
    xlabel('x'); ylabel('y'); zlabel('f');
    colormap turbo;
    view(45, 30);

    % Lock the colour axis to the ground-truth range so colour differences
    % between panels are due to actual value differences, not rescaling.
    caxis(displayRange);
end

end
