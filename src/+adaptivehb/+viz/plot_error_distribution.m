function fig = plot_error_distribution(dataset, solverResults)
%PLOT_ERROR_DISTRIBUTION Visualise residual histograms for each solver.
%   FIG = PLOT_ERROR_DISTRIBUTION(DATASET, SOLVERRESULTS) creates one
%   subplot per solver showing the histogram of residuals (fTrue - pred),
%   normalised as a probability density function (PDF).
%
%   Inputs:
%       DATASET       - struct with fTrue (N-by-1).
%       SOLVERRESULTS - 1-by-K struct array with fields prediction and name.
%
%   Output:
%       FIG - figure handle.
%
%   See also adaptivehb.viz.plot_solution_surface,
%            adaptivehb.viz.plot_convergence_curve.

% --- Input validation -------------------------------------------------
arguments
    dataset (1, 1) struct
    solverResults (1, :) struct
end

% Count the number of solvers (one subplot each).
nSolvers = numel(solverResults);

% Create a new figure with white background.
fig = figure('Name', 'Error Distributions', 'Color', 'w');

for i = 1:nSolvers
    % Place this solver's histogram in the i-th subplot (single row).
    subplot(1, nSolvers, i);

    % Compute the point-wise residual: e = fTrue - prediction.
    % Positive = under-prediction, negative = over-prediction.
    residual = dataset.fTrue - solverResults(i).prediction;

    % Draw a histogram normalised as a probability density function (PDF).
    % 'pdf' normalisation ensures the area under the histogram equals 1,
    % making histograms comparable across solvers with different bin counts.
    histogram(residual, 'Normalization', 'pdf');

    title(sprintf('%s residuals', solverResults(i).name));
    xlabel('Error');    % residual value on the x-axis
    ylabel('PDF');      % probability density on the y-axis
    grid on;            % add grid lines for easier reading
end

end
