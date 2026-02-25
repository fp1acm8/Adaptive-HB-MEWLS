function fig = plot_error_distribution(dataset, solverResults)
%PLOT_ERROR_DISTRIBUTION Visualise error histograms for each solver.
%   FIG = PLOT_ERROR_DISTRIBUTION(DATASET, SOLVERRESULTS) creates one
%   subplot per solver showing the histogram of errors (f_true - QI),
%   normalised as a probability density function (PDF).
%
%   Inputs:
%       DATASET       - struct with f_true (N-by-1).
%       SOLVERRESULTS - 1-by-K struct array with fields QI and name.
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

    % Compute the point-wise error: e = f_true - QI.
    % Positive = under-prediction, negative = over-prediction.
    error = dataset.f_true - solverResults(i).QI;

    % Draw a histogram normalised as a probability density function (PDF).
    % 'pdf' normalisation ensures the area under the histogram equals 1,
    % making histograms comparable across solvers with different bin counts.
    histogram(error, 'Normalization', 'pdf');

    title(sprintf('%s errors', solverResults(i).name));
    xlabel('Error');    % error value on the x-axis
    ylabel('PDF');      % probability density on the y-axis
    grid on;            % add grid lines for easier reading
end

end
