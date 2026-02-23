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

arguments
    dataset (1, 1) struct
    solverResults (1, :) struct
end

nSolvers = numel(solverResults);
fig = figure('Name', 'Error Distributions', 'Color', 'w');

for i = 1:nSolvers
    subplot(1, nSolvers, i);
    residual = dataset.fTrue - solverResults(i).prediction;
    histogram(residual, 'Normalization', 'pdf');
    title(sprintf('%s residuals', solverResults(i).name));
    xlabel('Error');
    ylabel('PDF');
    grid on;
end

end
