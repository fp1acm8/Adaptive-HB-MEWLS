function fig = plot_error_distribution(dataset, solverResults)
%PLOT_ERROR_DISTRIBUTION Visualise residual histograms for each solver.

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
