function run_comparison(configPath)
%RUN_COMPARISON Top-level orchestrator for solver comparisons.
%   RUN_COMPARISON(CONFIGPATH) loads a JSON configuration, prepares the
%   requested dataset (including optional perturbations) and invokes the
%   solvers specified in the configuration using identical inputs.
%   Numerical outputs are aggregated into a table with metric deltas and
%   stored alongside a suite of diagnostic plots inside reports/<timestamp>.
%
%   When CONFIGPATH is omitted the default configuration located at
%   config/comparison_config.json is used.
%
%   Example:
%       run_comparison();
%       run_comparison('config/my_experiment.json');
%
%   The function assumes the repository root is the current working folder.
%
%   See also adaptivehb.io.load_dataset, adaptivehb.io.apply_noise,
%            adaptivehb.solvers.leastSquaresSolver,
%            adaptivehb.solvers.mewlsSolver.

% Ensure the package folder is visible on the MATLAB path. We avoid
% permanently modifying the environment by restoring the path on cleanup.
rootDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(rootDir);
srcDir = fullfile(repoRoot, 'src');
oldPath = path;
cleanupPath = onCleanup(@() path(oldPath)); %#ok<NASGU>
addpath(srcDir);

if nargin == 0 || isempty(configPath)
    configPath = fullfile('config', 'comparison_config.json');
end

if ~isfile(configPath)
    error('adaptivehb:comparison:MissingConfig', ...
        'Configuration file not found: %s', configPath);
end

configText = fileread(configPath);
config = jsondecode(configText);

% Load dataset and apply optional noise.
dataset = adaptivehb.io.load_dataset(config.dataset.path, config.dataset);
dataset = adaptivehb.io.apply_noise(dataset, config.noise);

% Dispatch solvers.
solverSpecs = config.solvers;
solverResults = repmat(struct(), 1, numel(solverSpecs));

for i = 1:numel(solverSpecs)
    solverFcn = str2func(solverSpecs(i).function);
    if isfield(solverSpecs(i), 'parameters')
        solverParams = solverSpecs(i).parameters;
    else
        solverParams = struct();
    end

    try
        solverResults(i) = solverFcn(dataset, solverParams);
        solverResults(i).name = solverSpecs(i).name;
    catch ME
        warning('adaptivehb:comparison:SolverFailed', ...
            'Solver "%s" failed: %s', solverSpecs(i).name, ME.message);
        solverResults(i).name = solverSpecs(i).name;
        solverResults(i).prediction = NaN(size(dataset.fObserved));
        solverResults(i).coefficients = [];
        solverResults(i).metrics = struct('rmse', NaN, 'maxAbsError', NaN, 'mae', NaN);
        solverResults(i).convergence = table();
    end
end

% Aggregate metrics.
metricNames = string(config.metrics(:));
metricsMatrix = zeros(numel(solverResults), numel(metricNames));
for i = 1:numel(solverResults)
    for j = 1:numel(metricNames)
        key = metricNames(j);
        metricsMatrix(i, j) = solverResults(i).metrics.(key);
    end
end

metricsTable = array2table(metricsMatrix, 'VariableNames', cellstr(metricNames));
solverNames = string({solverSpecs.name});
metricsTable.solver = solverNames';
metricsTable = movevars(metricsTable, 'solver', 'Before', 1);

% Compute deltas relative to the first solver.
referenceRow = metricsMatrix(1, :);
deltaMatrix = metricsMatrix - referenceRow;
deltaTable = array2table(deltaMatrix, 'VariableNames', strcat('delta_', cellstr(metricNames)));
deltaTable.solver = solverNames';
deltaTable = movevars(deltaTable, 'solver', 'Before', 1);

summaryTable = [metricsTable, deltaTable(:, 2:end)];

% Prepare report directory.
timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
reportDir = fullfile('reports', ['comparison_' timestamp]);
if ~exist(reportDir, 'dir')
    mkdir(reportDir);
end

% Persist tables.
writetable(summaryTable, fullfile(reportDir, 'metrics.csv'));
save(fullfile(reportDir, 'metrics_table.mat'), 'summaryTable');

metricsStruct = table2struct(summaryTable);
jsonText = jsonencode(metricsStruct, 'PrettyPrint', true);
fid = fopen(fullfile(reportDir, 'metrics.json'), 'w');
fprintf(fid, '%s', jsonText);
fclose(fid);

% Generate plots.
figSurface = adaptivehb.viz.plot_solution_surface(dataset, solverResults);
saveas(figSurface, fullfile(reportDir, 'solution_surfaces.png'));
close(figSurface);

figError = adaptivehb.viz.plot_error_distribution(dataset, solverResults);
saveas(figError, fullfile(reportDir, 'error_distributions.png'));
close(figError);

figConv = adaptivehb.viz.plot_convergence_curve(solverResults);
saveas(figConv, fullfile(reportDir, 'convergence_curves.png'));
close(figConv);

% Save dataset snapshot for reproducibility.
save(fullfile(reportDir, 'dataset_snapshot.mat'), 'dataset');

fprintf('Comparison complete. Results available in %s\n', reportDir);

end
