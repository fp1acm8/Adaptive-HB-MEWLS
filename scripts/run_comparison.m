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
%   Nota: l'implementazione a B-spline gerarchiche adattive (HBS) che estende
%   questo confronto polinomiale richiede la libreria GeoPDEs per MATLAB.
%   Disponibile su: https://rafavzqz.github.io/geopdes/
%   Le funzioni GeoPDEs utilizzate nel codice HBS sono:
%     adaptivity_initialize_laplace, adaptivity_refine,
%     hspace_subdivision_matrix, op_gradgradu_gradgradv_hier,
%     sp_eval, hmsh_plot_cells, sp_get_cells.
%
%   See also adaptivehb.io.load_dataset, adaptivehb.io.apply_noise,
%            adaptivehb.solvers.leastSquaresSolver,
%            adaptivehb.solvers.mewlsSolver.

% =====================================================================
%  SECTION 1: Bootstrap — add src/ to the MATLAB path
% =====================================================================

% Determine the absolute path of this script file (scripts/ folder).
rootDir = fileparts(mfilename('fullpath'));

% Go one level up to reach the repository root.
repoRoot = fileparts(rootDir);

% The +adaptivehb package lives inside src/.
srcDir = fullfile(repoRoot, 'src');

% Save the current path so we can restore it when the function exits.
% onCleanup ensures restoration even if an error is thrown.
oldPath = path;
cleanupPath = onCleanup(@() path(oldPath)); %#ok<NASGU>

% Make the package visible to MATLAB for this session.
addpath(srcDir);

% =====================================================================
%  SECTION 2: Load and parse the JSON configuration
% =====================================================================

% If no config path was provided, fall back to the default config.
if nargin == 0 || isempty(configPath)
    configPath = fullfile('config', 'comparison_config.json');
end

% Verify that the configuration file exists before attempting to read it.
if ~isfile(configPath)
    error('adaptivehb:comparison:MissingConfig', ...
        'Configuration file not found: %s', configPath);
end

% Read the entire JSON file into a character vector and decode it into
% a MATLAB struct. The struct mirrors the JSON structure with fields:
%   config.dataset, config.noise, config.solvers, config.metrics.
configText = fileread(configPath);
config = jsondecode(configText);

% =====================================================================
%  SECTION 3: Prepare the dataset (load + optional noise injection)
% =====================================================================

% Load the point-cloud dataset from disk. The loader auto-detects whether
% the file has 3 columns (clean) or 5 columns (noisy with outlier flags).
% Coordinates are normalised to [0,1]^2 unless config.dataset.normalize=false.
% I dati sono nella forma (u,v,f(u,v)) come nel codice HBS originale.
dataset = adaptivehb.io.load_dataset(config.dataset.path, config.dataset);

% Inject synthetic noise/outliers according to the noise config section.
% If config.noise.type is "none", the dataset is returned unchanged.
dataset = adaptivehb.io.apply_noise(dataset, config.noise);

% =====================================================================
%  SECTION 4: Dispatch solvers
% =====================================================================

% Extract the solver specification array from the config.
% Each entry has: name, function (fully qualified), optional parameters.
solverSpecs = config.solvers;

% Pre-allocate the struct array that will hold each solver's output.
solverResults = repmat(struct(), 1, numel(solverSpecs));

for i = 1:numel(solverSpecs)
    % Convert the function name string (e.g. "adaptivehb.solvers.mewlsSolver")
    % into a callable function handle.
    solverFcn = str2func(solverSpecs(i).function);

    % Extract solver-specific parameters (method_data) if provided;
    % otherwise empty struct (defaults will be used inside the solver).
    if isfield(solverSpecs(i), 'parameters')
        solverParams = solverSpecs(i).parameters;
    else
        solverParams = struct();
    end

    % Call the solver inside a try-catch so that a single failing solver
    % does not abort the entire comparison. On failure, we fill the result
    % with NaN values and issue a warning.
    try
        solverResults(i) = solverFcn(dataset, solverParams);
        % Override the name field with the display name from the config
        % (the solver may set its own default name internally).
        solverResults(i).name = solverSpecs(i).name;
    catch ME
        warning('adaptivehb:comparison:SolverFailed', ...
            'Solver "%s" failed: %s', solverSpecs(i).name, ME.message);
        % Populate a placeholder result with NaN metrics so downstream
        % aggregation and plotting still work without crashing.
        solverResults(i).name = solverSpecs(i).name;
        solverResults(i).QI = NaN(size(dataset.f));
        solverResults(i).QI_coeff = [];
        solverResults(i).metrics = struct('rmse', NaN, 'maxAbsError', NaN, 'mae', NaN);
        solverResults(i).convergence = table();
    end
end

% =====================================================================
%  SECTION 5: Aggregate metrics into a summary table
% =====================================================================

% Read the list of metric names from the config (e.g. ["rmse","maxAbsError","mae"]).
metricNames = string(config.metrics(:));

% Build a matrix of size (nSolvers x nMetrics) from each solver's metrics struct.
metricsMatrix = zeros(numel(solverResults), numel(metricNames));
for i = 1:numel(solverResults)
    for j = 1:numel(metricNames)
        key = metricNames(j);
        % Dynamic field access: solverResults(i).metrics.(key)
        metricsMatrix(i, j) = solverResults(i).metrics.(key);
    end
end

% Convert the numeric matrix into a MATLAB table with named columns.
metricsTable = array2table(metricsMatrix, 'VariableNames', cellstr(metricNames));

% Add a 'solver' column with the display names and move it to the front.
solverNames = string({solverSpecs.name});
metricsTable.solver = solverNames';
metricsTable = movevars(metricsTable, 'solver', 'Before', 1);

% =====================================================================
%  SECTION 6: Compute deltas relative to the first (reference) solver
% =====================================================================

% The first solver in the config is treated as the reference baseline.
% Delta = metric_i - metric_reference; negative delta means improvement.
referenceRow = metricsMatrix(1, :);
deltaMatrix = metricsMatrix - referenceRow;

% Build a table with delta columns prefixed by "delta_".
deltaTable = array2table(deltaMatrix, 'VariableNames', strcat('delta_', cellstr(metricNames)));
deltaTable.solver = solverNames';
deltaTable = movevars(deltaTable, 'solver', 'Before', 1);

% Merge absolute metrics and deltas into one wide summary table.
% Skip the duplicate 'solver' column from deltaTable (column 1).
summaryTable = [metricsTable, deltaTable(:, 2:end)];

% =====================================================================
%  SECTION 7: Create the timestamped report directory
% =====================================================================

% Generate a timestamp string like "20240315_143022" for the folder name.
timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
reportDir = fullfile('reports', ['comparison_' timestamp]);

% Create the directory if it does not exist yet.
if ~exist(reportDir, 'dir')
    mkdir(reportDir);
end

% =====================================================================
%  SECTION 8: Persist metric tables (CSV, MAT, JSON)
% =====================================================================

% CSV — human-readable, easy to import in Excel or Python.
writetable(summaryTable, fullfile(reportDir, 'metrics.csv'));

% MAT — preserves full MATLAB table type for programmatic reloading.
save(fullfile(reportDir, 'metrics_table.mat'), 'summaryTable');

% JSON — machine-readable, useful for CI pipelines or web dashboards.
metricsStruct = table2struct(summaryTable);
jsonText = jsonencode(metricsStruct, 'PrettyPrint', true);
fid = fopen(fullfile(reportDir, 'metrics.json'), 'w');
fprintf(fid, '%s', jsonText);
fclose(fid);

% =====================================================================
%  SECTION 9: Generate and save diagnostic plots
% =====================================================================

% 3-D scatter plot comparing ground truth vs each solver's QI approximation.
figSurface = adaptivehb.viz.plot_solution_surface(dataset, solverResults);
saveas(figSurface, fullfile(reportDir, 'solution_surfaces.png'));
close(figSurface);

% Error histograms (PDF-normalised) for each solver.
figError = adaptivehb.viz.plot_error_distribution(dataset, solverResults);
saveas(figError, fullfile(reportDir, 'error_distributions.png'));
close(figError);

% RMSE convergence curves (RMSE vs polynomial degree) for each solver.
figConv = adaptivehb.viz.plot_convergence_curve(solverResults);
saveas(figConv, fullfile(reportDir, 'convergence_curves.png'));
close(figConv);

% =====================================================================
%  SECTION 10: Save dataset snapshot for full reproducibility
% =====================================================================

% Store the exact dataset struct (after normalisation and noise injection)
% so a future user can reload it and reproduce the solver results exactly.
save(fullfile(reportDir, 'dataset_snapshot.mat'), 'dataset');

% Print the path to the report directory so the user knows where to look.
fprintf('Comparison complete. Results available in %s\n', reportDir);

end
