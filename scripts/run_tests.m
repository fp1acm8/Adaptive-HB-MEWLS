function run_tests()
%RUN_TESTS Verification suite for the Adaptive-HB-MEWLS framework.
%   RUN_TESTS() executes a battery of non-interactive checks that validate
%   every module in the +adaptivehb package.  Each test prints PASS or FAIL
%   and a final summary is shown at the end.
%
%   Usage (from the repository root):
%       >> addpath('scripts');
%       >> run_tests
%
%   The script does NOT require any toolbox beyond base MATLAB (R2020b+).
%
%   Exit status:
%       The function returns normally when all tests pass.  When any test
%       fails, a warning is emitted but execution continues so that the
%       full report is always available.

% =====================================================================
%  Bootstrap: add src/ to the MATLAB path with auto-restore on exit
% =====================================================================
rootDir  = fileparts(mfilename('fullpath'));   % scripts/ folder
repoRoot = fileparts(rootDir);                 % repository root
srcDir   = fullfile(repoRoot, 'src');          % +adaptivehb package folder
oldPath  = path;                               % snapshot current path
cleanup  = onCleanup(@() path(oldPath));       %#ok<NASGU> restore on exit
addpath(srcDir);

% Cell array to collect all test results (structs with name/passed/message).
results = {};

fprintf('=== Adaptive-HB-MEWLS test suite ===\n\n');

% =====================================================================
%  Test 1: Package visibility
%  Verifies that +adaptivehb is on the MATLAB path and its functions
%  are resolvable.
% =====================================================================
results{end+1} = run_one('Package visibility', @test_package_visible);

% =====================================================================
%  Test 2: load_dataset — 3-column (clean) file
%  Loads glacier.txt, checks that the output struct has all required
%  fields and that normalisation maps coordinates to [0,1].
% =====================================================================
results{end+1} = run_one('load_dataset (3-col)', ...
    @() test_load_dataset(fullfile(repoRoot, 'data', 'glacier.txt')));

% =====================================================================
%  Test 3: load_dataset — 5-column (noisy) file
%  Loads a pre-generated noisy dataset and verifies that outlier flags
%  and the distinction between fTrue/fObserved are preserved.
% =====================================================================
fiveColFile = fullfile(repoRoot, 'data', ...
    'outl_black_forest_noise_gauss_frac0.05_int0.10_seed42.txt');
results{end+1} = run_one('load_dataset (5-col)', ...
    @() test_load_dataset_5col(fiveColFile));

% =====================================================================
%  Test 4: io.apply_noise
%  Applies Gaussian noise to a clean dataset and checks that fObserved
%  was actually modified and outliers were flagged.
% =====================================================================
results{end+1} = run_one('io.apply_noise', ...
    @() test_io_apply_noise(fullfile(repoRoot, 'data', 'glacier.txt')));

% =====================================================================
%  Test 5: data.apply_noise
%  Tests the vector-level noise injector with synthetic data.
% =====================================================================
results{end+1} = run_one('data.apply_noise', @test_data_apply_noise);

% =====================================================================
%  Test 6: data.format_noise_filename
%  Verifies filename generation for both noisy and clean metadata.
% =====================================================================
results{end+1} = run_one('data.format_noise_filename', @test_format_noise_filename);

% =====================================================================
%  Test 7: polynomialDesignMatrix
%  Checks that the design matrix has the correct dimensions for degree 2
%  (6 columns = (2+1)(2+2)/2) and the constant column is all ones.
% =====================================================================
results{end+1} = run_one('polynomialDesignMatrix', @test_design_matrix);

% =====================================================================
%  Test 8: compute_metrics
%  Validates RMSE, maxAbsError, MAE against manually computed values.
% =====================================================================
results{end+1} = run_one('compute_metrics', @test_compute_metrics);

% =====================================================================
%  Test 9: leastSquaresSolver
%  Runs the OLS solver on synthetic data and checks the output struct.
% =====================================================================
results{end+1} = run_one('leastSquaresSolver', ...
    @() test_solver(@adaptivehb.solvers.leastSquaresSolver, 'LeastSquares'));

% =====================================================================
%  Test 10: mewlsSolver
%  Runs the MEWLS solver on synthetic data and checks the output struct.
% =====================================================================
results{end+1} = run_one('mewlsSolver', ...
    @() test_solver(@adaptivehb.solvers.mewlsSolver, 'MEWLS'));

% =====================================================================
%  Test 11: Visualization
%  Creates and closes each plot type to verify no runtime errors.
% =====================================================================
results{end+1} = run_one('Visualization', @test_visualizations);

% =====================================================================
%  Test 12: run_comparison end-to-end
%  Runs the full pipeline with the default config and checks that all
%  expected output files were created in the reports/ directory.
% =====================================================================
results{end+1} = run_one('run_comparison (e2e)', ...
    @() test_run_comparison(repoRoot));

% =====================================================================
%  Summary: count and display pass/fail totals
% =====================================================================
fprintf('\n=== Summary ===\n');
nPass = 0; nFail = 0;
for k = 1:numel(results)
    r = results{k};
    if r.passed
        nPass = nPass + 1;
        fprintf('  PASS  %s\n', r.name);
    else
        nFail = nFail + 1;
        fprintf('  FAIL  %s : %s\n', r.name, r.message);
    end
end
fprintf('\n%d passed, %d failed out of %d tests.\n', nPass, nFail, nPass + nFail);

% Emit a warning (not an error) so the full report is always visible.
if nFail > 0
    warning('adaptivehb:tests:SomeTestsFailed', ...
        '%d test(s) failed.  See output above for details.', nFail);
end

end


% =====================================================================
%  Test runner helper: run_one
%  Wraps a test function in try-catch and returns a result struct.
% =====================================================================
function result = run_one(name, fcn)
try
    fcn();  % execute the test function
    result = struct('name', name, 'passed', true, 'message', '');
    fprintf('  PASS  %s\n', name);
catch ME
    result = struct('name', name, 'passed', false, 'message', ME.message);
    fprintf('  FAIL  %s : %s\n', name, ME.message);
end
end


% =====================================================================
%  Individual test functions
%  Each function uses assert() to verify expected conditions.
%  If an assert fails, MATLAB throws an error caught by run_one().
% =====================================================================

function test_package_visible()
% Verify that the +adaptivehb package is on the path and its primary
% function can be resolved by MATLAB's which() command.
assert(exist('adaptivehb.io.load_dataset', 'file') > 0 || ...
       ~isempty(which('adaptivehb.io.load_dataset')), ...
    'Package +adaptivehb not found on path.');
end

function test_load_dataset(filePath)
% Load a 3-column dataset and verify the output struct fields.
ds = adaptivehb.io.load_dataset(filePath, struct('normalize', true));
assert(isstruct(ds), 'Expected struct output.');

% All five required fields must be present.
requiredFields = {'originalXY','xy','fTrue','fObserved','isOutlier'};
for k = 1:numel(requiredFields)
    assert(isfield(ds, requiredFields{k}), ...
        sprintf('Missing field: %s', requiredFields{k}));
end

% xy must have exactly 2 columns (x, y).
assert(size(ds.xy, 2) == 2, 'xy must have 2 columns.');

% Number of samples must be consistent across fields.
assert(numel(ds.fTrue) == size(ds.xy, 1), 'Dimension mismatch.');

% After normalisation, all coordinates must lie in [0,1].
assert(all(ds.xy(:) >= -eps & ds.xy(:) <= 1+eps), ...
    'Normalised coordinates outside [0,1].');
end

function test_load_dataset_5col(filePath)
% Load a 5-column (noisy) dataset and verify outlier-specific fields.
ds = adaptivehb.io.load_dataset(filePath);

% At least some points should be flagged as outliers.
assert(any(ds.isOutlier), 'Expected some outliers in 5-column dataset.');

% fTrue and fObserved must differ (the noise separates them).
assert(~isequal(ds.fTrue, ds.fObserved), ...
    'fTrue and fObserved should differ in noisy dataset.');
end

function test_io_apply_noise(filePath)
% Apply Gaussian noise to a clean dataset and verify the effect.
ds = adaptivehb.io.load_dataset(filePath);

noiseCfg = struct('type', 'gaussian', 'standardDeviation', 0.1, ...
    'outlierFraction', 0.1, 'outlierInflation', 0.5, 'seed', 99);
dsNoisy = adaptivehb.io.apply_noise(ds, noiseCfg);

% fObserved must have changed after noise injection.
assert(~isequal(ds.fObserved, dsNoisy.fObserved), ...
    'Noise should change fObserved.');

% Some points must be flagged as outliers.
assert(any(dsNoisy.isOutlier), 'Expected some outliers after noise.');
end

function test_data_apply_noise()
% Test the vector-level noise injector with synthetic sinusoidal data.
rng(1);
x = linspace(0, 1, 100)';
y = linspace(0, 1, 100)';
f = sin(pi*x) .* cos(pi*y);

settings = struct('addNoise', true, 'outlierFraction', 0.1, ...
    'outlierIntensity', 0.2, 'seed', 42, 'noiseType', 'gauss');
[aug, meta] = adaptivehb.data.apply_noise(x, y, f, settings);

assert(meta.noiseApplied, 'noiseApplied should be true.');
assert(meta.numOutliers > 0, 'Expected some outliers.');
assert(~isequal(aug.fTrue, aug.fNoisy), 'fNoisy should differ from fTrue.');
assert(numel(aug.x) == 100, 'Output size mismatch.');
end

function test_format_noise_filename()
% Test filename generation for a noisy metadata struct.
meta = struct('noiseApplied', true, 'noiseLabel', 'noise_gauss_frac0.10_int0.20_seed42');
name = adaptivehb.data.format_noise_filename('testdata', meta);
assert(contains(name, 'noise_gauss'), 'Filename should encode noise type.');
assert(endsWith(name, '.txt'), 'Default extension should be .txt.');

% Test filename generation for a clean (no noise) metadata struct.
metaClean = struct('noiseApplied', false);
nameClean = adaptivehb.data.format_noise_filename('testdata', metaClean);
assert(strcmp(nameClean, 'testdata.txt'), 'Clean filename should be basename.txt.');
end

function test_design_matrix()
% Verify design matrix dimensions for degree 2.
xy = rand(50, 2);
[A, terms] = adaptivehb.solvers.polynomialDesignMatrix(xy, 2);

% For degree 2: M = (2+1)(2+2)/2 = 6 columns.
assert(size(A, 1) == 50, 'Row count must match input points.');
assert(size(A, 2) == 6, 'Degree-2 should produce 6 columns.');
assert(size(terms, 1) == 6, 'terms should have 6 rows.');

% The first column (constant term, p=q=0) must be all ones.
assert(all(abs(A(:, 1) - 1) < eps), 'First column should be ones.');
end

function test_compute_metrics()
% Validate the three metrics against hand-computed expected values.
fTrue = [1; 2; 3; 4; 5];
pred  = [1.1; 2.0; 2.8; 4.2; 4.9];
m = adaptivehb.solvers.compute_metrics(fTrue, pred);

assert(numel(m) == 3, 'compute_metrics must return 3 values.');

residual = fTrue - pred;
% RMSE = sqrt(mean(residual.^2))
expectedRMSE = sqrt(mean(residual.^2));
assert(abs(m(1) - expectedRMSE) < 1e-12, 'RMSE mismatch.');
% maxAbsError = max(|residual|)
assert(abs(m(2) - max(abs(residual))) < 1e-12, 'maxAbsError mismatch.');
% MAE = mean(|residual|)
assert(abs(m(3) - mean(abs(residual))) < 1e-12, 'MAE mismatch.');
end

function test_solver(solverFcn, expectedName)
% Run a solver on synthetic data and verify the output struct fields.
rng(0);
N = 200;
xy = rand(N, 2);
fTrue = sin(pi*xy(:,1)) .* cos(pi*xy(:,2));

% Build a dataset struct with no noise (fObserved = fTrue).
ds = struct('xy', xy, 'fTrue', fTrue, 'fObserved', fTrue, ...
    'isOutlier', false(N,1));
params = struct('maxDegree', 2);
result = solverFcn(ds, params);

% Verify required output fields.
assert(isfield(result, 'prediction'), 'Missing field: prediction.');
assert(isfield(result, 'metrics'), 'Missing field: metrics.');
assert(isfield(result, 'convergence'), 'Missing field: convergence.');

% Prediction must have N elements.
assert(numel(result.prediction) == N, 'prediction length mismatch.');

% RMSE must be non-negative.
assert(result.metrics.rmse >= 0, 'RMSE must be non-negative.');

% Convergence table must have 2 rows (one per degree: 1 and 2).
assert(height(result.convergence) == 2, 'convergence should have 2 rows for maxDegree=2.');
end

function test_visualizations()
% Create and immediately close each plot type to verify no runtime errors.
rng(0);
N = 50;
xy = rand(N, 2);
fTrue = xy(:,1) + xy(:,2);

% Build a minimal dataset and solver result for plotting.
ds = struct('xy', xy, 'fTrue', fTrue, 'fObserved', fTrue, ...
    'isOutlier', false(N,1));
sr = struct('name', "TestSolver", 'prediction', fTrue + 0.01*randn(N,1), ...
    'convergence', table((1:2)', rand(2,1), rand(2,1), rand(2,1), ...
    'VariableNames', {'degree','rmse','maxAbsError','mae'}));

% Test each visualisation function: must return a valid figure handle.
fig1 = adaptivehb.viz.plot_solution_surface(ds, sr);
assert(ishandle(fig1), 'plot_solution_surface did not return a valid handle.');
close(fig1);

fig2 = adaptivehb.viz.plot_error_distribution(ds, sr);
assert(ishandle(fig2), 'plot_error_distribution did not return a valid handle.');
close(fig2);

fig3 = adaptivehb.viz.plot_convergence_curve(sr);
assert(ishandle(fig3), 'plot_convergence_curve did not return a valid handle.');
close(fig3);
end

function test_run_comparison(repoRoot)
% Run the full comparison pipeline with the default configuration
% and verify that all expected output files are created.

configPath = fullfile(repoRoot, 'config', 'comparison_config.json');
assert(isfile(configPath), 'Default config not found.');

% Execute run_comparison — it creates a reports/<timestamp>/ subfolder.
run_comparison(configPath);

% Verify that the reports/ directory was created.
reportBase = fullfile(repoRoot, 'reports');
assert(isfolder(reportBase), 'reports/ directory was not created.');

% Find the most recently created comparison report directory.
listing = dir(fullfile(reportBase, 'comparison_*'));
assert(~isempty(listing), 'No comparison report directory found.');

% Check that all seven expected output files exist in the latest report.
latestReport = fullfile(reportBase, listing(end).name);
expectedFiles = {'metrics.csv', 'metrics.json', 'metrics_table.mat', ...
    'solution_surfaces.png', 'error_distributions.png', ...
    'convergence_curves.png', 'dataset_snapshot.mat'};
for k = 1:numel(expectedFiles)
    assert(isfile(fullfile(latestReport, expectedFiles{k})), ...
        sprintf('Missing output: %s', expectedFiles{k}));
end
end
