function run_quick_demo()
%RUN_QUICK_DEMO Non-interactive demonstration of the MEWLS framework.
%   RUN_QUICK_DEMO() loads a sample dataset, injects Gaussian noise,
%   runs both the OLS and MEWLS solvers, and prints a side-by-side
%   comparison of their error metrics to the command window.
%
%   This script is intended as a 30-second orientation for new users who
%   want to see the framework in action without editing any configuration.
%
%   Usage (from the repository root):
%       >> addpath('scripts');
%       >> run_quick_demo
%
%   See also run_comparison, run_tests.

% ---- bootstrap paths ----
rootDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(rootDir);
srcDir   = fullfile(repoRoot, 'src');
oldPath  = path;
cleanup  = onCleanup(@() path(oldPath)); %#ok<NASGU>
addpath(srcDir);

fprintf('=== Adaptive-HB-MEWLS Quick Demo ===\n\n');

% 1. Load a sample dataset.
dataPath = fullfile(repoRoot, 'data', 'glacier.txt');
fprintf('Loading dataset: %s\n', dataPath);
dataset = adaptivehb.io.load_dataset(dataPath, struct('normalize', true));
fprintf('  Samples: %d\n', size(dataset.xy, 1));

% 2. Inject Gaussian noise with 5%% outliers.
noiseCfg = struct('type', 'gaussian', 'standardDeviation', 0.05, ...
    'outlierFraction', 0.05, 'outlierInflation', 0.2, 'seed', 42);
dataset = adaptivehb.io.apply_noise(dataset, noiseCfg);
fprintf('  Noise applied: Gaussian (sigma=%.2f, outliers=%.0f%%)\n', ...
    noiseCfg.standardDeviation, noiseCfg.outlierFraction * 100);

% 3. Run solvers.
maxDegree = 3;
fprintf('\nRunning solvers with maxDegree = %d ...\n', maxDegree);

lsResult = adaptivehb.solvers.leastSquaresSolver(dataset, ...
    struct('maxDegree', maxDegree));

mewlsResult = adaptivehb.solvers.mewlsSolver(dataset, ...
    struct('maxDegree', maxDegree, 'entropyScale', 15.0));

% 4. Print metrics comparison.
fprintf('\n%-18s %12s %12s %12s\n', 'Solver', 'RMSE', 'MaxAbsErr', 'MAE');
fprintf('%s\n', repmat('-', 1, 56));
fprintf('%-18s %12.6f %12.6f %12.6f\n', 'LeastSquares', ...
    lsResult.metrics.rmse, lsResult.metrics.maxAbsError, lsResult.metrics.mae);
fprintf('%-18s %12.6f %12.6f %12.6f\n', 'MEWLS', ...
    mewlsResult.metrics.rmse, mewlsResult.metrics.maxAbsError, mewlsResult.metrics.mae);

% 5. Show delta.
fprintf('\nDelta (MEWLS - LS):\n');
fprintf('  RMSE:        %+.6f\n', mewlsResult.metrics.rmse - lsResult.metrics.rmse);
fprintf('  MaxAbsError: %+.6f\n', mewlsResult.metrics.maxAbsError - lsResult.metrics.maxAbsError);
fprintf('  MAE:         %+.6f\n', mewlsResult.metrics.mae - lsResult.metrics.mae);

fprintf('\nDemo complete.\n');
fprintf('For full reports with plots, run: run_comparison()\n');

end
