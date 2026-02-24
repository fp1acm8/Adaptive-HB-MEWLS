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

% =====================================================================
%  Bootstrap: add src/ to the MATLAB path so +adaptivehb is resolvable
% =====================================================================

% Resolve the absolute path to this script's parent folder (scripts/).
rootDir  = fileparts(mfilename('fullpath'));

% Go one level up to reach the repository root directory.
repoRoot = fileparts(rootDir);

% Build the path to src/, which contains the +adaptivehb package tree.
srcDir   = fullfile(repoRoot, 'src');

% Snapshot the current MATLAB path before we modify it.
oldPath  = path;

% Register an onCleanup object that restores the original path when this
% function exits (normally or via error).  The %#ok<NASGU> pragma tells
% the MATLAB linter to suppress the "unused variable" warning — the
% variable must exist for the cleanup callback to fire on scope exit.
cleanup  = onCleanup(@() path(oldPath)); %#ok<NASGU>

% Prepend src/ to the path so all adaptivehb.* calls resolve correctly.
addpath(srcDir);

% Print a banner to the command window so the user knows the demo started.
fprintf('=== Adaptive-HB-MEWLS Quick Demo ===\n\n');

% =====================================================================
%  Step 1: Load a sample dataset
% =====================================================================

% Build the full path to the glacier.txt dataset shipped with the repo.
% This is a 3-column (x, y, f) file with no noise or outliers.
dataPath = fullfile(repoRoot, 'data', 'glacier.txt');

% Inform the user which file is being loaded (aids troubleshooting).
fprintf('Loading dataset: %s\n', dataPath);

% Call the IO loader with normalisation enabled.  Normalisation maps
% the x-y coordinates into [0,1]² via min-max scaling, which improves
% the numerical conditioning of the polynomial design matrix.
dataset = adaptivehb.io.load_dataset(dataPath, struct('normalize', true));

% Report the sample count so the user can verify the file was parsed.
fprintf('  Samples: %d\n', size(dataset.xy, 1));

% =====================================================================
%  Step 2: Inject Gaussian noise with 5% outliers
% =====================================================================

% Define a noise configuration struct.  Fields:
%   type              - "gaussian" adds N(0, sigma) noise to all samples
%   standardDeviation - sigma of the Gaussian perturbation
%   outlierFraction   - 5% of samples will receive inflated noise
%   outlierInflation  - outlier values are multiplied by 1 + 0.2 = 1.2
%   seed              - RNG seed for reproducibility (same seed → same noise)
noiseCfg = struct('type', 'gaussian', 'standardDeviation', 0.05, ...
    'outlierFraction', 0.05, 'outlierInflation', 0.2, 'seed', 42);

% Apply the noise to the dataset struct in-place.  After this call,
% dataset.fObserved differs from dataset.fTrue and dataset.isOutlier
% flags the corrupted samples.
dataset = adaptivehb.io.apply_noise(dataset, noiseCfg);

% Print a summary of the noise settings for the user's reference.
% The double %% is needed because fprintf interprets single % as a
% format specifier; %% produces a literal percent sign in the output.
fprintf('  Noise applied: Gaussian (sigma=%.2f, outliers=%.0f%%)\n', ...
    noiseCfg.standardDeviation, noiseCfg.outlierFraction * 100);

% =====================================================================
%  Step 3: Run both solvers (OLS and MEWLS) on the noisy dataset
% =====================================================================

% Maximum polynomial degree to fit.  The solvers sweep d = 1..maxDegree,
% building the design matrix Φ_d at each degree and recording metrics.
% degree 3 ⇒ 10 basis functions: M = (3+1)(3+2)/2 = 10.
maxDegree = 3;

% Inform the user before the (potentially slow) solver calls.
fprintf('\nRunning solvers with maxDegree = %d ...\n', maxDegree);

% --- OLS solver (Ordinary Least Squares) ---
% This is the unweighted baseline: minimises ||Φc - f||² via c = Φ\f.
% No entropy weights; all data points contribute equally.
lsResult = adaptivehb.solvers.leastSquaresSolver(dataset, ...
    struct('maxDegree', maxDegree));

% --- MEWLS solver (Maximum Entropy Weighted Least Squares) ---
% This solver computes entropy-based weights w_i = exp(-λ ||x_i - c̄||²)
% and solves the weighted system (Φ'WΦ)c = Φ'Wf.  The entropyScale
% parameter (λ = 15.0 here) controls how quickly the weights decay with
% distance from the centroid.  Ref: Brugnano et al. (2024), Section 3.
mewlsResult = adaptivehb.solvers.mewlsSolver(dataset, ...
    struct('maxDegree', maxDegree, 'entropyScale', 15.0));

% =====================================================================
%  Step 4: Print a formatted metrics comparison table
% =====================================================================

% Print the table header with fixed-width columns for alignment.
% %-18s: left-aligned solver name (18 chars wide)
% %12s:  right-aligned metric labels (12 chars wide each)
fprintf('\n%-18s %12s %12s %12s\n', 'Solver', 'RMSE', 'MaxAbsErr', 'MAE');

% Print a horizontal separator line (56 dashes matching the total width).
fprintf('%s\n', repmat('-', 1, 56));

% Print the OLS (LeastSquares) metrics row.
% metrics.rmse:        root mean squared error √(1/N · Σe²)
% metrics.maxAbsError: L∞ norm of the error vector max|e_i|
% metrics.mae:         mean absolute error 1/N · Σ|e_i|
fprintf('%-18s %12.6f %12.6f %12.6f\n', 'LeastSquares', ...
    lsResult.metrics.rmse, lsResult.metrics.maxAbsError, lsResult.metrics.mae);

% Print the MEWLS metrics row (same format for visual comparison).
fprintf('%-18s %12.6f %12.6f %12.6f\n', 'MEWLS', ...
    mewlsResult.metrics.rmse, mewlsResult.metrics.maxAbsError, mewlsResult.metrics.mae);

% =====================================================================
%  Step 5: Show the delta (MEWLS − OLS) for each metric
% =====================================================================

% The delta tells the user whether MEWLS improved or worsened each
% metric compared to the OLS baseline.  Negative delta = improvement
% (lower error), positive = regression.  The %+.6f format forces a
% sign character (+/-) before the number for clarity.
fprintf('\nDelta (MEWLS - LS):\n');
fprintf('  RMSE:        %+.6f\n', mewlsResult.metrics.rmse - lsResult.metrics.rmse);
fprintf('  MaxAbsError: %+.6f\n', mewlsResult.metrics.maxAbsError - lsResult.metrics.maxAbsError);
fprintf('  MAE:         %+.6f\n', mewlsResult.metrics.mae - lsResult.metrics.mae);

% =====================================================================
%  Footer: point the user to the full pipeline for richer output
% =====================================================================

% Signal that the demo completed without errors.
fprintf('\nDemo complete.\n');

% Suggest the next step: run_comparison produces plots, CSV, JSON, MAT.
fprintf('For full reports with plots, run: run_comparison()\n');

end
