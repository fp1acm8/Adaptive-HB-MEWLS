function result = mewlsSolver(dataset, params)
%MEWLSSOLVER Maximum Entropy Weighted Least Squares polynomial fitter.
%   RESULT = MEWLSSOLVER(DATASET, PARAMS) fits a polynomial surface using
%   Weighted Least Squares (WLS) where sample weights are derived from a
%   maximum-entropy inspired criterion.  Like the OLS baseline
%   (leastSquaresSolver), the solver sweeps polynomial degrees from 1 to
%   PARAMS.maxDegree and records the full convergence history.
%
%   Algorithm
%   ---------
%   1. Compute the centroid of the input point cloud.
%   2. Assign each sample an entropy weight:
%          w_i = exp(-entropyScale * ||p_i - centroid||^2)
%      Low entropyScale -> nearly uniform weights (similar to OLS).
%      High entropyScale -> weights concentrated near the centroid.
%   3. Penalise samples flagged as outliers by multiplying their weight
%      by outlierPenalty (default 0.5), reducing their influence.
%   4. Normalise: w = w / sum(w), then clamp to eps for numerical safety.
%   5. Solve the WLS problem: coeffs = (A' W A) \ (A' W f).
%
%   Inputs:
%       DATASET  - struct from adaptivehb.io.load_dataset with fields
%                  xy (N-by-2), fObserved (N-by-1), fTrue, isOutlier.
%       PARAMS   - struct with optional fields:
%                    maxDegree      : max polynomial degree (default 3)
%                    entropyScale   : controls weight decay rate (default 10).
%                                    Typical range: 1 (mild) to 50 (strong).
%                    outlierPenalty : multiplicative factor applied to
%                                    outlier weights (default 0.5, i.e.
%                                    outliers contribute half as much).
%
%   Output:
%       RESULT - struct with fields:
%           name         - "MEWLS"
%           prediction   - N-by-1 fitted values (final degree)
%           coefficients - polynomial coefficients (final degree)
%           metrics      - struct with rmse, maxAbsError, mae
%           convergence  - table with degree-by-degree error history
%           weights      - N-by-1 normalised sample weights
%
%   See also adaptivehb.solvers.leastSquaresSolver,
%            adaptivehb.solvers.polynomialDesignMatrix,
%            adaptivehb.solvers.compute_metrics.

arguments
    dataset (1, 1) struct
    params (1, 1) struct = struct()
end

if ~isfield(params, 'maxDegree') || isempty(params.maxDegree)
    params.maxDegree = 3;
end
if ~isfield(params, 'entropyScale') || isempty(params.entropyScale)
    params.entropyScale = 10.0;
end
if ~isfield(params, 'outlierPenalty') || isempty(params.outlierPenalty)
    params.outlierPenalty = 0.5;
end

maxDegree = max(1, params.maxDegree);
scale = params.entropyScale;
outlierPenalty = params.outlierPenalty;

weights = compute_entropy_weights(dataset.xy, scale, dataset.isOutlier, outlierPenalty);
W = diag(weights);

metricsHistory = zeros(maxDegree, 3);
degreeList = 1:maxDegree;

for d = degreeList
    [A, ~] = adaptivehb.solvers.polynomialDesignMatrix(dataset.xy, d);
    % Solve the weighted least squares normal equations.
    coeffs = (A' * W * A) \ (A' * W * dataset.fObserved);
    prediction = A * coeffs;
    metricsHistory(d, :) = adaptivehb.solvers.compute_metrics(dataset.fTrue, prediction);
    coeffHistory{d} = coeffs; %#ok<AGROW>
    predictionHistory{d} = prediction; %#ok<AGROW>
end

finalPrediction = predictionHistory{end};
metricNames = {"rmse", "maxAbsError", "mae"};
metricsStruct = cell2struct(num2cell(metricsHistory(end, :)), metricNames, 2);

convergence = table();
convergence.degree = degreeList';
convergence.rmse = metricsHistory(:, 1);
convergence.maxAbsError = metricsHistory(:, 2);
convergence.mae = metricsHistory(:, 3);

result = struct(...
    'name', "MEWLS", ...
    'prediction', finalPrediction, ...
    'coefficients', coeffHistory{end}, ...
    'metrics', metricsStruct, ...
    'convergence', convergence, ...
    'weights', weights);

end

function weights = compute_entropy_weights(xy, scale, isOutlier, outlierPenalty)
%COMPUTE_ENTROPY_WEIGHTS Build normalised entropy-based sample weights.
%   Weights decay exponentially with squared distance from the centroid.
%   Samples flagged as outliers are down-weighted by outlierPenalty.
centroid = mean(xy, 1);
dist2 = sum((xy - centroid) .^ 2, 2);
raw = exp(-scale * dist2);
raw(isOutlier) = raw(isOutlier) * outlierPenalty;
weights = raw / sum(raw);
weights(weights <= 0) = eps;
end
