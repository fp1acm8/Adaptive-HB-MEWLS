function result = mewlsSolver(dataset, params)
%MEWLSSOLVER Entropy-weighted polynomial fitter inspired by MEWLS.
%   RESULT = MEWLSSOLVER(DATASET, PARAMS) mirrors the polynomial sweep from
%   adaptivehb.solvers.leastSquaresSolver but reweights the samples with a
%   maximum-entropy inspired weight field. PARAMS.entropyScale controls how
%   rapidly weights decay with distance from the dataset centroid.

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

maxDegree = max(1, params.maxDegree);
scale = params.entropyScale;

weights = compute_entropy_weights(dataset.xy, scale, dataset.isOutlier);
W = diag(weights);

metricsHistory = zeros(maxDegree, 3);
degreeList = 1:maxDegree;

for d = degreeList
    [A, ~] = adaptivehb.solvers.polynomialDesignMatrix(dataset.xy, d);
    coeffs = (A' * W * A) \ (A' * W * dataset.fObserved);
    prediction = A * coeffs;
    metricsHistory(d, :) = compute_metrics(dataset.fTrue, prediction);
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

function weights = compute_entropy_weights(xy, scale, isOutlier)
centroid = mean(xy, 1);
dist2 = sum((xy - centroid) .^ 2, 2);
raw = exp(-scale * dist2);
raw(isOutlier) = raw(isOutlier) * 0.5; % penalise tagged outliers
weights = raw / sum(raw);
weights(weights <= 0) = eps;
end

function metrics = compute_metrics(fTrue, prediction)
residual = fTrue - prediction;
rmse = sqrt(mean(residual .^ 2));
maxAbs = max(abs(residual));
mae = mean(abs(residual));
metrics = [rmse, maxAbs, mae];
end
