function result = leastSquaresSolver(dataset, params)
%LEASTSQUARESSOLVER Polynomial surrogate for the reference LS workflow.
%   RESULT = LEASTSQUARESSOLVER(DATASET, PARAMS) fits a polynomial surface
%   to the observed samples using ordinary least squares. PARAMS.maxDegree
%   controls the maximum polynomial degree evaluated during the sweep.
%
%   The returned struct exposes:
%       - name        : solver identifier
%       - prediction  : fitted values at the input points
%       - metrics     : struct with rmse/maxAbsError/mae
%       - convergence : table with per-degree error trends
%
%   The simplified solver provides a deterministic baseline that can be
%   compared against the entropy-weighted variant defined in
%   adaptivehb.solvers.mewlsSolver.

arguments
    dataset (1, 1) struct
    params (1, 1) struct = struct()
end

if ~isfield(params, 'maxDegree') || isempty(params.maxDegree)
    params.maxDegree = 3;
end

maxDegree = params.maxDegree;
maxDegree = max(1, maxDegree);

metricsHistory = zeros(maxDegree, 3);
degreeList = 1:maxDegree;

for d = degreeList
    [A, ~] = adaptivehb.solvers.polynomialDesignMatrix(dataset.xy, d);
    coeffs = A \ dataset.fObserved;
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
    'name', "LeastSquares", ...
    'prediction', finalPrediction, ...
    'coefficients', coeffHistory{end}, ...
    'metrics', metricsStruct, ...
    'convergence', convergence);

end

function metrics = compute_metrics(fTrue, prediction)
residual = fTrue - prediction;
rmse = sqrt(mean(residual .^ 2));
maxAbs = max(abs(residual));
mae = mean(abs(residual));
metrics = [rmse, maxAbs, mae];
end
