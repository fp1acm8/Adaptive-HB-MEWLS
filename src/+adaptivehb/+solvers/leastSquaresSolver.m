function result = leastSquaresSolver(dataset, params)
%LEASTSQUARESSOLVER Polynomial surrogate using Ordinary Least Squares.
%   RESULT = LEASTSQUARESSOLVER(DATASET, PARAMS) fits a polynomial surface
%   to the observed samples using ordinary least squares (OLS).  The solver
%   sweeps polynomial degrees from 1 to PARAMS.maxDegree and retains the
%   full convergence history so that the user can inspect how the error
%   evolves with increasing polynomial complexity.
%
%   Inputs:
%       DATASET  - struct produced by adaptivehb.io.load_dataset containing
%                  at least xy (N-by-2), fObserved (N-by-1), and fTrue.
%       PARAMS   - struct with optional fields:
%                    maxDegree : maximum polynomial degree (default 3).
%
%   Output:
%       RESULT - struct with the following fields:
%           name         - solver identifier string ("LeastSquares")
%           prediction   - N-by-1 fitted values at input points (degree = maxDegree)
%           coefficients - polynomial coefficient vector for the final degree
%           metrics      - struct with rmse, maxAbsError, mae for the final degree
%           convergence  - table with columns degree, rmse, maxAbsError, mae
%
%   This solver provides a deterministic baseline for comparison against
%   the entropy-weighted MEWLS variant (adaptivehb.solvers.mewlsSolver).
%
%   See also adaptivehb.solvers.mewlsSolver,
%            adaptivehb.solvers.polynomialDesignMatrix,
%            adaptivehb.solvers.compute_metrics.

arguments
    dataset (1, 1) struct
    params (1, 1) struct = struct()
end

if ~isfield(params, 'maxDegree') || isempty(params.maxDegree)
    params.maxDegree = 3;
end

maxDegree = max(1, params.maxDegree);

metricsHistory = zeros(maxDegree, 3);
degreeList = 1:maxDegree;

for d = degreeList
    [A, ~] = adaptivehb.solvers.polynomialDesignMatrix(dataset.xy, d);
    coeffs = A \ dataset.fObserved;
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
    'name', "LeastSquares", ...
    'prediction', finalPrediction, ...
    'coefficients', coeffHistory{end}, ...
    'metrics', metricsStruct, ...
    'convergence', convergence);

end
