function result = leastSquaresSolver(dataset, params)
%LEASTSQUARESSOLVER Polynomial surrogate using Ordinary Least Squares.
%   RESULT = LEASTSQUARESSOLVER(DATASET, PARAMS) fits a polynomial surface
%   to the observed samples using ordinary least squares (OLS).  The solver
%   sweeps polynomial degrees from 1 to PARAMS.maxDegree and retains the
%   full convergence history so that the user can inspect how the error
%   evolves with increasing polynomial complexity.
%
%   This solver is the unweighted baseline (W = I) for comparison against
%   the entropy-weighted MEWLS variant. In the notation of Brugnano et al.
%   (2024), setting all weights w_i = 1/N reduces the MEWLS formulation
%   to standard OLS. The comparison quantifies the improvement obtained by
%   the MEWLS weighting strategy.
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
%   See also adaptivehb.solvers.mewlsSolver,
%            adaptivehb.solvers.polynomialDesignMatrix,
%            adaptivehb.solvers.compute_metrics.

% --- Input validation (MATLAB R2020b+ arguments block) ---------------
% dataset must be a scalar struct; params defaults to empty struct so the
% caller can omit it entirely.
arguments
    dataset (1, 1) struct
    params (1, 1) struct = struct()
end

% --- Default parameter handling ---------------------------------------

% maxDegree (d): highest total polynomial degree d used in the sweep.
% For degree d the bivariate polynomial space Pi_d^2 has
% M = (d+1)(d+2)/2 basis functions. Default d=3 -> M=10.
if ~isfield(params, 'maxDegree') || isempty(params.maxDegree)
    params.maxDegree = 3;
end

% Clamp to at least 1 (degree 0 = constant is trivial).
maxDegree = max(1, params.maxDegree);

% --- Polynomial degree sweep ------------------------------------------
% Fit polynomials of total degree d = 1, 2, ..., maxDegree and record
% error metrics at each step to produce a convergence curve.
metricsHistory = zeros(maxDegree, 3); % each row: [RMSE, maxAbsErr, MAE]
degreeList = 1:maxDegree;             % degrees to evaluate

for d = degreeList
    % Build the design matrix Phi_d in R^{N x M_d} for the bivariate
    % polynomial space Pi_d^2.  Each column is x^p * y^q with p+q <= d.
    % Ref: same Phi matrix used in MEWLS (Brugnano et al. 2024, Sec. 3),
    % but here solved without the weight matrix W.
    [A, ~] = adaptivehb.solvers.polynomialDesignMatrix(dataset.xy, d);

    % Solve the OLS problem via MATLAB backslash operator:
    %   c = Phi \ f
    % This minimises ||Phi*c - f||_2^2 using QR factorisation internally.
    % Equivalent to the MEWLS normal equations with W = I (uniform weights).
    coeffs = A \ dataset.fObserved;

    % Evaluate the fitted polynomial at all data points: f_hat = Phi * c.
    prediction = A * coeffs;

    % Compute RMSE, maxAbsError, MAE between true values and prediction.
    metricsHistory(d, :) = adaptivehb.solvers.compute_metrics(dataset.fTrue, prediction);

    % Store coefficient vector and prediction for this degree.
    coeffHistory{d} = coeffs; %#ok<AGROW>
    predictionHistory{d} = prediction; %#ok<AGROW>
end

% --- Assemble output struct -------------------------------------------

% Use results from the highest polynomial degree as the final output.
finalPrediction = predictionHistory{end};

% Pack the three metric values into a named struct for easy access.
metricNames = {"rmse", "maxAbsError", "mae"};
metricsStruct = cell2struct(num2cell(metricsHistory(end, :)), metricNames, 2);

% Build a convergence table (degree vs metrics) for plotting/analysis.
convergence = table();
convergence.degree = degreeList';       % polynomial degrees 1..maxDegree
convergence.rmse = metricsHistory(:, 1);
convergence.maxAbsError = metricsHistory(:, 2);
convergence.mae = metricsHistory(:, 3);

% Return all results. Note: unlike mewlsSolver, no 'weights' field is
% included because OLS uses uniform (implicit) weights.
result = struct(...
    'name', "LeastSquares", ...
    'prediction', finalPrediction, ...       % N-by-1 fitted values
    'coefficients', coeffHistory{end}, ...   % M-by-1 polynomial coeffs
    'metrics', metricsStruct, ...            % struct with rmse/maxAbsError/mae
    'convergence', convergence);             % table: degree x metrics

end
