function result = mewlsSolver(dataset, params)
%MEWLSSOLVER Maximum Entropy Weighted Least Squares polynomial fitter.
%   RESULT = MEWLSSOLVER(DATASET, PARAMS) fits a polynomial surface using
%   Weighted Least Squares (WLS) where sample weights are derived from a
%   maximum-entropy inspired criterion.  Like the OLS baseline
%   (leastSquaresSolver), the solver sweeps polynomial degrees from 1 to
%   PARAMS.maxDegree and records the full convergence history.
%
%   Reference
%   ---------
%   Based on the MEWLS formulation described in:
%       Brugnano, Ferretti, Lucci, Pucci (2024),
%       "Adaptive refinement with a locally weighted least-squares
%        approach based on maximum entropy."
%   The weight function and WLS assembly follow the notation of that paper.
%
%   Algorithm (see Brugnano et al. 2024, Section 3)
%   ---------
%   1. Compute the centroid c of the input point cloud.
%   2. Assign each sample an entropy weight:
%          w_i = exp(-lambda * ||p_i - c||^2)       [Eq. in Sec. 3]
%      Low lambda  -> nearly uniform weights (approaches OLS).
%      High lambda -> weights concentrated near the centroid.
%   3. Penalise samples flagged as outliers by multiplying their weight
%      by outlierPenalty (default 0.5), reducing their influence.
%   4. Normalise: w = w / sum(w), then clamp to eps for numerical safety.
%   5. Build the design matrix Phi for polynomial space Pi_d^2.
%   6. Solve the WLS normal equations:
%          (Phi' W Phi) c = Phi' W f                [Eq. in Sec. 3]
%
%   Inputs:
%       DATASET  - struct from adaptivehb.io.load_dataset with fields
%                  xy (N-by-2), fObserved (N-by-1), fTrue, isOutlier.
%       PARAMS   - struct with optional fields:
%                    maxDegree      : max polynomial degree (default 3)
%                    entropyScale   : lambda, controls weight decay rate
%                                    (default 10). Range: 1..50.
%                    outlierPenalty : multiplicative factor applied to
%                                    outlier weights (default 0.5).
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

% --- Input validation (MATLAB R2020b+ arguments block) ---------------
% dataset must be a scalar struct; params defaults to empty struct so that
% all parameters below become optional.
arguments
    dataset (1, 1) struct
    params (1, 1) struct = struct()
end

% --- Default parameter handling ---------------------------------------

% maxDegree (d): highest total polynomial degree used in the sweep over
% the bivariate polynomial space Pi_d^2 (Brugnano et al. 2024, Sec. 3).
% Default d=3 gives M=(3+1)(3+2)/2 = 10 basis functions.
if ~isfield(params, 'maxDegree') || isempty(params.maxDegree)
    params.maxDegree = 3;
end

% entropyScale (lambda in the paper): controls how fast the MEWLS weights
% decay with distance from the centroid.
%   lambda ~ 1-5   : mild decay, weights nearly uniform (close to OLS)
%   lambda ~ 10-20 : moderate decay (recommended starting range)
%   lambda > 30    : aggressive decay, strong centroid bias
% Ref: Brugnano et al. 2024, weight function w_i = exp(-lambda * d_i^2).
if ~isfield(params, 'entropyScale') || isempty(params.entropyScale)
    params.entropyScale = 10.0;
end

% outlierPenalty (alpha_out): multiplicative factor for outlier weights.
%   0.0 -> outliers completely ignored
%   0.5 -> outliers contribute at half weight (default)
%   1.0 -> no penalty (outlier flag has no effect)
if ~isfield(params, 'outlierPenalty') || isempty(params.outlierPenalty)
    params.outlierPenalty = 0.5;
end

% Clamp maxDegree to at least 1 (degree 0 = constant is trivial).
maxDegree = max(1, params.maxDegree);

% Store lambda and penalty in local variables for readability.
scale = params.entropyScale;            % lambda in the paper
outlierPenalty = params.outlierPenalty;  % alpha_out

% --- Entropy weight computation (Brugnano et al. 2024, Sec. 3) -------
% Build the MEWLS weight vector w (N-by-1):
%   w_i = exp(-lambda * ||x_i - centroid||^2)
% then penalise outliers and normalise so sum(w) = 1.
weights = compute_entropy_weights(dataset.xy, scale, dataset.isOutlier, outlierPenalty);

% Construct the diagonal weight matrix W = diag(w_1, ..., w_N).
% This appears as W in the WLS normal equations: (Phi' W Phi) c = Phi' W f.
W = diag(weights);

% --- Polynomial degree sweep ------------------------------------------
% Fit polynomials of total degree d = 1, 2, ..., maxDegree and record
% error metrics at each step. This provides a convergence curve showing
% how the approximation improves with richer polynomial spaces Pi_d^2.
metricsHistory = zeros(maxDegree, 3); % each row: [RMSE, maxAbsErr, MAE]
degreeList = 1:maxDegree;             % degrees to evaluate

for d = degreeList
    % Build the design matrix Phi_d in R^{N x M_d} for the bivariate
    % polynomial space Pi_d^2 where M_d = (d+1)(d+2)/2 columns.
    % Each column is one monomial x^p * y^q with p+q <= d.
    % Ref: Brugnano et al. 2024, polynomial basis construction.
    [A, ~] = adaptivehb.solvers.polynomialDesignMatrix(dataset.xy, d);

    % Solve the WLS normal equations (Brugnano et al. 2024, Sec. 3):
    %   c = (Phi' W Phi) \ (Phi' W f)
    % This minimises the weighted residual: sum_i w_i*(f_i - p(x_i,y_i))^2
    % where p is the polynomial of degree d and w_i are the MEWLS weights.
    coeffs = (A' * W * A) \ (A' * W * dataset.fObserved);

    % Evaluate the fitted polynomial at all data points: f_hat = Phi * c.
    prediction = A * coeffs;

    % Compute RMSE, maxAbsError, MAE between true values and prediction.
    % These quantify the approximation quality for this degree.
    metricsHistory(d, :) = adaptivehb.solvers.compute_metrics(dataset.fTrue, prediction);

    % Store coefficient vector and prediction for this degree so we can
    % retrieve the final-degree results after the loop ends.
    coeffHistory{d} = coeffs; %#ok<AGROW>
    predictionHistory{d} = prediction; %#ok<AGROW>
end

% --- Assemble output struct -------------------------------------------

% Use results from the highest polynomial degree as the final output.
finalPrediction = predictionHistory{end};

% Pack the three metric values from the last degree into a named struct
% for convenient access: result.metrics.rmse, result.metrics.maxAbsError, etc.
metricNames = {"rmse", "maxAbsError", "mae"};
metricsStruct = cell2struct(num2cell(metricsHistory(end, :)), metricNames, 2);

% Build a convergence table so the caller can plot RMSE vs degree or
% inspect how each metric evolves across the polynomial degree sweep.
convergence = table();
convergence.degree = degreeList';       % polynomial degrees 1..maxDegree
convergence.rmse = metricsHistory(:, 1);
convergence.maxAbsError = metricsHistory(:, 2);
convergence.mae = metricsHistory(:, 3);

% Return all results in a single struct. The 'name' field is overwritten
% by run_comparison with the name from the JSON config, but we provide a
% default here for standalone usage.
result = struct(...
    'name', "MEWLS", ...
    'prediction', finalPrediction, ...       % N-by-1 fitted values
    'coefficients', coeffHistory{end}, ...   % M-by-1 polynomial coeffs
    'metrics', metricsStruct, ...            % struct with rmse/maxAbsError/mae
    'convergence', convergence, ...          % table: degree x metrics
    'weights', weights);                     % N-by-1 MEWLS weights (for diagnostics)

end


% =====================================================================
%  Local function: compute_entropy_weights
%  Implements the MEWLS weight function (Brugnano et al. 2024, Sec. 3).
% =====================================================================
function weights = compute_entropy_weights(xy, scale, isOutlier, outlierPenalty)
%COMPUTE_ENTROPY_WEIGHTS Build normalised entropy-based sample weights.
%   Implements: w_i = exp(-lambda * ||x_i - c||^2)
%   where c = centroid, lambda = scale.

% Step 1: Compute the centroid (mean) of the 2D point cloud.
% This is the reference centre for the distance-based weighting.
centroid = mean(xy, 1);  % 1-by-2 vector [mean_x, mean_y]

% Step 2: Squared Euclidean distance of each point from the centroid.
% d_i^2 = (x_i - c_x)^2 + (y_i - c_y)^2
% We keep the squared form because the weight function uses d^2 directly.
dist2 = sum((xy - centroid) .^ 2, 2);  % N-by-1

% Step 3: Entropy-inspired weight (Brugnano et al. 2024, Eq. in Sec. 3):
%   w_i = exp(-lambda * d_i^2)
% Points close to the centroid receive weight ~ 1; distant points receive
% exponentially smaller weights controlled by lambda (= scale).
raw = exp(-scale * dist2);  % N-by-1, all positive

% Step 4: Down-weight outlier samples by the penalty factor.
% Multiplying by outlierPenalty < 1 reduces their influence on the WLS fit.
raw(isOutlier) = raw(isOutlier) * outlierPenalty;

% Step 5: Normalise weights to sum to 1 (discrete probability distribution).
% This makes the scale of W independent of the number of samples N.
weights = raw / sum(raw);

% Step 6: Clamp any zero or near-zero weights to machine epsilon.
% Prevents singular or ill-conditioned (Phi' W Phi) in the normal equations.
weights(weights <= 0) = eps;
end
