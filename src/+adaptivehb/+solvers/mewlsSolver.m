function result = mewlsSolver(dataset, method_data)
%MEWLSSOLVER Maximum Entropy Weighted Least Squares polynomial fitter.
%   RESULT = MEWLSSOLVER(DATASET, METHOD_DATA) fits a polynomial surface
%   using Weighted Least Squares (WLS) where sample weights are derived
%   from a maximum-entropy inspired criterion.  Like the OLS baseline
%   (leastSquaresSolver), the solver sweeps polynomial degrees from 1 to
%   METHOD_DATA.degree and records the full convergence history.
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
%      by alpha_out (default 0.5), reducing their influence.
%   4. Normalise: weight = weight / sum(weight), then clamp to eps.
%   5. Build the design matrix Phi for polynomial space Pi_d^2.
%   6. Solve the WLS normal equations:
%          (Phi' W Phi) QI_coeff = Phi' W f          [Eq. in Sec. 3]
%
%   Inputs:
%       DATASET     - struct from adaptivehb.io.load_dataset with fields
%                     data (N-by-2), f (N-by-1), f_true, is_outlier.
%       METHOD_DATA - struct with optional fields:
%                       degree    : max polynomial degree (default 3)
%                       lambda    : controls weight decay rate (default 10).
%                                   Range: 1..50.
%                       alpha_out : multiplicative factor applied to
%                                   outlier weights (default 0.5).
%
%   Output:
%       RESULT - struct with fields:
%           name         - "MEWLS"
%           QI           - N-by-1 fitted values (final degree)
%           QI_coeff     - polynomial coefficients (final degree)
%           metrics      - struct with rmse, maxAbsError, mae
%           convergence  - table with degree-by-degree error history
%           weight       - N-by-1 normalised sample weights
%
%   Note: the adaptive hierarchical B-spline (HBS) implementation that
%   extends this polynomial solver requires the GeoPDEs library for MATLAB.
%   Available at: https://rafavzqz.github.io/geopdes/
%   The equivalent function in the HBS code is:
%       getcoeff_weighted_least_squares_pen(hspace,hmsh,data,f,lambda,weight)
%
%   See also adaptivehb.solvers.leastSquaresSolver,
%            adaptivehb.solvers.polynomialDesignMatrix,
%            adaptivehb.solvers.compute_metrics.

% --- Input validation (MATLAB R2020b+ arguments block) ---------------
% dataset must be a scalar struct; method_data defaults to empty struct so
% all parameters below become optional.
arguments
    dataset (1, 1) struct
    method_data (1, 1) struct = struct()
end

% --- Default parameter handling ---------------------------------------

% degree (d): highest total polynomial degree used in the sweep over
% the bivariate polynomial space Pi_d^2 (Brugnano et al. 2024, Sec. 3).
% Default d=3 gives M=(3+1)(3+2)/2 = 10 basis functions.
if ~isfield(method_data, 'degree') || isempty(method_data.degree)
    method_data.degree = 3;
end

% lambda: controls how fast the MEWLS weights decay with distance from
% the centroid.
%   lambda ~ 1-5   : mild decay, weights nearly uniform (close to OLS)
%   lambda ~ 10-20 : moderate decay (recommended starting range)
%   lambda > 30    : aggressive decay, strong centroid bias
% Ref: Brugnano et al. 2024, weight function w_i = exp(-lambda * d_i^2).
if ~isfield(method_data, 'lambda') || isempty(method_data.lambda)
    method_data.lambda = 10.0;
end

% alpha_out: multiplicative factor for outlier weights.
%   0.0 -> outliers completely ignored
%   0.5 -> outliers contribute at half weight (default)
%   1.0 -> no penalty (outlier flag has no effect)
if ~isfield(method_data, 'alpha_out') || isempty(method_data.alpha_out)
    method_data.alpha_out = 0.5;
end

% Clamp degree to at least 1 (degree 0 = constant is trivial).
degree = max(1, method_data.degree);

% Store lambda and penalty in local variables.
lambda    = method_data.lambda;    % weight decay scale parameter
alpha_out = method_data.alpha_out; % outlier weight penalty factor

% --- Extract data vectors for readability --------------------------------
data       = dataset.data;       % N-by-2 coordinate matrix
f          = dataset.f;          % N-by-1 observed values
f_true     = dataset.f_true;     % N-by-1 true values
is_outlier = dataset.is_outlier; % N-by-1 logical outlier mask

% --- Degenerate input guard -------------------------------------------
% If every sample is an outlier the weight vector collapses and the
% normal equations become meaningless.  Warn the caller rather than
% silently returning garbage.
if all(is_outlier)
    warning('adaptivehb:mewls:allOutliers', ...
        ['All %d data points are flagged as outliers. ' ...
         'MEWLS results may be unreliable. ' ...
         'Consider revising the outlier detection threshold.'], numel(is_outlier));
end

% --- Entropy weight computation (Brugnano et al. 2024, Sec. 3) -------
% Build the MEWLS weight vector (N-by-1):
%   w_i = exp(-lambda * ||x_i - centroid||^2)
% then penalise outliers and normalise so sum(weight) = 1.
weight = compute_entropy_weights(data, lambda, is_outlier, alpha_out);

% Construct the diagonal weight matrix W = diag(w_1, ..., w_N).
% This appears as W in the WLS normal equations: (Phi' W Phi) QI_coeff = Phi' W f.
W = diag(weight);

% --- Polynomial degree sweep ------------------------------------------
% Fit polynomials of total degree d = 1, 2, ..., degree and record
% error metrics at each step. This provides a convergence curve showing
% how the approximation improves with richer polynomial spaces Pi_d^2.
metricsHistory   = zeros(degree, 3); % each row: [RMSE, maxAbsErr, MAE]
degreeList       = 1:degree;         % degrees to evaluate
QI_coeff_history = cell(1, degree);  % polynomial coefficients per degree
QI_history       = cell(1, degree);  % QI values per degree

for d = degreeList
    % Build the design matrix Phi_d in R^{N x M_d} for the bivariate
    % polynomial space Pi_d^2 where M_d = (d+1)(d+2)/2 columns.
    % Each column is one monomial x^p * y^q with p+q <= d.
    % Ref: Brugnano et al. 2024, polynomial basis construction.
    [A, ~] = adaptivehb.solvers.polynomialDesignMatrix(data, d);

    % Solve the WLS normal equations (Brugnano et al. 2024, Sec. 3):
    %   QI_coeff = (Phi' W Phi) \ (Phi' W f)
    % This minimises the weighted residual:
    %   sum_i w_i * (f_i - p(x_i,y_i))^2
    % Analogo a: S = col_matrix*diag(weight)*col_matrix'; QI_coeff = S\rhs
    AtWA = A' * W * A;

    % Check the condition number before solving to detect ill-conditioned
    % normal equations.  High lambda or near-degenerate data can cause
    % A'WA to become singular.  condest is O(n^2) — cheaper than cond().
    condNum = condest(AtWA);
    if condNum > 1e10
        warning('adaptivehb:mewls:illConditioned', ...
            ['Normal equations are ill-conditioned (cond ~ %.2e) at degree %d. ' ...
             'Consider reducing lambda (currently %.1f) or normalizing coordinates.'], ...
            condNum, d, lambda);
    end

    QI_coeff = AtWA \ (A' * W * f);

    % Evaluate the fitted polynomial at all data points: QI = Phi * QI_coeff.
    QI = A * QI_coeff;

    % Compute RMSE, maxAbsError, MAE between true values and QI.
    metricsHistory(d, :) = adaptivehb.solvers.compute_metrics(f_true, QI);

    % Store coefficient vector and QI for this degree.
    QI_coeff_history{d} = QI_coeff;
    QI_history{d} = QI;
end

% --- Assemble output struct -------------------------------------------

% Use results from the highest polynomial degree as the final output.
QI_final = QI_history{end};

% Pack the three metric values from the last degree into a named struct.
metricNames = {"rmse", "maxAbsError", "mae"};
metricsStruct = cell2struct(num2cell(metricsHistory(end, :)), metricNames, 2);

% Build a convergence table so the caller can plot RMSE vs degree.
convergence = table();
convergence.degree = degreeList';       % polynomial degrees 1..degree
convergence.rmse = metricsHistory(:, 1);
convergence.maxAbsError = metricsHistory(:, 2);
convergence.mae = metricsHistory(:, 3);

% Return all results in a single struct.
result = struct(...
    'name', "MEWLS", ...
    'QI', QI_final, ...                     % N-by-1 fitted values
    'QI_coeff', QI_coeff_history{end}, ...  % M-by-1 polynomial coefficients
    'metrics', metricsStruct, ...           % struct with rmse/maxAbsError/mae
    'convergence', convergence, ...         % table: degree x metrics
    'weight', weight);                      % N-by-1 MEWLS weights (per diagnostica)

end


% =====================================================================
%  Local function: compute_entropy_weights
%  Implements the MEWLS weight function (Brugnano et al. 2024, Sec. 3).
% =====================================================================
function weight = compute_entropy_weights(data, lambda, is_outlier, alpha_out)
%COMPUTE_ENTROPY_WEIGHTS Build normalised entropy-based sample weights.
%   Implements: w_i = exp(-lambda * ||x_i - c||^2)
%   where c = centroid, lambda = scale parameter.

% Step 1: Compute the centroid (mean) of the 2D point cloud.
% This is the reference centre for the distance-based weighting.
centroid = mean(data, 1);  % 1-by-2 vector [mean_x, mean_y]

% Step 2: Squared Euclidean distance of each point from the centroid.
% d_i^2 = (x_i - c_x)^2 + (y_i - c_y)^2
% We keep the squared form because the weight function uses d^2 directly.
dist2 = sum((data - centroid) .^ 2, 2);  % N-by-1

% Step 3: Entropy-inspired weight (Brugnano et al. 2024, Eq. in Sec. 3):
%   w_i = exp(-lambda * d_i^2)
% Points close to the centroid receive weight ~ 1; distant points receive
% exponentially smaller weights controlled by lambda.
raw = exp(-lambda * dist2);  % N-by-1, all positive

% Step 4: Down-weight outlier samples by the penalty factor alpha_out.
% Multiplying by alpha_out < 1 reduces their influence on the WLS fit.
raw(is_outlier) = raw(is_outlier) * alpha_out;

% Step 5: Normalise weights to sum to 1 (discrete probability distribution).
% This makes the scale of W independent of the number of samples N.
weight = raw / sum(raw);

% Step 6: Clamp any zero or near-zero weights to machine epsilon.
% Prevents singular or ill-conditioned (Phi' W Phi) in the normal equations.
weight(weight <= 0) = eps;
end
