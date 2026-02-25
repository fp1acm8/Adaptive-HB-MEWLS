function result = leastSquaresSolver(dataset, method_data)
%LEASTSQUARESSOLVER Polynomial surrogate using Ordinary Least Squares.
%   RESULT = LEASTSQUARESSOLVER(DATASET, METHOD_DATA) fits a polynomial
%   surface to the observed samples using ordinary least squares (OLS).
%   The solver sweeps polynomial degrees from 1 to METHOD_DATA.degree and
%   retains the full convergence history so that the user can inspect how
%   the error evolves with increasing polynomial complexity.
%
%   This solver is the unweighted baseline (W = I) for comparison against
%   the entropy-weighted MEWLS variant. In the notation of Brugnano et al.
%   (2024), setting all weights w_i = 1/N reduces the MEWLS formulation
%   to standard OLS. The comparison quantifies the improvement obtained by
%   the MEWLS weighting strategy.
%
%   Inputs:
%       DATASET     - struct produced by adaptivehb.io.load_dataset containing
%                     at least data (N-by-2), f (N-by-1), and f_true.
%       METHOD_DATA - struct with optional fields:
%                       degree : maximum polynomial degree (default 3).
%
%   Output:
%       RESULT - struct with the following fields:
%           name         - solver identifier string ("LeastSquares")
%           QI           - N-by-1 fitted values at input points (degree = method_data.degree)
%           QI_coeff     - polynomial coefficient vector for the final degree
%           metrics      - struct with rmse, maxAbsError, mae for the final degree
%           convergence  - table with columns degree, rmse, maxAbsError, mae
%
%   Nota: l'implementazione a B-spline gerarchiche adattive (HBS) che
%   estende questo solver polinomiale richiede la libreria GeoPDEs per MATLAB.
%   Disponibile su: https://rafavzqz.github.io/geopdes/
%   La funzione analoga nel codice HBS e':
%       getcoeff_weighted_least_squares_pen(hspace,hmsh,data,f,0,weight)
%
%   See also adaptivehb.solvers.mewlsSolver,
%            adaptivehb.solvers.polynomialDesignMatrix,
%            adaptivehb.solvers.compute_metrics.

% --- Input validation (MATLAB R2020b+ arguments block) ---------------
% dataset must be a scalar struct; method_data defaults to empty struct so
% the caller can omit it entirely.
arguments
    dataset (1, 1) struct
    method_data (1, 1) struct = struct()
end

% --- Default parameter handling ---------------------------------------

% degree (d): highest total polynomial degree d used in the sweep.
% For degree d the bivariate polynomial space Pi_d^2 has
% M = (d+1)(d+2)/2 basis functions. Default d=3 -> M=10.
if ~isfield(method_data, 'degree') || isempty(method_data.degree)
    method_data.degree = 3;
end

% Clamp to at least 1 (degree 0 = constant is trivial).
degree = max(1, method_data.degree);

% --- Extract data vectors for readability (stile vecchio codice) ------
data = dataset.data;    % N-by-2 coordinate matrix
f    = dataset.f;       % N-by-1 observed values
f_true = dataset.f_true; % N-by-1 true values (per il calcolo degli errori)

% --- Polynomial degree sweep ------------------------------------------
% Fit polynomials of total degree d = 1, 2, ..., degree and record
% error metrics at each step to produce a convergence curve.
metricsHistory = zeros(degree, 3); % each row: [RMSE, maxAbsErr, MAE]
degreeList = 1:degree;             % degrees to evaluate
QI_coeff_history = cell(1, degree);    % polynomial coefficients per degree
QI_history       = cell(1, degree);    % QI values per degree

for d = degreeList
    % Build the design matrix Phi_d in R^{N x M_d} for the bivariate
    % polynomial space Pi_d^2.  Each column is x^p * y^q with p+q <= d.
    % Ref: same Phi matrix used in MEWLS (Brugnano et al. 2024, Sec. 3),
    % but here solved without the weight matrix W.
    [A, ~] = adaptivehb.solvers.polynomialDesignMatrix(data, d);

    % Solve the OLS problem via MATLAB backslash operator:
    %   QI_coeff = Phi \ f
    % This minimises ||Phi*QI_coeff - f||_2^2 using QR factorisation.
    % Equivalent to the MEWLS normal equations with W = I (uniform weights).
    QI_coeff = A \ f;

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

% Pack the three metric values into a named struct for easy access.
metricNames = {"rmse", "maxAbsError", "mae"};
metricsStruct = cell2struct(num2cell(metricsHistory(end, :)), metricNames, 2);

% Build a convergence table (degree vs metrics) for plotting/analysis.
convergence = table();
convergence.degree = degreeList';       % polynomial degrees 1..degree
convergence.rmse = metricsHistory(:, 1);
convergence.maxAbsError = metricsHistory(:, 2);
convergence.mae = metricsHistory(:, 3);

% Return all results. Note: unlike mewlsSolver, no 'weight' field is
% included because OLS uses uniform (implicit) weights.
result = struct(...
    'name', "LeastSquares", ...
    'QI', QI_final, ...              % N-by-1 fitted values
    'QI_coeff', QI_coeff_history{end}, ... % M-by-1 polynomial coefficients
    'metrics', metricsStruct, ...    % struct with rmse/maxAbsError/mae
    'convergence', convergence);     % table: degree x metrics

end
