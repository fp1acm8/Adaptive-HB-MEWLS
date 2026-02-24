function metrics = compute_metrics(fTrue, prediction)
%COMPUTE_METRICS Compute standard approximation error metrics.
%   METRICS = COMPUTE_METRICS(FTRUE, PREDICTION) returns a 1-by-3 vector
%   [RMSE, MAXABSERROR, MAE] that summarises the discrepancy between the
%   true signal FTRUE and the solver output PREDICTION.
%
%   The three metrics are:
%       RMSE         - Root Mean Squared Error: sqrt(mean((fTrue - pred).^2))
%       maxAbsError  - Maximum Absolute Error : max(|fTrue - pred|)
%       MAE          - Mean Absolute Error    : mean(|fTrue - pred|)
%
%   Both inputs must be column vectors of equal length.
%
%   Example:
%       m = adaptivehb.solvers.compute_metrics(fTrue, prediction);
%       fprintf('RMSE = %.4f, MaxAbs = %.4f, MAE = %.4f\n', m(1), m(2), m(3));
%
%   See also adaptivehb.solvers.leastSquaresSolver,
%            adaptivehb.solvers.mewlsSolver.

% Point-wise error vector: e_i = f_true_i - f_hat_i.
% Positive means under-prediction, negative means over-prediction.
residual = fTrue - prediction;

% RMSE = sqrt( 1/N * sum(e_i^2) ).
% Penalises large errors more than small ones due to the squaring.
% This is the L2-norm of the residual divided by sqrt(N).
rmse = sqrt(mean(residual .^ 2));

% Maximum Absolute Error = ||e||_inf = max_i |e_i|.
% Reports the single worst-case point error across the entire dataset.
maxAbs = max(abs(residual));

% MAE = 1/N * sum(|e_i|).
% Average absolute deviation — less sensitive to outliers than RMSE.
mae = mean(abs(residual));

% Return as a 1-by-3 row vector in the order [RMSE, maxAbsError, MAE].
% Both leastSquaresSolver and mewlsSolver unpack this vector into the
% metricsHistory matrix and subsequently into a named struct.
metrics = [rmse, maxAbs, mae];

end
