function metrics = compute_metrics(f_true, QI)
%COMPUTE_METRICS Compute standard approximation error metrics.
%   METRICS = COMPUTE_METRICS(F_TRUE, QI) returns a 1-by-3 vector
%   [RMSE, MAXABSERROR, MAE] that summarises the discrepancy between the
%   true signal F_TRUE and the solver output QI.
%
%   The three metrics are:
%       RMSE         - Root Mean Squared Error: sqrt(mean((f_true - QI).^2))
%       maxAbsError  - Maximum Absolute Error : max(|f_true - QI|)
%       MAE          - Mean Absolute Error    : mean(|f_true - QI|)
%
%   Both inputs must be column vectors of equal length.
%
%   Example:
%       m = adaptivehb.solvers.compute_metrics(f_true, QI);
%       fprintf('RMSE = %.4f, MaxAbs = %.4f, MAE = %.4f\n', m(1), m(2), m(3));
%
%   See also adaptivehb.solvers.leastSquaresSolver,
%            adaptivehb.solvers.mewlsSolver.

% Point-wise error vector: e_i = f_true_i - QI_i.
% Positive means under-prediction, negative means over-prediction.
error = f_true - QI;

% RMSE = sqrt( 1/N * sum(e_i^2) ).
% Penalises large errors more than small ones due to the squaring.
% This is the L2-norm of the error vector divided by sqrt(N).
rmse = sqrt(mean(error .^ 2));

% Maximum Absolute Error = ||e||_inf = max_i |e_i|.
% Reports the single worst-case point error across the entire dataset.
maxAbs = max(abs(error));

% MAE = 1/N * sum(|e_i|).
% Average absolute deviation — less sensitive to outliers than RMSE.
mae = mean(abs(error));

% Return as a 1-by-3 row vector in the order [RMSE, maxAbsError, MAE].
% Both leastSquaresSolver and mewlsSolver unpack this vector into the
% metricsHistory matrix and subsequently into a named struct.
metrics = [rmse, maxAbs, mae];

end
