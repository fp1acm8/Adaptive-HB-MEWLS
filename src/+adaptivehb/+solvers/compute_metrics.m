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

residual = fTrue - prediction;
rmse     = sqrt(mean(residual .^ 2));
maxAbs   = max(abs(residual));
mae      = mean(abs(residual));
metrics  = [rmse, maxAbs, mae];

end
