function summary = summarize_errors(predictions, truth, isOutlier)
%SUMMARIZE_ERRORS Compute scalar error metrics and table information.
%   SUMMARY = SUMMARIZE_ERRORS(PREDICTIONS, TRUTH, ISOUTLIER) returns a
%   struct with pointwise absolute errors, global metrics and metadata that
%   can be used to assemble a textual table. ISOUTLIER is optional; if
%   omitted or empty the dataset is treated as outlier-free.
%
%   Fields in SUMMARY:
%       pointwise       Absolute error at each sample
%       metrics         Struct with fields maxAll, rmseAll, maxClean,
%                       rmseClean (NaN if not applicable)
%       headers         Cell array of column headers for table printing
%       values          Numeric row matching HEADERS
%       formatSpec      sprintf format string for VALUES
%
%   See also adaptivehb.core.prepare_dataset, adaptivehb.core.run_adaptivity.

if nargin < 3 || isempty(isOutlier)
    isOutlier = false(size(predictions));
else
    isOutlier = logical(isOutlier);
end

predictions = predictions(:);
truth = truth(:);

if numel(predictions) ~= numel(truth)
    error('summarize_errors:DimensionMismatch', ...
        'PREDICTIONS and TRUTH must have the same number of entries.');
end

if numel(isOutlier) ~= numel(truth)
    error('summarize_errors:MaskMismatch', ...
        'ISOUTLIER must match the number of samples.');
end

errors = abs(truth - predictions);
maxAll = max(errors);
rmseAll = sqrt(mean(errors.^2));

cleanMask = ~isOutlier;
if any(cleanMask)
    cleanErrors = errors(cleanMask);
    maxClean = max(cleanErrors);
    rmseClean = sqrt(mean(cleanErrors.^2));
else
    maxClean = NaN;
    rmseClean = NaN;
end

hasOutliers = any(isOutlier);

if hasOutliers
    headers = {'max_err_all', 'rmse_all', 'max_err_clean', 'rmse_clean'};
    values = [maxAll, rmseAll, maxClean, rmseClean];
    formatSpec = '%14.4e %14.4e %14.4e %14.4e';
else
    headers = {'max_err', 'rmse'};
    values = [maxAll, rmseAll];
    formatSpec = '%10.4e %10.4e';
end

metrics = struct( ...
    'maxAll', maxAll, ...
    'rmseAll', rmseAll, ...
    'maxClean', maxClean, ...
    'rmseClean', rmseClean, ...
    'hasOutliers', hasOutliers);

summary = struct( ...
    'pointwise', errors, ...
    'metrics', metrics, ...
    'headers', {headers}, ...
    'values', values, ...
    'formatSpec', formatSpec);
end
