function dataset = load_dataset(filePath, options)
%LOAD_DATASET Load a point cloud dataset from disk.
%   DATASET = LOAD_DATASET(FILEPATH, OPTIONS) reads either a clean dataset
%   (three columns: x, y, value) or a dataset with noise/outliers (five
%   columns: x, y, trueValue, observedValue, isOutlier). The returned
%   struct contains the raw points, the normalized coordinates, and helper
%   vectors used by the comparison harness.
%
%   OPTIONS.normalize (logical) controls whether coordinates are mapped to
%   the unit square. When false, the original coordinates are preserved.
%
%   The dataset struct exposes the fields:
%       - originalXY : coordinates as found on disk
%       - xy         : coordinates (normalized if requested)
%       - fTrue      : true signal when available, otherwise observed
%       - fObserved  : noisy observations used as solver input
%       - isOutlier  : logical mask for outliers (false when unavailable)
%
%   This helper is intentionally self-contained to keep the orchestrator
%   independent from the original demo scripts.
%
%   See also adaptivehb.io.apply_noise.

arguments
    filePath (1, :) char
    options (1, 1) struct = struct()
end

if ~isfield(options, 'normalize') || isempty(options.normalize)
    options.normalize = true;
end

if ~isfile(filePath)
    error('adaptivehb:io:FileNotFound', ...
        'Dataset file not found: %s', filePath);
end

raw = readmatrix(filePath);
if isempty(raw)
    error('adaptivehb:io:EmptyDataset', ...
        'Dataset %s does not contain any samples.', filePath);
end

nCols = size(raw, 2);
if ~ismember(nCols, [3, 5])
    error('adaptivehb:io:UnsupportedFormat', ...
        'Dataset %s must have either 3 or 5 columns. Found %d.', ...
        filePath, nCols);
end

originalXY = raw(:, 1:2);

switch nCols
    case 3
        fTrue = raw(:, 3);
        fObserved = fTrue;
        isOutlier = false(size(fTrue));
    case 5
        fTrue = raw(:, 3);
        fObserved = raw(:, 4);
        isOutlier = logical(raw(:, 5));
end

xy = originalXY;
if options.normalize
    xy = normalize_coordinates(originalXY);
end

if isempty(fTrue) || all(isnan(fTrue))
    % In legacy datasets the true field might not be recorded. Fall back to
    % the observed values so downstream metrics remain defined.
    fTrue = fObserved;
end

dataset = struct(...
    'filePath', filePath, ...
    'originalXY', originalXY, ...
    'xy', xy, ...
    'fTrue', fTrue(:), ...
    'fObserved', fObserved(:), ...
    'isOutlier', isOutlier(:));

end

function xyNorm = normalize_coordinates(xy)
minVals = min(xy, [], 1);
maxVals = max(xy, [], 1);
rangeVals = maxVals - minVals;
rangeVals(rangeVals == 0) = 1;
xyNorm = (xy - minVals) ./ rangeVals;
end
