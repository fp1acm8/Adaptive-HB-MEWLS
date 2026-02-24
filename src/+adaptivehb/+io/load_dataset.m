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

% --- Input validation -------------------------------------------------
% filePath must be a character vector; options defaults to empty struct.
arguments
    filePath (1, :) char
    options (1, 1) struct = struct()
end

% Default: normalise coordinates to [0,1]^2 unless explicitly disabled.
% Normalisation is important for polynomial fitting because it prevents
% numerical instability caused by large coordinate values raised to high
% powers in the design matrix Phi.
if ~isfield(options, 'normalize') || isempty(options.normalize)
    options.normalize = true;
end

% --- File existence check ---------------------------------------------
if ~isfile(filePath)
    error('adaptivehb:io:FileNotFound', ...
        'Dataset file not found: %s', filePath);
end

% --- Read the file into a numeric matrix ------------------------------
% readmatrix auto-detects the delimiter (space, tab, comma).
raw = readmatrix(filePath);

% Guard against empty files.
if isempty(raw)
    error('adaptivehb:io:EmptyDataset', ...
        'Dataset %s does not contain any samples.', filePath);
end

% --- Validate the number of columns ----------------------------------
% Only two formats are supported:
%   3 columns: x, y, f              (clean dataset)
%   5 columns: x, y, fTrue, fObs, isOutlier (noisy dataset)
nCols = size(raw, 2);
if ~ismember(nCols, [3, 5])
    error('adaptivehb:io:UnsupportedFormat', ...
        'Dataset %s must have either 3 or 5 columns. Found %d.', ...
        filePath, nCols);
end

% --- Parse columns based on format -----------------------------------
% Extract the 2D coordinates (always in columns 1-2).
originalXY = raw(:, 1:2);

switch nCols
    case 3
        % Clean dataset: column 3 is the only value; use it as both
        % the true signal and the observed signal (no noise present).
        fTrue = raw(:, 3);
        fObserved = fTrue;
        % No outliers in a clean dataset.
        isOutlier = false(size(fTrue));

    case 5
        % Noisy dataset:
        %   col 3 = true (uncorrupted) values
        %   col 4 = observed (noisy) values used for fitting
        %   col 5 = outlier flag (0 or 1)
        fTrue = raw(:, 3);
        fObserved = raw(:, 4);
        isOutlier = logical(raw(:, 5));  % convert 0/1 to logical
end

% --- Optional coordinate normalisation --------------------------------
% Map coordinates to [0,1]^2 using min-max scaling per axis.
% This is critical for polynomial fitting: without normalisation, terms
% like x^3 can vary over many orders of magnitude.
xy = originalXY;
if options.normalize
    xy = normalize_coordinates(originalXY);
end

% --- Handle legacy datasets without a true-value column ---------------
% Some older datasets may have NaN or empty values in the true column.
% In that case, fall back to observed values so that downstream metric
% computations (fTrue - prediction) remain well-defined.
if isempty(fTrue) || all(isnan(fTrue))
    fTrue = fObserved;
end

% --- Pack everything into a single output struct ----------------------
% All fields are forced to column vectors via (:) for consistency.
dataset = struct(...
    'filePath', filePath, ...       % path to the source file (for traceability)
    'originalXY', originalXY, ...   % raw coordinates as read from disk
    'xy', xy, ...                   % normalised coordinates (or raw if normalize=false)
    'fTrue', fTrue(:), ...          % true signal (N-by-1)
    'fObserved', fObserved(:), ...  % observed signal, possibly noisy (N-by-1)
    'isOutlier', isOutlier(:));     % logical outlier mask (N-by-1)

end


% =====================================================================
%  Local function: normalize_coordinates
%  Maps each coordinate axis independently to [0, 1].
% =====================================================================
function xyNorm = normalize_coordinates(xy)
% Compute per-axis minimum and maximum.
minVals = min(xy, [], 1);   % 1-by-2: [min_x, min_y]
maxVals = max(xy, [], 1);   % 1-by-2: [max_x, max_y]

% Range of each axis. Guard against zero range (constant coordinate)
% by replacing zeros with 1, which maps constant axes to 0.
rangeVals = maxVals - minVals;
rangeVals(rangeVals == 0) = 1;

% Apply min-max scaling: xy_norm = (xy - min) / range.
xyNorm = (xy - minVals) ./ rangeVals;
end
