function dataset = load_dataset(filePath, options)
%LOAD_DATASET Load a point cloud dataset from disk.
%   DATASET = LOAD_DATASET(FILEPATH, OPTIONS) reads either a clean dataset
%   (three columns: x, y, value) or a dataset with noise/outliers (five
%   columns: x, y, f_true, f, is_outlier). The returned struct contains the
%   raw points, the normalized coordinates, and helper vectors used by the
%   comparison harness.
%
%   OPTIONS.normalize (logical) controls whether coordinates are mapped to
%   the unit square. When false, the original coordinates are preserved.
%
%   The dataset struct exposes the fields:
%       - originalXY  : coordinates as found on disk
%       - data        : coordinates (normalized if requested), N-by-2
%       - f_true      : true signal when available, otherwise observed
%       - f           : noisy observations used as solver input, N-by-1
%       - is_outlier  : logical mask for outliers (false when unavailable)
%
%   The coordinate layout is compatible with the [0,1]x[0,1] domain used
%   by the GeoPDEs library (geo_square.txt) in the hierarchical B-spline
%   component. See: https://rafavzqz.github.io/geopdes/
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
%   5 columns: x, y, f_true, f, is_outlier (noisy dataset)
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
        f_true = raw(:, 3);
        f = f_true;
        % No outliers in a clean dataset.
        is_outlier = false(size(f_true));

    case 5
        % Noisy dataset:
        %   col 3 = true (uncorrupted) values
        %   col 4 = observed (noisy) values used for fitting
        %   col 5 = outlier flag (0 or 1)
        f_true = raw(:, 3);
        f = raw(:, 4);

        % Validate that column 5 contains only binary 0/1 values before
        % converting to logical.  A non-binary value most likely means
        % the columns are in the wrong order or the file format is wrong.
        outlierCol = raw(:, 5);
        if ~all(ismember(outlierCol, [0, 1]))
            error('adaptivehb:io:load_dataset:invalidOutlierColumn', ...
                ['Column 5 of "%s" must contain only binary outlier flags ' ...
                 '(0 = inlier, 1 = outlier). ' ...
                 'Found non-binary values; check the column order ' ...
                 '[x, y, f_true, f, is_outlier].'], filePath);
        end

        is_outlier = logical(outlierCol);  % convert 0/1 to logical
end

% --- Optional coordinate normalisation --------------------------------
% Map coordinates to [0,1]^2 using min-max scaling per axis.
% This is critical for polynomial fitting: without normalisation, terms
% like x^3 can vary over many orders of magnitude.
data = originalXY;
if options.normalize
    data = normalize_coordinates(originalXY);
end

% --- Handle legacy datasets without a true-value column ---------------
% Some older datasets may have NaN or empty values in the true column.
% In that case, fall back to observed values so that downstream metric
% computations (f_true - QI) remain well-defined.
if isempty(f_true) || all(isnan(f_true))
    f_true = f;
end

% --- Pack everything into a single output struct ----------------------
% All fields are forced to column vectors via (:) for consistency.
dataset = struct(...
    'filePath', filePath, ...       % path to the source file (for traceability)
    'originalXY', originalXY, ...   % raw coordinates as read from disk
    'data', data, ...               % normalised coordinates (or raw if normalize=false), N-by-2
    'f_true', f_true(:), ...        % true signal (N-by-1)
    'f', f(:), ...                  % observed signal, possibly noisy (N-by-1)
    'is_outlier', is_outlier(:));   % logical outlier mask (N-by-1)

end


% =====================================================================
%  Local function: normalize_coordinates
%  Maps each coordinate axis independently to [0, 1].
% =====================================================================
function dataNorm = normalize_coordinates(data)
% Compute per-axis minimum and maximum.
minVals = min(data, [], 1);   % 1-by-2: [min_x, min_y]
maxVals = max(data, [], 1);   % 1-by-2: [max_x, max_y]

% Range of each axis. Guard against zero range (constant coordinate)
% by replacing zeros with 1, which maps constant axes to 0.
rangeVals = maxVals - minVals;
rangeVals(rangeVals == 0) = 1;

% Apply min-max scaling: data_norm = (data - min) / range.
dataNorm = (data - minVals) ./ rangeVals;
end
