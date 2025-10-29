function dataset = prepare_dataset(config)
%PREPARE_DATASET Load and normalize scattered data for adaptive solvers.
%   DATASET = PREPARE_DATASET(CONFIG) reads the dataset described by CONFIG
%   and returns a struct with normalized coordinates, raw coordinates and
%   metadata about outliers. CONFIG must contain at least the field
%   `filename`. Optional fields:
%       - outlierFlag (logical or numeric 0/1)
%       - verbose (logical) to enable command window output.
%
%   The returned DATASET struct contains:
%       points          Normalized coordinates in [0,1]^d
%       rawPoints       Original coordinates
%       fTrue           Ground truth values
%       fObserved       Values to be used for fitting (noisy measurements)
%       isOutlier       Logical mask of outliers (empty if not provided)
%       bounds          Struct with min/max for each coordinate
%       info            Metadata (filename, number of points, etc.)
%
%   This helper encapsulates the dataset preparation logic that previously
%   lived inside the solver scripts.
%
%   See also adaptivehb.core.run_adaptivity, adaptivehb.core.summarize_errors.

if nargin == 0 || ~isstruct(config)
    error('prepare_dataset:InvalidInput', 'CONFIG must be a struct.');
end

if ~isfield(config, 'filename')
    error('prepare_dataset:MissingField', 'CONFIG must contain a filename field.');
end

filename = config.filename;
if ~(isstring(filename) || ischar(filename))
    error('prepare_dataset:InvalidFilename', 'CONFIG.filename must be a string or char vector.');
end
filename = char(filename);

if isfield(config, 'outlierFlag')
    outlierFlag = logical(config.outlierFlag);
else
    outlierFlag = false;
end

if isfield(config, 'verbose')
    verbose = logical(config.verbose);
else
    verbose = false;
end

fid = fopen(filename, 'r');
if fid == -1
    error('prepare_dataset:FileNotFound', 'Unable to open dataset file ''%s''.', filename);
end
raw_vector = fscanf(fid, '%f');
fclose(fid);

total_elements = numel(raw_vector);

if outlierFlag
    expected_cols = 5;
else
    expected_cols = 3;
end

if mod(total_elements, expected_cols) ~= 0
    error('prepare_dataset:InvalidFormat', ...
        'Dataset format mismatch: %d elements cannot be reshaped into %d columns.', ...
        total_elements, expected_cols);
end

num_points = total_elements / expected_cols;
raw_matrix = reshape(raw_vector, expected_cols, num_points);
rawPoints = raw_matrix(1:2, :)';

if expected_cols == 3
    fTrue = raw_matrix(3, :)';
    fObserved = fTrue;
    isOutlier = false(size(fTrue));
else
    fTrue = raw_matrix(3, :)';
    fObserved = raw_matrix(4, :)';
    isOutlier = logical(raw_matrix(5, :)');
end

mins = min(rawPoints, [], 1);
maxs = max(rawPoints, [], 1);
span = maxs - mins;
span(span == 0) = 1;  % Avoid division by zero if coordinate is constant

points = (rawPoints - mins) ./ span;

bounds.x = [mins(1), maxs(1)];
bounds.y = [mins(2), maxs(2)];

if verbose
    fprintf('Dataset loaded: %s\n', filename);
    fprintf('  Points: %d | Columns: %d | Outliers: %d (%.2f%%%%)\n', ...
        num_points, expected_cols, nnz(isOutlier), 100 * nnz(isOutlier) / num_points);
end

info = struct( ...
    'filename', filename, ...
    'numPoints', num_points, ...
    'hasOutliers', any(isOutlier), ...
    'expectedColumns', expected_cols);

dataset = struct( ...
    'points', points, ...
    'rawPoints', rawPoints, ...
    'fTrue', fTrue, ...
    'fObserved', fObserved, ...
    'isOutlier', isOutlier, ...
    'bounds', bounds, ...
    'info', info);
end
