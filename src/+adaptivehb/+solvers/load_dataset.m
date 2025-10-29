function dataset = load_dataset(dataset_config)
%LOAD_DATASET Load and normalise scattered data from a text file.
%   DATASET = LOAD_DATASET(DATASET_CONFIG) reads the dataset specified by
%   DATASET_CONFIG.full_path. The helper returns a struct containing the
%   normalised coordinates, measurements and bookkeeping information. When
%   DATASET_CONFIG.has_outliers is provided, the loader enforces the
%   expected column layout. Otherwise it attempts to detect whether the
%   dataset has three (clean) or five (noisy) columns.
%
%   The returned struct includes:
%       data        - normalised coordinates in [0, 1]^2
%       data0       - original coordinates before normalisation
%       f_true      - ground truth values
%       f_noise     - noisy observations used for fitting
%       is_outlier  - logical mask (all false when dataset has 3 columns)
%       expected_cols - number of columns detected
%       filename    - absolute path to the dataset file
%
%   See also FULLFILE, FSCANF.

arguments
    dataset_config (1, 1) struct
end

if ~isfield(dataset_config, 'full_path')
    error('adaptivehb:solvers:MissingDatasetPath', ...
        'dataset_config.full_path is required');
end

file_path = dataset_config.full_path;
if ~isfile(file_path)
    error('adaptivehb:solvers:DatasetNotFound', ...
        'Dataset file not found: %s', file_path);
end

fid = fopen(file_path, 'r');
cleaner = onCleanup(@() fclose(fid));
raw_vector = fscanf(fid, '%f');

if isempty(raw_vector)
    error('adaptivehb:solvers:EmptyDataset', ...
        'Dataset file %s does not contain numeric data.', file_path);
end

total_elements = numel(raw_vector);

if isfield(dataset_config, 'has_outliers')
    has_outliers = logical(dataset_config.has_outliers);
    expected_cols = ternary(has_outliers, 5, 3);
else
    expected_cols = detect_column_count(total_elements);
end

if mod(total_elements, expected_cols) ~= 0
    error('adaptivehb:solvers:InvalidDatasetFormat', ...
        ['Total number of elements (%d) is not divisible by the detected ' ...
         'column count (%d).'], total_elements, expected_cols);
end

M = total_elements / expected_cols;
raw = reshape(raw_vector, expected_cols, M);

data0 = raw(1:2, :)';
[f_true, f_noise, is_outlier] = parse_measurements(raw, expected_cols);

a1 = min(data0(:, 1));
b1 = max(data0(:, 1));
a2 = min(data0(:, 2));
b2 = max(data0(:, 2));

if (b1 - a1) <= 0 || (b2 - a2) <= 0
    error('adaptivehb:solvers:DegenerateDomain', ...
        'Dataset coordinates must span a non-zero range in both axes.');
end

data = [(data0(:, 1) - a1) / (b1 - a1), (data0(:, 2) - a2) / (b2 - a2)];

dataset = struct(
    'data', data,
    'data0', data0,
    'f_true', f_true,
    'f_noise', f_noise,
    'is_outlier', is_outlier,
    'expected_cols', expected_cols,
    'filename', file_path,
    'relative_path', get_field_or(dataset_config, 'relative_path', file_path),
    'num_points', M,
    'domain', struct('x', [a1, b1], 'y', [a2, b2]));
end

function expected_cols = detect_column_count(total_elements)
    candidates = [5, 3];
    divisible = candidates(mod(total_elements, candidates) == 0);
    if isempty(divisible)
        error('adaptivehb:solvers:UnknownDatasetFormat', ...
            ['Unable to infer the dataset layout. The total number of ' ...
             'elements (%d) is not compatible with 3- or 5-column formats.'], ...
             total_elements);
    end
    % Prefer the format with more columns when ambiguous (e.g., 15 elements)
    expected_cols = max(divisible);
end

function [f_true, f_noise, is_outlier] = parse_measurements(raw, expected_cols)
    switch expected_cols
        case 3
            f_true = raw(3, :)';
            f_noise = f_true;
            is_outlier = false(size(f_true));
        case 5
            f_true = raw(3, :)';
            f_noise = raw(4, :)';
            is_outlier = logical(raw(5, :)');
        otherwise
            error('adaptivehb:solvers:UnsupportedColumnCount', ...
                'Only 3- and 5-column datasets are supported.');
    end
end

function value = ternary(condition, when_true, when_false)
    if condition
        value = when_true;
    else
        value = when_false;
    end
end

function value = get_field_or(structure, field, default)
    if isfield(structure, field)
        value = structure.(field);
    else
        value = default;
    end
end
