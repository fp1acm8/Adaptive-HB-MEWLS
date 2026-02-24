%% generate_noisy_dataset.m
% Interactive wrapper for adding synthetic noise to an existing dataset
% using adaptivehb.data.apply_noise. The script walks the user through
% noise configuration, applies the perturbation with full reproducibility,
% and saves the augmented data to a file whose name encodes the noise
% parameters (type, fraction, intensity, seed).
%
% Usage:
%   1. Run the script from the repository root:
%        >> cd /path/to/Adaptive-HB-MEWLS
%        >> addpath('scripts');
%        >> generate_noisy_dataset
%   2. Follow the on-screen prompts to select the input dataset,
%      output directory, noise type, and corruption parameters.
%
% The output file has five columns:
%   x  y  fTrue  fNoisy  isOutlier
%
% See also adaptivehb.data.apply_noise, adaptivehb.data.format_noise_filename.

% Clear the workspace, close all figures, and clear the command window
% to start with a clean environment.
clearvars; close all; clc;

% Use long format for numeric display during the interactive session.
format long;

% =====================================================================
%  Bootstrap: add the src/ package folder to the MATLAB path
% =====================================================================

% Resolve the absolute path of this script file.
thisFile = mfilename('fullpath');

% Go one level up from scripts/ to reach the repository root, then
% append 'src' to locate the +adaptivehb package.
projectRoot = fileparts(thisFile);
srcPath = fullfile(projectRoot, '..', 'src');

% Abort with a clear message if the src folder cannot be found.
if ~isfolder(srcPath)
    error('Cannot find the src folder: %s', srcPath);
end

% Make the +adaptivehb package visible in this MATLAB session.
addpath(srcPath);

%% ================== PARAMETERS ==================
% The following section collects all user inputs interactively.
% Each input() call reads a string from the command window.

% --- Input dataset path ---
% The user must provide a path to a text file with exactly 3 columns: x y f.
input_path = input('Path to the input dataset (format: x y f): ', 's');
if isempty(input_path)
    error('The input dataset path is required.');
end

% Verify that the file exists on disk before attempting to load it.
if ~isfile(input_path)
    error('File not found: %s', input_path);
end

% Load the dataset as a numeric matrix using MATLAB's built-in load().
% Expected shape: N rows x 3 columns (x, y, f).
data = load(input_path);
if size(data, 2) ~= 3
    error('The file must contain exactly 3 columns: x y f');
end

% --- Output directory ---
% If left empty, defaults to the same directory as the input file.
output_dir = input('Output directory: ', 's');
if isempty(output_dir)
    output_dir = fileparts(input_path);
end

% Validate that the output directory exists.
if ~isfolder(output_dir)
    error('Invalid output directory: %s', output_dir);
end

% --- Output filename prefix ---
% An optional prefix (e.g. "outl") prepended to the generated filename.
% Example: prefix "outl" -> file "outl_dataset_noise_gauss_...txt"
output_prefix = input('Output filename prefix (e.g. outl): ', 's');
if isempty(output_prefix)
    output_prefix = '';
end

% --- Noise toggle ---
% 1 = add noise, 0 = skip. Default to 1 if the user enters non-numeric text.
add_noise = input('Add noise/outliers? (1=yes, 0=no): ', 's');
add_noise = str2double(add_noise);
if isnan(add_noise)
    add_noise = 1;
end

% Pack the toggle into a settings struct expected by adaptivehb.data.apply_noise.
settings = struct('addNoise', add_noise == 1);

% --- Noise parameters (only asked if noise is enabled) ---
if settings.addNoise
    % outlierFraction: what fraction of the N points should be corrupted.
    % E.g. 0.10 means 10% of points will receive noise.
    fraction = input('Fraction of points to corrupt (e.g. 0.10): ', 's');
    settings.outlierFraction = str2double(fraction);
    if isnan(settings.outlierFraction)
        settings.outlierFraction = 0.10;  % default: 10%
    end

    % outlierIntensity: amplitude of the injected noise relative to the
    % range of f values. E.g. 0.20 means noise magnitude ~ 20% of range(f).
    intensity = input('Relative noise intensity (e.g. 0.20): ', 's');
    settings.outlierIntensity = str2double(intensity);
    if isnan(settings.outlierIntensity)
        settings.outlierIntensity = 0.20;  % default: 20%
    end

    % seed: RNG seed for reproducibility. Same seed = same noise pattern.
    seed = input('RNG seed (e.g. 42): ', 's');
    settings.seed = str2double(seed);
    if isnan(settings.seed)
        settings.seed = 42;  % default seed
    end

    % noiseType: "gauss" for Gaussian noise (randn), "spike" for uniform noise (rand).
    noise_type = input('Noise type (gauss/spike): ', 's');
    if isempty(noise_type)
        noise_type = 'gauss';  % default: Gaussian
    end
    settings.noiseType = noise_type;
else
    % Inform the user that no noise will be applied.
    fprintf('Noise not applied: the dataset will remain unchanged.\n');
end

% =====================================================================
%  Apply noise using the adaptivehb.data.apply_noise function
% =====================================================================

% Extract the three columns from the loaded data matrix.
x = data(:, 1);      % x coordinates
y = data(:, 2);      % y coordinates
f_true = data(:, 3); % true function values (clean)

% Call the vector-level noise injector. Returns:
%   augmented : struct with fields x, y, fTrue, fNoisy, isOutlier
%   metadata  : struct with bookkeeping info (numOutliers, noiseLabel, etc.)
[augmented, metadata] = adaptivehb.data.apply_noise(x, y, f_true, settings);

% If the user chose not to add noise, exit early.
if ~metadata.noiseApplied
    fprintf('No noise applied: no additional file generated.\n');
    return;
end

% =====================================================================
%  Compose the output filename and write the augmented dataset
% =====================================================================

% Extract the base name of the input file (without path or extension).
[~, input_name, ~] = fileparts(input_path);

% Prepend the prefix with an underscore separator if provided.
prefix = output_prefix;
if ~isempty(prefix)
    prefix = [prefix '_'];
end

% Build a descriptive filename that encodes the noise parameters.
% Example: "outl_glacier_noise_gauss_frac0.10_int0.20_seed42.txt"
filename = adaptivehb.data.format_noise_filename(input_name, metadata, ...
    struct('Prefix', prefix));
output_path = fullfile(output_dir, filename);

% Assemble the 5-column output matrix:
%   col 1: x, col 2: y, col 3: fTrue, col 4: fNoisy, col 5: isOutlier (0/1).
data_to_write = [augmented.x, augmented.y, augmented.fTrue, augmented.fNoisy, augmented.isOutlier];

% Open the output file for writing. Abort with a clear error if it fails.
fid = fopen(output_path, 'w');
if fid == -1
    error('Cannot open file for writing: %s', output_path);
end

% onCleanup ensures the file is closed even if an error occurs during writing.
cleaner = onCleanup(@() fclose(fid));

% Write each row in scientific notation for the numeric columns and
% integer format for the outlier flag.
for i = 1:size(data_to_write, 1)
    fprintf(fid, '%.8e %.8e %.8e %.8e %d\n', data_to_write(i, 1), data_to_write(i, 2), ...
        data_to_write(i, 3), data_to_write(i, 4), data_to_write(i, 5));
end

% Explicitly clear the cleanup object to close the file immediately.
clear cleaner;

% =====================================================================
%  Print a summary to the command window
% =====================================================================
fprintf('Noisy dataset saved to: %s\n', output_path);
fprintf('Number of outliers: %d (%.2f%%)\n', metadata.numOutliers, metadata.actualOutlierFraction * 100);
