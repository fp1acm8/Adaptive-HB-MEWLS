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

clearvars; close all; clc;
format long;

% Ensure the package folder is on the MATLAB path.
thisFile = mfilename('fullpath');
projectRoot = fileparts(thisFile);
srcPath = fullfile(projectRoot, '..', 'src');
if ~isfolder(srcPath)
    error('Cannot find the src folder: %s', srcPath);
end
addpath(srcPath);

%% ================== PARAMETERS ==================
input_path = input('Path to the input dataset (format: x y f): ', 's');
if isempty(input_path)
    error('The input dataset path is required.');
end

if ~isfile(input_path)
    error('File not found: %s', input_path);
end

data = load(input_path);
if size(data, 2) ~= 3
    error('The file must contain exactly 3 columns: x y f');
end

output_dir = input('Output directory: ', 's');
if isempty(output_dir)
    output_dir = fileparts(input_path);
end

if ~isfolder(output_dir)
    error('Invalid output directory: %s', output_dir);
end

output_prefix = input('Output filename prefix (e.g. outl): ', 's');
if isempty(output_prefix)
    output_prefix = '';
end

add_noise = input('Add noise/outliers? (1=yes, 0=no): ', 's');
add_noise = str2double(add_noise);
if isnan(add_noise)
    add_noise = 1;
end

settings = struct('addNoise', add_noise == 1);

if settings.addNoise
    fraction = input('Fraction of points to corrupt (e.g. 0.10): ', 's');
    settings.outlierFraction = str2double(fraction);
    if isnan(settings.outlierFraction)
        settings.outlierFraction = 0.10;
    end

    intensity = input('Relative noise intensity (e.g. 0.20): ', 's');
    settings.outlierIntensity = str2double(intensity);
    if isnan(settings.outlierIntensity)
        settings.outlierIntensity = 0.20;
    end

    seed = input('RNG seed (e.g. 42): ', 's');
    settings.seed = str2double(seed);
    if isnan(settings.seed)
        settings.seed = 42;
    end

    noise_type = input('Noise type (gauss/spike): ', 's');
    if isempty(noise_type)
        noise_type = 'gauss';
    end
    settings.noiseType = noise_type;
else
    fprintf('Noise not applied: the dataset will remain unchanged.\n');
end

x = data(:, 1);
y = data(:, 2);
f_true = data(:, 3);

[augmented, metadata] = adaptivehb.data.apply_noise(x, y, f_true, settings);

if ~metadata.noiseApplied
    fprintf('No noise applied: no additional file generated.\n');
    return;
end

[~, input_name, ~] = fileparts(input_path);
prefix = output_prefix;
if ~isempty(prefix)
    prefix = [prefix '_'];
end

filename = adaptivehb.data.format_noise_filename(input_name, metadata, ...
    struct('Prefix', prefix));
output_path = fullfile(output_dir, filename);

data_to_write = [augmented.x, augmented.y, augmented.fTrue, augmented.fNoisy, augmented.isOutlier];

fid = fopen(output_path, 'w');
if fid == -1
    error('Cannot open file for writing: %s', output_path);
end
cleaner = onCleanup(@() fclose(fid));

for i = 1:size(data_to_write, 1)
    fprintf(fid, '%.8e %.8e %.8e %.8e %d\n', data_to_write(i, 1), data_to_write(i, 2), ...
        data_to_write(i, 3), data_to_write(i, 4), data_to_write(i, 5));
end

clear cleaner;

fprintf('Noisy dataset saved to: %s\n', output_path);
fprintf('Number of outliers: %d (%.2f%%)\n', metadata.numOutliers, metadata.actualOutlierFraction * 100);
