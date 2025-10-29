%% generate_noisy_dataset.m
% Wrapper script to generate a noisy version of an existing dataset using
% adaptivehb.data.apply_noise. This script keeps the interactive workflow of
% the original add_noise_to_dataset.m while relying on reusable functions
% that can also operate in-memory.

clearvars; close all; clc;
format long;

% Ensure the package folder is on the MATLAB path
thisFile = mfilename('fullpath');
projectRoot = fileparts(thisFile);
srcPath = fullfile(projectRoot, '..', 'src');
if ~isfolder(srcPath)
    error('Impossibile trovare la cartella src: %s', srcPath);
end
addpath(srcPath);

%% ================== PARAMETRI ==================
input_path = input('Path del dataset di input (formato x y f): ', 's');
if isempty(input_path)
    error('Il path del dataset di input è obbligatorio.');
end

if ~isfile(input_path)
    error('File non trovato: %s', input_path);
end

data = load(input_path);
if size(data, 2) ~= 3
    error('Il file deve contenere esattamente 3 colonne: x y f');
end

output_dir = input('Directory di output: ', 's');
if isempty(output_dir)
    output_dir = fileparts(input_path);
end

if ~isfolder(output_dir)
    error('Directory di output non valida: %s', output_dir);
end

output_prefix = input('Prefisso per il file di output (es. outl): ', 's');
if isempty(output_prefix)
    output_prefix = '';
end

add_noise = input('Aggiungere rumore/outliers? (1=si, 0=no): ', 's');
add_noise = str2double(add_noise);
if isnan(add_noise)
    add_noise = 1;
end

settings = struct('addNoise', add_noise == 1);

if settings.addNoise
    fraction = input('Frazione punti da corrompere (es. 0.10): ', 's');
    settings.outlierFraction = str2double(fraction);
    if isnan(settings.outlierFraction)
        settings.outlierFraction = 0.10;
    end

    intensity = input('Intensità rumore relativa (es. 0.20): ', 's');
    settings.outlierIntensity = str2double(intensity);
    if isnan(settings.outlierIntensity)
        settings.outlierIntensity = 0.20;
    end

    seed = input('Seed RNG (es. 42): ', 's');
    settings.seed = str2double(seed);
    if isnan(settings.seed)
        settings.seed = 42;
    end

    noise_type = input('Tipo di rumore (gauss/spike): ', 's');
    if isempty(noise_type)
        noise_type = 'gauss';
    end
    settings.noiseType = noise_type;
else
    fprintf('Rumore non applicato: il dataset resterà invariato.\n');
end

x = data(:, 1);
y = data(:, 2);
f_true = data(:, 3);

[augmented, metadata] = adaptivehb.data.apply_noise(x, y, f_true, settings);

if ~metadata.noiseApplied
    fprintf('Nessun rumore applicato: nessun file aggiuntivo generato.\n');
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
    error('Impossibile aprire il file per scrittura: %s', output_path);
end
cleaner = onCleanup(@() fclose(fid));

for i = 1:size(data_to_write, 1)
    fprintf(fid, '%.8e %.8e %.8e %.8e %d\n', data_to_write(i, 1), data_to_write(i, 2), ...
        data_to_write(i, 3), data_to_write(i, 4), data_to_write(i, 5));
end

clear cleaner;

fprintf('Dataset rumoroso salvato in: %s\n', output_path);
fprintf('Numero di outliers: %d (%.2f%%)\n', metadata.numOutliers, metadata.actualOutlierFraction * 100);
