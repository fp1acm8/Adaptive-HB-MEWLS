%% add_noise_to_dataset.m
% Aggiunge rumore e outliers a un dataset esistente (formato: x y f)
% Salva il nuovo dataset in formato: x y f_true f_noisy is_outlier

clearvars; close all; clc;
format long;

%% ================== PARAMETRI ==================

% Path del dataset di input (deve avere 3 colonne: x y f)
input_path = '/MATLAB Drive/hierarchical_Bsplines_least_squares_with_GeoPDEs/data/black_forest.txt';

% Directory di output + nome file
output_dir = '/MATLAB Drive/hierarchical_Bsplines_least_squares_with_GeoPDEs/data';
output_filename_prefix = 'outl';

% Parametri rumore
add_noise = true;
outlier_fraction = 0.10;         % Frazione punti da corrompere (0.0 - 1.0)
outlier_intensity = 0.20;         % Intensità del rumore (in percentuale del range)
seed = 42;                       % Per riproducibilità
noise_type = 'gauss';            % 'gauss' (normale) o 'spike' (uniforme)

%% ================== CARICAMENTO DATI ==================

if ~isfile(input_path)
    error('File non trovato: %s', input_path);
end

data = load(input_path);
if size(data, 2) ~= 3
    error('Il file deve contenere esattamente 3 colonne: x y f');
end

x = data(:,1);
y = data(:,2);
f_true = data(:,3);

M = numel(f_true);
range_f = max(f_true) - min(f_true);

%% ================== AGGIUNTA RUMORE ==================

if add_noise
    rng(seed);
    num_outliers = round(outlier_fraction * M);
    outlier_idx = randperm(M, num_outliers);
    
    is_outlier = zeros(M,1);
    is_outlier(outlier_idx) = 1;

    f_noisy = f_true;

    switch lower(noise_type)
        case 'gauss'
            noise = outlier_intensity * range_f * randn(num_outliers, 1);
        case 'spike'
            noise = outlier_intensity * range_f * (2*rand(num_outliers, 1) - 1);
        otherwise
            error('Tipo di rumore non supportato: %s', noise_type);
    end

    f_noisy(outlier_idx) = f_noisy(outlier_idx) + noise;

    %% ================== SALVATAGGIO FILE ==================

    % Costruzione nome file output
    [~, input_name, ~] = fileparts(input_path);
    noise_info = sprintf('_noise_%s_frac%.2f_int%.2f_seed%d', ...
        noise_type, outlier_fraction, outlier_intensity, seed);
    output_name = sprintf('%s_%s%s.txt', output_filename_prefix, input_name, noise_info);
    output_path = fullfile(output_dir, output_name);

    % Salva in formato: x y f_true f_noisy is_outlier
    fid = fopen(output_path, 'w');
    if fid == -1
        error('Impossibile aprire il file per scrittura: %s', output_path);
    end

    for i = 1:M
        fprintf(fid, '%.8e %.8e %.8e %.8e %d\n', x(i), y(i), f_true(i), f_noisy(i), is_outlier(i));
    end

    fclose(fid);
    fprintf('Dataset rumoroso salvato in: %s\n', output_path);
    fprintf('Numero di outliers: %d (%.2f%%)\n', num_outliers, 100 * outlier_fraction);
else
    fprintf('Rumore non aggiunto. Setta `add_noise = true` per generare il dataset rumoroso.\n');
end
