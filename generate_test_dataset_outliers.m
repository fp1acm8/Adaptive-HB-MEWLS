%% generate_test_datasets.m
% Script per generare i dataset degli esempi 11 e 13 dal paper THB-splines
% Salva i dati in formato txt compatibile con gli script HB-splines
% Supporta generazione di versioni rumorose con outliers marcati
% Riutilizzabile: Input utente per parametri e opzioni

clearvars; close all; clc;
format long;

%% ================== INPUT UTENTE PER PARAMETRI ==================
% Parametri generali
n_points = input('Dimensione griglia (es. 150 per 150x150): ', 's');
n_points = str2double(n_points);
if isnan(n_points) || n_points <= 0
    n_points = 150;  % Default
end
fprintf('Generazione dataset con griglia %dx%d\n', n_points, n_points);

% Directory di salvataggio
output_dir = input('Directory di output (es. /MATLAB Drive/hierarchical_Bsplines_least_squares_with_GeoPDEs/data/): ', 's');
if isempty(output_dir)
    output_dir = '/MATLAB Drive/hierarchical_Bsplines_least_squares_with_GeoPDEs/data/';
end

% Esempi da generare
examples = input('Esempi da generare (11, 13, o entrambi separati da spazio): ', 's');
examples = str2num(examples);  %#ok<ST2NM> % Converti a numeri
if isempty(examples)
    examples = [11 13];  % Default: entrambi
end

% Opzione per rumore
add_noise = input('Aggiungere rumore/outliers? (1=si, 0=no): ', 's');
add_noise = str2double(add_noise);
if add_noise == 1
    outlier_fraction = input('Frazione outliers (es. 0.05 per 5%): ', 's');
    outlier_fraction = str2double(outlier_fraction);
    if isnan(outlier_fraction)
        outlier_fraction = 0.05;
    end
    
    outlier_intensity = input('Intensità rumore (es. 0.1 per ±10% range): ', 's');
    outlier_intensity = str2double(outlier_intensity);
    if isnan(outlier_intensity)
        outlier_intensity = 0.1;
    end
    
    seed = input('Seed per riproducibilità (es. 42): ', 's');
    seed = str2double(seed);
    if isnan(seed)
        seed = 42;
    end
    
    noise_type = input('Tipo di rumore (gauss=gaussiano, spike=estremi): ', 's');
    if isempty(noise_type)
        noise_type = 'gauss';
    end
else
    add_noise = 0;
end

%% ================== GENERAZIONE DATASET ==================
% Inizializza variabili specifiche per visualizzazione
generate11 = any(examples == 11);
generate13 = any(examples == 13);

if generate11
    %% ================== ESEMPIO 11 ==================
    fprintf('\n=== ESEMPIO 11 ===\n');
    fprintf('Dominio: [-1, 1] × [-1, 1]\n');
    fprintf('Funzione: f = (2/3)*exp(-((10x-3)²+(10y-3)²)) + (2/3)*exp(-((10x+3)²+(10y+3)²)) + (2/3)*exp(-((10x)²+(10y)²))\n');
    
    % Definizione del dominio
    x11 = linspace(-1, 1, n_points);
    y11 = linspace(-1, 1, n_points);
    [X11, Y11] = meshgrid(x11, y11);
    
    % Definizione della funzione
    f11 = (2/3) * exp(-((10*X11 - 3).^2 + (10*Y11 - 3).^2)) + ...
           (2/3) * exp(-((10*X11 + 3).^2 + (10*Y11 + 3).^2)) + ...
           (2/3) * exp(-((10*X11).^2 + (10*Y11).^2));
    
    % Genera e salva versione pulita
    data_ex = [X11(:), Y11(:), f11(:)];
    filename_clean11 = fullfile(output_dir, 'example11_dataset.txt');
    save_dataset(filename_clean11, data_ex, false);  % Senza mark
    
    fprintf('Dataset pulito salvato in: %s\n', filename_clean11);
    fprintf('Numero di punti: %d\n', size(data_ex, 1));
    fprintf('Range valori f: [%.6f, %.6f]\n', min(f11(:)), max(f11(:)));
    
    % Genera e salva versione rumorosa (se richiesta)
    if add_noise
        rng(seed);  % Per riproducibilità
        M = numel(f11);
        range_f = max(f11(:)) - min(f11(:));
        num_outliers = round(outlier_fraction * M);
        outlier_idx = randperm(M, num_outliers);
        is_outlier11 = zeros(M, 1);
        is_outlier11(outlier_idx) = 1;
        
        f11_noisy = f11(:);
        if strcmp(noise_type, 'gauss')
            noise = outlier_intensity * range_f * randn(num_outliers, 1);
        elseif strcmp(noise_type, 'spike')
            noise = outlier_intensity * range_f * (2 * rand(num_outliers, 1) - 1);
        else
            error('Tipo di rumore non supportato');
        end
        f11_noisy(outlier_idx) = f11_noisy(outlier_idx) + noise;
        
        data_noisy = [X11(:), Y11(:), f11(:), f11_noisy, is_outlier11];
        
        % Nome file con info rumore
        noise_info = sprintf('_noise_%s_frac%.2f_int%.2f_seed%d', noise_type, outlier_fraction, outlier_intensity, seed);
        filename_noisy11 = fullfile(output_dir, ['example11_dataset' noise_info '.txt']);
        save_dataset(filename_noisy11, data_noisy, true);  % Con mark
        
        fprintf('Dataset rumoroso salvato in: %s\n', filename_noisy11);
        fprintf('Numero di outliers: %d (%.2f%%)\n', num_outliers, outlier_fraction*100);
        fprintf('Range valori f_noisy: [%.6f, %.6f]\n', min(f11_noisy), max(f11_noisy));
    end
end

if generate13
    %% ================== ESEMPIO 13 ==================
    fprintf('\n=== ESEMPIO 13 ===\n');
    fprintf('Dominio: [0, 2] × [0, 1]\n');
    fprintf('Funzione: piecewise function con condizioni multiple\n');
    
    % Definizione del dominio
    x13 = linspace(0, 2, n_points);
    y13 = linspace(0, 1, n_points);
    [X13, Y13] = meshgrid(x13, y13);
    
    % Inizializzazione della funzione
    f13 = zeros(size(X13));
    
    % Definizione della funzione q(x,y)
    q = (X13 - 3/2).^2 + (Y13 - 1/2).^2;
    
    % Applicazione delle condizioni piecewise
    cond1 = (Y13 - X13) > 1/2; f13(cond1) = 1;
    cond2 = (Y13 - X13) >= 0 & (Y13 - X13) <= 1/2; f13(cond2) = 2 * (Y13(cond2) - X13(cond2));
    cond3 = q <= 1/16; f13(cond3) = 0.5 * cos(4 * pi * sqrt(q(cond3))) + 0.5;
    
    % Genera e salva versione pulita
    data_ex = [X13(:), Y13(:), f13(:)];
    filename_clean13 = fullfile(output_dir, 'example13_dataset.txt');
    save_dataset(filename_clean13, data_ex, false);  % Senza mark
    
    fprintf('Dataset pulito salvato in: %s\n', filename_clean13);
    fprintf('Numero di punti: %d\n', size(data_ex, 1));
    fprintf('Range valori f: [%.6f, %.6f]\n', min(f13(:)), max(f13(:)));
    
    % Genera e salva versione rumorosa (se richiesta)
    if add_noise
        rng(seed);  % Per riproducibilità
        M = numel(f13);
        range_f = max(f13(:)) - min(f13(:));
        num_outliers = round(outlier_fraction * M);
        outlier_idx = randperm(M, num_outliers);
        is_outlier13 = zeros(M, 1);
        is_outlier13(outlier_idx) = 1;
        
        f13_noisy = f13(:);
        if strcmp(noise_type, 'gauss')
            noise = outlier_intensity * range_f * randn(num_outliers, 1);
        elseif strcmp(noise_type, 'spike')
            noise = outlier_intensity * range_f * (2 * rand(num_outliers, 1) - 1);
        else
            error('Tipo di rumore non supportato');
        end
        f13_noisy(outlier_idx) = f13_noisy(outlier_idx) + noise;
        
        data_noisy = [X13(:), Y13(:), f13(:), f13_noisy, is_outlier13];
        
        % Nome file con info rumore
        noise_info = sprintf('_noise_%s_frac%.2f_int%.2f_seed%d', noise_type, outlier_fraction, outlier_intensity, seed);
        filename_noisy13 = fullfile(output_dir, ['example13_dataset' noise_info '.txt']);
        save_dataset(filename_noisy13, data_noisy, true);  % Con mark
        
        fprintf('Dataset rumoroso salvato in: %s\n', filename_noisy13);
        fprintf('Numero di outliers: %d (%.2f%%)\n', num_outliers, outlier_fraction*100);
        fprintf('Range valori f_noisy: [%.6f, %.6f]\n', min(f13_noisy), max(f13_noisy));
    end
end

%% ================== VISUALIZZAZIONE ==================
if generate11
    % Visualizza Esempio 11 (pulito e rumoroso se presente)
    figure('Name', 'Dataset Esempio 11', 'Position', [100, 100, 800, 600]);
    
    subplot(2, 2, 1);
    surf(X11, Y11, f11, 'EdgeColor', 'none');
    camlight;
    title('Esempio 11 Pulito - Superficie');
    xlabel('x'); ylabel('y'); zlabel('f(x,y)');
    colorbar;
    
    subplot(2, 2, 2);
    contour(X11, Y11, f11, 20);
    title('Esempio 11 Pulito - Contour');
    xlabel('x'); ylabel('y');
    colorbar;
    
    if add_noise
        subplot(2, 2, 3);
        surf(X11, Y11, reshape(f11_noisy, size(f11)), 'EdgeColor', 'none');
        camlight;
        hold on;
        outlier_points = [X11(:), Y11(:), f11_noisy, is_outlier11];
        outlier_points = outlier_points(is_outlier11 == 1, :);
        scatter3(outlier_points(:,1), outlier_points(:,2), outlier_points(:,3), 50, 'r', 'filled', 'Marker', 'o');
        title('Esempio 11 Rumoroso - Superficie con Outliers Rossi');
        xlabel('x'); ylabel('y'); zlabel('f_noisy');
        colorbar;
        
        subplot(2, 2, 4);
        contour(X11, Y11, reshape(f11_noisy, size(f11)), 20);
        hold on;
        scatter(outlier_points(:,1), outlier_points(:,2), 50, 'r', 'filled');
        title('Esempio 11 Rumoroso - Contour con Outliers Rossi');
        xlabel('x'); ylabel('y');
        colorbar;
    end
    
    sgtitle('Dataset Esempio 11 - Pulito e Rumoroso', 'FontSize', 14, 'FontWeight', 'bold');
end

if generate13
    % Visualizza Esempio 13 (pulito e rumoroso se presente)
    figure('Name', 'Dataset Esempio 13', 'Position', [100, 100, 800, 600]);
    
    subplot(2, 2, 1);
    surf(X13, Y13, f13, 'EdgeColor', 'none');
    camlight;
    title('Esempio 13 Pulito - Superficie');
    xlabel('x'); ylabel('y'); zlabel('f(x,y)');
    colorbar;
    
    subplot(2, 2, 2);
    contour(X13, Y13, f13, 20);
    title('Esempio 13 Pulito - Contour');
    xlabel('x'); ylabel('y');
    colorbar;
    
    if add_noise
        subplot(2, 2, 3);
        surf(X13, Y13, reshape(f13_noisy, size(f13)), 'EdgeColor', 'none');
        camlight;
        hold on;
        outlier_points = [X13(:), Y13(:), f13_noisy, is_outlier13];
        outlier_points = outlier_points(is_outlier13 == 1, :);
        scatter3(outlier_points(:,1), outlier_points(:,2), outlier_points(:,3), 50, 'r', 'filled', 'Marker', 'o');
        title('Esempio 13 Rumoroso - Superficie con Outliers Rossi');
        xlabel('x'); ylabel('y'); zlabel('f_noisy');
        colorbar;
        
        subplot(2, 2, 4);
        contour(X13, Y13, reshape(f13_noisy, size(f13)), 20);
        hold on;
        scatter(outlier_points(:,1), outlier_points(:,2), 50, 'r', 'filled');
        title('Esempio 13 Rumoroso - Contour con Outliers Rossi');
        xlabel('x'); ylabel('y');
        colorbar;
    end
    
    sgtitle('Dataset Esempio 13 - Pulito e Rumoroso', 'FontSize', 14, 'FontWeight', 'bold');
end

%% ================== VERIFICA FORMATO ==================
fprintf('\n=== VERIFICA FORMATO ===\n');

% Verifica file pulito e rumoroso per Esempio 11
if generate11
    fprintf('Prime 5 righe del file pulito %s:\n', filename_clean11);
    fid_test = fopen(filename_clean11, 'r');
    for i = 1:5
        line = fgetl(fid_test);
        fprintf('%s\n', line);
    end
    fclose(fid_test);
    
    if add_noise
        fprintf('\nPrime 5 righe del file rumoroso %s:\n', filename_noisy11);
        fid_test = fopen(filename_noisy11, 'r');
        for i = 1:5
            line = fgetl(fid_test);
            fprintf('%s\n', line);
        end
        fclose(fid_test);
    end
end

% Verifica file pulito e rumoroso per Esempio 13
if generate13
    fprintf('Prime 5 righe del file pulito %s:\n', filename_clean13);
    fid_test = fopen(filename_clean13, 'r');
    for i = 1:5
        line = fgetl(fid_test);
        fprintf('%s\n', line);
    end
    fclose(fid_test);
    
    if add_noise
        fprintf('\nPrime 5 righe del file rumoroso %s:\n', filename_noisy13);
        fid_test = fopen(filename_noisy13, 'r');
        for i = 1:5
            line = fgetl(fid_test);
            fprintf('%s\n', line);
        end
        fclose(fid_test);
    end
end

%% ================== STATISTICHE ==================
fprintf('\n=== STATISTICHE ===\n');

if generate11
    fprintf('Esempio 11 Pulito:\n');
    fprintf('  Media: %.6f\n', mean(f11(:)));
    fprintf('  Std Dev: %.6f\n', std(f11(:)));
    fprintf('  Min: %.6f\n', min(f11(:)));
    fprintf('  Max: %.6f\n', max(f11(:)));
    if add_noise
        fprintf('Esempio 11 Rumoroso:\n');
        fprintf('  Media: %.6f\n', mean(f11_noisy));
        fprintf('  Std Dev: %.6f\n', std(f11_noisy));
        fprintf('  Min: %.6f\n', min(f11_noisy));
        fprintf('  Max: %.6f\n', max(f11_noisy));
    end
end

if generate13
    fprintf('Esempio 13 Pulito:\n');
    fprintf('  Media: %.6f\n', mean(f13(:)));
    fprintf('  Std Dev: %.6f\n', std(f13(:)));
    fprintf('  Min: %.6f\n', min(f13(:)));
    fprintf('  Max: %.6f\n', max(f13(:)));
    if add_noise
        fprintf('Esempio 13 Rumoroso:\n');
        fprintf('  Media: %.6f\n', mean(f13_noisy));
        fprintf('  Std Dev: %.6f\n', std(f13_noisy));
        fprintf('  Min: %.6f\n', min(f13_noisy));
        fprintf('  Max: %.6f\n', max(f13_noisy));
    end
end

fprintf('\n=== GENERAZIONE COMPLETATA ===\n');
fprintf('I file sono pronti per essere utilizzati con gli script HB-splines.\n');
fprintf('Formato file pulito: x y f(x,y)\n');
fprintf('Formato file rumoroso: x y f_true f_noisy is_outlier (1=outlier, 0=no)\n');

%% ================== FUNZIONE HELPER PER SALVATAGGIO ==================
function save_dataset(filename, data, with_mark)
    fid = fopen(filename, 'w');
    if fid == -1
        error('Impossibile creare il file %s', filename);
    end
    
    if with_mark
        % Salva x y f_true f_noisy is_outlier
        for i = 1:size(data, 1)
            fprintf(fid, '%.8e %.8e %.8e %.8e %d\n', data(i, 1), data(i, 2), data(i, 3), data(i, 4), data(i, 5));
        end
    else
        % Salva x y f
        for i = 1:size(data, 1)
            fprintf(fid, '%.8e %.8e %.8e\n', data(i, 1), data(i, 2), data(i, 3));
        end
    end
    fclose(fid);
end
