%% generate_test_datasets.m
% Script per generare i dataset degli esempi 11 e 13 dal paper THB-splines
% Salva i dati in formato txt compatibile con gli script HB-splines

clearvars; close all; clc;
format long;

%% ================== PARAMETRI GENERALI ==================
n_points = 150;  % Griglia uniforme 150×150
fprintf('Generazione dataset con griglia %dx%d\n', n_points, n_points);

%% ================== ESEMPIO 11 ==================
fprintf('\n=== ESEMPIO 11 ===\n');
fprintf('Dominio: [-1, 1] × [-1, 1]\n');
fprintf('Funzione: f = (2/3)*exp(-((10x-3)²+(10y-3)²)) + (2/3)*exp(-((10x+3)²+(10y+3)²)) + (2/3)*exp(-((10x)²+(10y)²))\n');

% Definizione del dominio
x11 = linspace(-1, 1, n_points);
y11 = linspace(-1, 1, n_points);
[X11, Y11] = meshgrid(x11, y11);

% Definizione della funzione dell'esempio 11
f11 = (2/3) * exp(-((10*X11 - 3).^2 + (10*Y11 - 3).^2)) + ...
      (2/3) * exp(-((10*X11 + 3).^2 + (10*Y11 + 3).^2)) + ...
      (2/3) * exp(-((10*X11).^2 + (10*Y11).^2));

% Conversione in formato colonna per salvataggio
data11 = [X11(:), Y11(:), f11(:)];

% Salvataggio del dataset
filename11 = '/MATLAB Drive/hierarchical_Bsplines_least_squares_with_GeoPDEs/data/example11_dataset.txt';
fid11 = fopen(filename11, 'w');
if fid11 == -1
    error('Impossibile creare il file %s', filename11);
end

% Scrivi i dati nel formato richiesto (x y f(x,y))
for i = 1:size(data11, 1)
    fprintf(fid11, '%.8e %.8e %.8e\n', data11(i, 1), data11(i, 2), data11(i, 3));
end
fclose(fid11);

fprintf('Dataset salvato in: %s\n', filename11);
fprintf('Numero di punti: %d\n', size(data11, 1));
fprintf('Range valori f: [%.6f, %.6f]\n', min(f11(:)), max(f11(:)));

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
% Condizione 1: f = 1 se y - x > 1/2
cond1 = (Y13 - X13) > 1/2;
f13(cond1) = 1;

% Condizione 2: f = 2(y - x) se 0 <= y - x <= 1/2
cond2 = (Y13 - X13) >= 0 & (Y13 - X13) <= 1/2;
f13(cond2) = 2 * (Y13(cond2) - X13(cond2));

% Condizione 3: f = 1/2 * cos(4π√q) + 1/2 se q <= 1/16
cond3 = q <= 1/16;
f13(cond3) = 0.5 * cos(4 * pi * sqrt(q(cond3))) + 0.5;

% Le altre regioni rimangono a 0 (condizione "otherwise")

% Conversione in formato colonna per salvataggio
data13 = [X13(:), Y13(:), f13(:)];

% Salvataggio del dataset
filename13 = '/MATLAB Drive/hierarchical_Bsplines_least_squares_with_GeoPDEs/data/example13_dataset.txt';
fid13 = fopen(filename13, 'w');
if fid13 == -1
    error('Impossibile creare il file %s', filename13);
end

% Scrivi i dati nel formato richiesto (x y f(x,y))
for i = 1:size(data13, 1)
    fprintf(fid13, '%.8e %.8e %.8e\n', data13(i, 1), data13(i, 2), data13(i, 3));
end
fclose(fid13);

fprintf('Dataset salvato in: %s\n', filename13);
fprintf('Numero di punti: %d\n', size(data13, 1));
fprintf('Range valori f: [%.6f, %.6f]\n', min(f13(:)), max(f13(:)));

%% ================== VISUALIZZAZIONE ==================
% Visualizzazione dei dataset generati
figure('Name', 'Dataset Esempio 11', 'Position', [100, 100, 800, 600]);

subplot(2, 2, 1);
surf(X11, Y11, f11, 'EdgeColor', 'none');
camlight;
title('Esempio 11 - Superficie');
xlabel('x'); ylabel('y'); zlabel('f(x,y)');
colorbar;

subplot(2, 2, 2);
contour(X11, Y11, f11, 20);
title('Esempio 11 - Contour');
xlabel('x'); ylabel('y');
colorbar;

subplot(2, 2, 3);
surf(X13, Y13, f13, 'EdgeColor', 'none');
camlight;
title('Esempio 13 - Superficie');
xlabel('x'); ylabel('y'); zlabel('f(x,y)');
colorbar;

subplot(2, 2, 4);
contour(X13, Y13, f13, 20);
title('Esempio 13 - Contour');
xlabel('x'); ylabel('y');
colorbar;

sgtitle('Dataset Generati - Esempi 11 e 13', 'FontSize', 14, 'FontWeight', 'bold');

%% ================== VERIFICA FORMATO ==================
fprintf('\n=== VERIFICA FORMATO ===\n');

% Verifica del formato del file (prime righe)
fprintf('Prime 5 righe del file %s:\n', filename11);
fid_test = fopen(filename11, 'r');
for i = 1:5
    line = fgetl(fid_test);
    fprintf('%s\n', line);
end
fclose(fid_test);

fprintf('\nPrime 5 righe del file %s:\n', filename13);
fid_test = fopen(filename13, 'r');
for i = 1:5
    line = fgetl(fid_test);
    fprintf('%s\n', line);
end
fclose(fid_test);

%% ================== STATISTICHE ==================
fprintf('\n=== STATISTICHE ===\n');

% Statistiche Esempio 11
fprintf('Esempio 11:\n');
fprintf('  Media: %.6f\n', mean(f11(:)));
fprintf('  Std Dev: %.6f\n', std(f11(:)));
fprintf('  Min: %.6f\n', min(f11(:)));
fprintf('  Max: %.6f\n', max(f11(:)));

% Statistiche Esempio 13
fprintf('Esempio 13:\n');
fprintf('  Media: %.6f\n', mean(f13(:)));
fprintf('  Std Dev: %.6f\n', std(f13(:)));
fprintf('  Min: %.6f\n', min(f13(:)));
fprintf('  Max: %.6f\n', max(f13(:)));

fprintf('\n=== GENERAZIONE COMPLETATA ===\n');
fprintf('I file sono pronti per essere utilizzati con gli script HB-splines.\n');
fprintf('Formato compatibile con black_forest.txt: x y f(x,y)\n');
