%% main_script_weighted_least_squares.m
% Script MATLAB per HB-splines + Weighted Least Squares con termine di penalizzazione
clearvars; close all; clc; format long;
DEBUG = true;  % Flag per output dettagliato

%% ================== PARAMETRI MODIFICABILI ==================
% File dataset (modifica qui per cambiare dataset)
filename = '/MATLAB Drive/hierarchical_Bsplines_least_squares_with_GeoPDEs/data/example13_dataset_noise_gauss_frac0.05_int0.10_seed18.txt';

% Input per tipo dataset
outlier_dataset = input('Dataset con outliers? (0=no (3 colonne), 1=si (5 colonne)): ');

% Parametri del problema (modifica se necessario)
problem_data = struct(...
    'geo_name',     'geo_square.txt', ...
    'nmnn_sides',   [], ...
    'drchlt_sides', [1 2 3 4], ...
    'c_diff',       @(x, y) ones(size(x)));

C = 100;  % Costante per soluzione esatta
normax2 = @(x,y) ((x-.5).^2+(y-.5).^2);
problem_data.uex = @(x,y) exp(-C*normax2(x,y));
problem_data.f = @(x,y) 4*C*(1-C*normax2(x,y)).*problem_data.uex(x,y);
problem_data.g = @(x, y, ind) zeros(size(x));
problem_data.h = @(x, y, ind) problem_data.uex(x,y);
problem_data.graduex = @(x,y) -2*C*cat (1, ...
            reshape (problem_data.uex(x,y).*(x-.5), [1, size(x)]), ...
            reshape (problem_data.uex(x,y).*(y-.5), [1, size(x)]));

% Parametri dello spazio discreto
method_data = struct(...
    'degree',      [2 2], ...     % B-spline biquadratic
    'regularity',  [1 1], ...     % Continuità C^1
    'nsub_coarse', [16 8], ...    % Suddivisione iniziale n×n elementi
    'nsub_refine', [2 2], ...     % Ogni raffinamento dimezza ogni spigolo
    'nquad',       [3 3], ...     % 3×3 punti di quadratura
    'space_type',  'standard', ...
    'truncated',   1);            % Truncated THB-splines

% Parametri adattività
adaptivity_data = struct(...
    'flag',           'functions', ...   % Marcatura basata su funzioni
    'mark_param',     0.5, ...           
    'mark_strategy',  'MS', ...   
    'max_level',      6, ...             % Profondità massima dell'albero
    'max_ndof',       50000, ...         % Numero max di dof consentiti
    'num_max_iter',   15, ...            % Iterazioni massime di adattività
    'max_nel',        50000, ...         % Numero max di elementi
    'tol',            5e-3, ...          % Tolleranza globale sull'errore
    'adm',            0);

% Parametri weighted least squares
lambda_pen = 1e-12;  % Parametro di penalizzazione

%% --------------------- CARICAMENTO DATI SPARSE -------------------------
% Caricamento dati
fid = fopen(filename, 'r');
if fid == -1
    error('File %s non trovato', filename);
end
raw_vector = fscanf(fid, '%f');  % Legge tutti i numeri in un vettore lineare
fclose(fid);

total_elements = length(raw_vector);

% Determina il numero di colonne atteso basato su outlier_dataset
if outlier_dataset == 0
    expected_cols = 3;
elseif outlier_dataset == 1
    expected_cols = 5;
else
    error('Valore non valido per outlier_dataset: deve essere 0 o 1');
end

% Verifica e reshape
if mod(total_elements, expected_cols) ~= 0
    error('Formato file non compatibile: numero totale elementi %d non divisibile per %d', total_elements, expected_cols);
end

M = total_elements / expected_cols;
raw = reshape(raw_vector, expected_cols, M);

% Assegna variabili in base al formato
data0 = raw(1:2, :)';
if expected_cols == 3
    % Dataset pulito: x y f
    f_true = raw(3, :)';
    f_noise = f_true;  % Nessun rumore
    is_outlier = zeros(size(f_true));  % Nessun outlier
elseif expected_cols == 5
    % Dataset rumoroso: x y f_true f_noise is_outlier
    f_true = raw(3, :)';
    f_noise = raw(4, :)';  % Usa f_noise per fitting
    is_outlier = raw(5, :)';
end

% Mappa le coordinate in [0,1]^2
a = min(data0(:,1)); b = max(data0(:,1));
c = min(data0(:,2)); d = max(data0(:,2));
data = [(data0(:,1)-a)/(b-a), (data0(:,2)-c)/(d-c)];
data0 = data;  % Salva originali per plot (corretto da codice originale)
f = f_noise;  % Usa f_noise per fitting
M = size(data,1);

fprintf('Dataset caricato: %s (formato: %d colonne)\n', filename, expected_cols);
fprintf('Numero di punti: %d\n', M);
if expected_cols == 5
    fprintf('Outliers marcati: %d (%.2f%%)\n', sum(is_outlier), 100*sum(is_outlier)/M);
end

%% ---------- INIZIALIZZAZIONE SPAZIO GERARCHICO ----------
[hmsh, hspace, geometry] = adaptivity_initialize_laplace(problem_data, method_data);

%% -------- PARAMETRI WEIGHTED LEAST SQUARES -------------------------
weight = ones(1,M) / M;  % Inizializzazione pesi uniformi

%% ================== CICLO ADATTATIVO ===========================
results = [];  % Matrice per salvare [livello, ndof, max_err_all, L2_err_all, max_err_no_out, L2_err_no_out, cond_num]
tol_sat = 0;   % Flag tolleranza soddisfatta
stuck = 0;     % Flag spazio invariato

while true
    % Estrai livello e grado di libertà correnti
    lev = hspace.nlevels;
    ndof = hspace.ndof;
    
    if DEBUG
        fprintf('\n=== Livello %d | ndof = %d ===\n', lev, ndof); 
    end
    
    % Calcola coefficienti con Weighted Least Squares + penalizzazione
    if DEBUG
        disp('Computing weighted least squares solution...')
    end
    [QI_coeff, cond_num] = getcoeff_weighted_least_squares_pen(...
        hspace, hmsh, data, f, lambda_pen, weight);
    
    if DEBUG
        disp('Solution computed...')
        disp('Evaluating solution at data points...')
    end
    
    % Valuta la soluzione sui punti dati
    QI = sp_eval_alt(QI_coeff, hspace, data);
    
    % Calcola errori rispetto a f_true
    error = abs(f_true - QI');
    
    % Metriche globali di errore (su tutti i punti)
    max_err_all = max(error);
    L2_err_all = norm(error,2)/sqrt(length(error));  % RMSE su tutti
    
    % Metriche escludendo outliers (is_outlier == 1)
    no_out_mask = (is_outlier == 0);
    if any(no_out_mask)
        error_no_out = error(no_out_mask);
        max_err_no_out = max(error_no_out);
        L2_err_no_out = norm(error_no_out,2)/sqrt(length(error_no_out));
    else
        max_err_no_out = NaN;
        L2_err_no_out = NaN;
    end
    
    % Salva risultati
    results = [results; lev, ndof, max_err_all, L2_err_all, max_err_no_out, L2_err_no_out, cond_num];
    
    if DEBUG
        fprintf('cond(S)=%.2e | max_err_all=%.2e | L2_err_all=%.2e\n', ...
                cond_num, max_err_all, L2_err_all);
        fprintf('max_err_no_out=%.2e | L2_err_no_out=%.2e\n', max_err_no_out, L2_err_no_out);
    end
    
    % Condizioni di arresto
    if lev > 1
        prev_err = results(end-1, 4);  % errore L2_all precedente
        improvement = (prev_err - L2_err_all) / prev_err;
        if improvement < 0  % miglioramento < 0%
            fprintf('> Miglioramento insufficiente (%.3f%%) -> uscita\n', improvement*100);
            break;
        end
    end
    
    if max_err_all <= adaptivity_data.tol
        disp('> Tolleranza globale soddisfatta'); 
        tol_sat = 1;
        break;
    end
    
    if ndof >= adaptivity_data.max_ndof
        warning('Raggiunto ndof max'); 
        break;
    end
    
    if lev >= adaptivity_data.max_level
        warning('Raggiunto livello max'); 
        break;
    end
    
    % Trova punti con errore superiore alla tolleranza (su error da f_true)
    ex_tol1 = find(error > adaptivity_data.tol);
    
    if DEBUG
        fprintf('→ marcati %d punti su %d (%.2f%%) con err > %.2f\n', ...
                numel(ex_tol1), M, 100*numel(ex_tol1)/M, adaptivity_data.tol);
    end
    
    if isempty(ex_tol1)
        disp('> Nessun punto da raffinare -> uscita'); 
        break;
    end
    
    % Aggiorna pesi per i punti con errore elevato
    %weight(ex_tol1) = weight(ex_tol1) .* (1 + error(ex_tol1)');
    
    % Trova supporti da raffinare
    marked = support_containing_point(hspace, hmsh, data(ex_tol1,:));
    
    % Salva stato prima del raffinamento
    hmsh_temp = hmsh;
    hspace_temp = hspace;
    
    % Raffina mesh e spazio
    %marked_enlarged = apply_offset_strategy(marked, hspace, hmsh, method_data.degree);
    %[hmsh,hspace] = adaptivity_refine(hmsh, hspace, marked_enlarged, adaptivity_data);
    [hmsh, hspace] = adaptivity_refine(hmsh, hspace, marked, adaptivity_data);
    
    if DEBUG
        fprintf('→ raffinati %d supporti | nuovo ndof = %d\n', ...
                numel(marked{lev}), hspace.ndof);
    end
    
    % Controlla se lo spazio è rimasto invariato
    if hspace.ndof == hspace_temp.ndof
        stuck = 1;
        disp('> Spazio invariato -> uscita');
        break;
    end
end

%% -------------------- RISULTATI --------------------------------
% Stampa tabella comparativa
fprintf('\n%s\n', repmat('=', 1, 80));
fprintf('RISULTATI FINALI - WEIGHTED LEAST SQUARES\n');
fprintf('%s\n', repmat('=', 1, 80));

if outlier_dataset == 0
    % Tabella standard senza outliers
    fprintf('Liv  DOF   Max Err    L2 Err     cond_num\n');
    fprintf('%s\n', repmat('-', 1, 50));
    fprintf('%3d %5d %10.4e %10.4e %10.2e\n', results(:, [1,2,3,4,7])');
else
    % Tabella completa con outliers
    fprintf('Liv  ndof   max_err_all    L2_err_all   max_err_no_out    L2_err_no_out   cond_num\n');
    fprintf('%s\n', repmat('-', 1, 80));
    fprintf('%3d %5d %14.4e %14.4e %14.4e %14.4e %10.2e\n', results');
end

fprintf('\n%s\n', repmat('=', 1, 50));
fprintf('ELABORAZIONE COMPLETATA\n');
fprintf('Dataset: %s\n', filename);
fprintf('Livelli finali: %d\n', size(results,1));
fprintf('DOF finali: %d\n', results(end,2));
fprintf('Errore finale (all): %.4e\n', results(end,4));
fprintf('%s\n', repmat('=', 1, 50));

%% -------------------- PLOT FINALI ------------------------------
% % Controlla se il dominio è quadrato per axis equal
is_square_domain = abs((b - a) - (d - c)) < 1e-6;  % Confronta lunghezze assi
% 
% Mesh adattata
figure('Name', 'Mesh');
hmsh_plot_cells(hmsh, 10, 1);
view(0, 90); 
if is_square_domain
    axis equal;
end
axis off; 
% 
% % Superficie della soluzione
% figure('Name', 'Superficie WLS');
% xplot = linspace(a, b, 151);
% yplot = linspace(c, d, 151);
% [Xplot, Yplot] = meshgrid(xplot, yplot);
% try
%     QIplot = sp_eval(QI_coeff', hspace_temp, geometry, [151 151]);
% catch
%     QIplot = sp_eval(QI_coeff', hspace, geometry, [151 151]);
% end
% QIplot = QIplot';
% mesh(Xplot, Yplot, QIplot, 'FaceLighting', 'phong', 'FaceColor', [.5 .5 .5], ...
%      'EdgeColor', 'none', 'AmbientStrength', 1.0, 'SpecularExponent', 15, ...
%      'SpecularStrength', 0);
% camlight(-37.5, 30+20);
% title('Superficie Weighted Least Squares');
% xlabel('X'); ylabel('Y'); zlabel('U');
% if is_square_domain
%     axis equal;
% end
% axis off;
% 
% % Grafico 3D dei pesi (invece di errore, visualizza pesi weight con outliers evidenziati se presenti)
% figure('Name', 'Weights 3D');
% scatter3(data0(:,1), data0(:,2), weight', 20, weight', 'filled');  % weight è 1xM, trasponi
% hold on;
% if outlier_dataset == 1
%     outlier_mask = (is_outlier == 1);
%     scatter3(data0(outlier_mask,1), data0(outlier_mask,2), weight(outlier_mask)', 30, 'r', 'filled', 'MarkerEdgeColor', 'k');
% end
% colorbar;
% title('Pesi Finali in 3D (Outliers Rossi se presenti)');
% xlabel('X'); ylabel('Y'); zlabel('Peso');
% 
% % Dati originali 3D (mostra f_noise con outliers evidenziati se presenti)
% figure('Name', 'Dati Originali');
% scatter3(data0(:,1), data0(:,2), f, 20, f, 'filled');
% hold on;
% if outlier_dataset == 1
%     scatter3(data0(outlier_mask,1), data0(outlier_mask,2), f(outlier_mask), 30, 'r', 'filled', 'MarkerEdgeColor', 'k');
% end
% colorbar;
% %title('Dati originali');
% xlabel('X'); ylabel('Y'); zlabel('f');
% 
% % Punti originali (proiezione sul piano xy)
% figure('Name', 'Punti Dati');
% plot(data0(:,1), data0(:,2), '.');
% axis([a b c d]);
% %title('Distribuzione punti dati');
% xlabel('X'); ylabel('Y');
% 
% % Distribuzione dei pesi
% figure('Name', 'Distribuzione Pesi');
% histogram(weight, 50);
% %title('Distribuzione finale dei pesi');
% xlabel('Peso'); ylabel('Frequenza');
% 
% % Grafico di convergenza
% if size(results, 1) > 1
%     figure('Name', 'Convergenza');
%     subplot(2,2,1);
%     semilogy(results(:,1), results(:,3), 'o-');
%     title('Errore Massimo vs Livello');
%     xlabel('Livello'); ylabel('Errore Max');
%     grid on;
% 
%     subplot(2,2,2);
%     semilogy(results(:,1), results(:,4), 's-');
%     title('Errore L2 vs Livello');
%     xlabel('Livello'); ylabel('Errore L2');
%     grid on;
% 
%     subplot(2,2,3);
%     semilogy(results(:,2), results(:,3), 'o-');
%     title('Errore Massimo vs DOF');
%     xlabel('DOF'); ylabel('Errore Max');
%     grid on;
% 
%     subplot(2,2,4);
%     semilogy(results(:,2), results(:,4), 's-');
%     title('Errore L2 vs DOF');
%     xlabel('DOF'); ylabel('Errore L2');
%     grid on;
% end

%% Strategia di raffinamento con offset ring (S2)
% function marked_enlarged = apply_offset_strategy(marked, hspace, hmsh, degree)
%     % Calcola offset ring size
%     d = degree(1);  % assumendo grado uguale nelle due direzioni
%     offset_size = ceil((d-1)/4);
%     if offset_size == 0
%         offset_size = 1;  % minimo 1 per degree=2
%     end
% 
%     marked_enlarged = marked;
% 
%     % Per ogni livello nella gerarchia
%     for lev = 1:length(marked)
%         if ~isempty(marked{lev})
%             % Espandi la regione marcata con offset ring
%             marked_enlarged{lev} = enlarge_marked_region(...
%                 marked{lev}, hspace, hmsh, lev, offset_size);
%         end
%     end
% end
% 
% function enlarged_region = enlarge_marked_region(marked_cells, hspace, hmsh, level, offset)
%     % Implementa l'enlargement con offset ring Questa funzione dovrebbe: 1.
%     % Identificare le celle adiacenti alle celle marcate 2. Espandere la
%     % regione di 'offset' celle in tutte le direzioni 3. Assicurarsi che la
%     % regione risultante sia abbastanza grande
%     %    da contenere almeno una funzione di base del livello più fine
% 
%     % Implementazione specifica dipende dalla struttura di hmsh
%     enlarged_region = marked_cells;  % placeholder
% end
%% -------------------- RISULTATI FINALI ------------------------------
%create_formatted_table_image(results, 'LS', 'black_forest_LS.png');
