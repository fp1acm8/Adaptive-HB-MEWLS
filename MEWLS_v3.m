%% main_script_MEWLS.m
% Script MATLAB per HB-splines + MEWLS 
clearvars; close all; clc; format long;
DEBUG = true;  % Flag per output dettagliato

%% ================== PARAMETRI MODIFICABILI ==================
% File dataset (modifica qui per cambiare dataset)
filename = '/MATLAB Drive/hierarchical_Bsplines_least_squares_with_GeoPDEs/data/example13_dataset_noise_gauss_frac0.05_int0.10_seed18.txt';

% Input per tipo dataset
outlier_dataset = input('Dataset con outliers? (0=no (3 colonne), 1=si (5 colonne)): ');

% Parametri del problema (modifica se necessario)
problem_data = struct(...
    'geo_name',     'geo_square.txt', ...       % File con nodo e connettività
    'nmnn_sides',   [], ...                      % Nessuna condizione Neumann
    'drchlt_sides', [1 2 3 4]);                  % Tutti i lati Dirichlet

C = 100;  % Costante per soluzione esatta
r2 = @(x,y) (x-0.5).^2 + (y-0.5).^2;
problem_data.uex = @(x,y) exp(-C * r2(x,y));    % Soluzione esatta
problem_data.c_diff = @(x,y) ones(size(x));      % Coeff. diffusione uniforme
problem_data.f = @(x,y) 4*C*(1 - C*r2(x,y)) .* problem_data.uex(x,y);  % Forcing function
problem_data.g = @(x,y,ind) zeros(size(x));      % Neumann nullo
problem_data.h = @(x,y,ind) problem_data.uex(x,y);  % Valore Dirichlet

% Parametri dello spazio discreto
method_data = struct(...
    'degree',      [2 2], ...     % B-spline biquadratic
    'regularity',  [1 1], ...     % Continuità C^1
    'nsub_coarse', [16 8], ...     % Suddivisione iniziale mesh
    'nsub_refine', [2 2], ...     % Ogni raffinamento dimezza ogni spigolo
    'nquad',       [3 3], ...     % 3×3 punti di quadratura
    'space_type',  'standard', ...
    'truncated',   1);            % Truncated THB-splines

% Parametri adattività
adaptivity_data = struct(...
    'flag',           'functions', ...   % Marcatura basata su valori di funzione
    'mark_param',     0.5, ...           % Non usato nella strategia entropy-based
    'mark_strategy', 'MS', ...   
    'max_level',      6, ...             % Profondità massima dell'albero
    'max_ndof',       50000, ...         % Numero max di dof consentiti
    'num_max_iter',   11, ...            % Iterazioni massime di adattività
    'max_nel',        50000, ...         % Numero max di elementi
    'tol',            5e-3, ...          % Tolleranza globale sull'errore
    'adm',            0);

% Parametri MEWLS & smoothing
lambda_smooth = 1e-12;  % Stabilizzazione Tikhonov (ridotto)
w_floor       = 1e-12; % Peso minimo per flooring (ridotto)
opts_me       = struct(...
    'r',        2, ...   % Fattore di riduzione dell'errore target
    'tol',      1e-4, ...  % Tolleranza per convergenza
    'max_iter', 100, ...   % Iterazioni massime MEWLS
    'min_iter', 3, ...     % Iterazioni minime per garantire stabilità
    'verbose',  DEBUG);

%% --------------------- CARICAMENTO DATI SPARSE -------------------------
% Caricamento dataset
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
a = min(data0); 
b = max(data0);
data = [(data0(:,1)-a(1))/(b(1)-a(1)), (data0(:,2)-a(2))/(b(2)-a(2))];
f = f_noise;  % Usa f_noise per fitting
M = size(data,1);

fprintf('Dataset caricato: %s (formato: %d colonne)\n', filename, expected_cols);
fprintf('Numero di punti: %d\n', M);
if expected_cols == 5
    fprintf('Outliers marcati: %d (%.2f%%)\n', sum(is_outlier), 100*sum(is_outlier)/M);
end

%% ---------- INIZIALIZZAZIONE SPAZIO GERARCHICO ----------
[hmsh,hspace,geometry] = adaptivity_initialize_laplace(problem_data,method_data);

%% ================== CICLO ADATTATIVO ===========================
results = [];  % Matrice per salvare [livello, ndof, max_err_all, L2_err_all, max_err_no_out, L2_err_no_out]
while true
    % Estrai livello e grado di libertà correnti
    lev  = hspace.nlevels;
    ndof = hspace.ndof;
    
    if DEBUG
        fprintf('\n=== Livello %d | ndof = %d ===\n', lev, ndof); 
    end
    
    % Calcola coefficienti e pesi con MEWLS (entropy-based)
    [coeff_col, w, lambda2, condS, converged] = getcoeff_MEWLS(...
        hspace, hmsh, data, f, lambda_smooth, opts_me, w_floor);
    
    if ~converged
        warning('MEWLS non è convergente al livello %d', lev);
    end
    
    % Valuta la spline sui punti dati e calcola errore punto a punto rispetto a f_true
    vals = evaluate_spline_at_points(coeff_col, hspace, data);
    err  = abs(f_true - vals);
    
    % Metriche globali di errore (su tutti i punti)
    max_err_all = max(err);                % Errore massimo su tutti
    L2_err_all  = sqrt(mean(err.^2));      % Errore L2 su tutti
    
    % Metriche escludendo outliers (is_outlier == 1)
    no_out_mask = (is_outlier == 0);
    if any(no_out_mask)
        err_no_out = err(no_out_mask);
        max_err_no_out = max(err_no_out);
        L2_err_no_out  = sqrt(mean(err_no_out.^2));
    else
        max_err_no_out = NaN;
        L2_err_no_out = NaN;
    end
    
    results = [results; lev, ndof, max_err_all, L2_err_all, max_err_no_out, L2_err_no_out];
    
    if DEBUG
        fprintf('cond(S)=%.2e | max_err_all=%.2e | L2_err_all=%.2e | conv=%d\n', ...
                condS, max_err_all, L2_err_all, converged);
        fprintf('max_err_no_out=%.2e | L2_err_no_out=%.2e\n', max_err_no_out, L2_err_no_out);
        fprintf('Peso min=%.2e | Peso max=%.2e | λ₂=%.2e\n', ...
                min(w), max(w), lambda2);
    end
    
    % Condizioni di arresto (usa max_err_all per la tolleranza)
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
    
    % Seleziona i punti da raffinare: error > tol (su err da f_true)
    mark_pt = err > adaptivity_data.tol;
    ex_idx  = find(mark_pt);
    
    if DEBUG
        fprintf('→ marcati %d punti su %d (%.2f%%) con err > %.2f\n', ...
                numel(ex_idx), M, 100*numel(ex_idx)/M, adaptivity_data.tol);
    end
    
    if isempty(ex_idx)
        disp('> Nessun supporto da raffinare -> uscita'); 
        break;
    end
    
    % Trova i supporti corrispondenti e raffina
    marked = support_containing_point(hspace, hmsh, data(ex_idx,:));
    [hmsh,hspace] = adaptivity_refine(hmsh, hspace, marked, adaptivity_data);
    
    if DEBUG
        fprintf('→ raffinati %d supporti | nuovo ndof = %d\n', ...
                numel(marked{lev}), hspace.ndof);
    end
end

%% -------------------- RISULTATI --------------------------------
% Stampa tabella comparativa: livello | ndof | errori
fprintf('\n%s\n', repmat('=', 1, 80));
fprintf('RISULTATI FINALI\n');
fprintf('%s\n', repmat('=', 1, 80));

if outlier_dataset == 1
    fprintf('Liv  DOF   max_err_all    L2_err_all   max_err_no_out    L2_err_no_out\n');
    fprintf('%s\n', repmat('-', 1, 70));
    fprintf('%3d %5d %14.4e %14.4e %14.4e %14.4e\n', results');
else
    fprintf('Liv  DOF   Max Err    L2 Err\n');
    fprintf('%s\n', repmat('-', 1, 35));
    fprintf('%3d %5d %10.4e %10.4e\n', results(:, [1,2,3,4])');
end

% % -------------------- PLOT FINALI ------------------------------
% is_square = abs(b(1)-a(1)) == abs(b(2)-a(2));
% outlier_mask = (is_outlier == 1);
% 
% figure('Name', 'Mesh');
% hmsh_plot_cells(hmsh, 10, 1);
% view(0, 90);
% if is_square, axis equal; end
% axis off;

% figure('Name', 'Surface MEWLS');
% [XX, YY] = meshgrid(linspace(a(1), b(1), 101), linspace(a(2), b(2), 101));
% Ugrid = sp_eval(coeff_col(:), hspace, geometry, [101 101]);
% surf(XX, YY, Ugrid, 'EdgeColor', 'none'); 
% camlight;
% if is_square, axis equal; end
% xlabel('X'); ylabel('Y'); zlabel('U');
% 
% % ----------------- 3D Pesi w -----------------
% figure('Name', 'Weights 3D');
% scatter3(data0(:,1), data0(:,2), w, 20, w, 'filled');
% hold on;
% if outlier_dataset == 1 && any(outlier_mask)
%     scatter3(data0(outlier_mask,1), data0(outlier_mask,2), w(outlier_mask), 30, 'r', 'filled', 'MarkerEdgeColor', 'k');
% end
% colorbar;
% xlabel('X'); ylabel('Y'); zlabel('Peso');
% 
% % --------------- Dati originali 3D ---------------
% figure('Name', 'Original Data');
% scatter3(data0(:,1), data0(:,2), f, 20, f, 'filled');
% hold on;
% if outlier_dataset == 1 && any(outlier_mask)
%     scatter3(data0(outlier_mask,1), data0(outlier_mask,2), f(outlier_mask), 30, 'r', 'filled', 'MarkerEdgeColor', 'k');
% end
% colorbar;
% xlabel('X'); ylabel('Y'); zlabel('f');
% 
% % Distribuzione dei pesi
% figure('Name', 'Weight Distribution');
% histogram(w, 50);
% title('Final MEWLS weight distribution');
% xlabel('Weight'); ylabel('Frequency');
% 
% % --------------- Grafico di convergenza ---------------
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


%% ================= FUNZIONI LOCALI CORRETTE =============================
function [coeff_col, w, lambda2, condS, converged] = getcoeff_MEWLS(...
    hspace, hmsh, data, f, lambda_s, opts, w_floor)
% Implementazione del MEWLS per B-spline gerarchiche
% 
% Input:
%   hspace, hmsh: strutture di spazio gerarchico
%   data: punti sparsi [M x 2]
%   f: valori target [M x 1]
%   lambda_s: parametro di smoothing Tikhonov
%   opts: parametri MEWLS
%   w_floor: peso minimo per flooring
%
% Output:
%   coeff_col: coefficienti della spline [N x 1]
%   w: pesi finali [M x 1]
%   lambda2: parametro Lagrange finale
%   condS: numero di condizione della matrice del sistema
%   converged: flag di convergenza
    M = size(data, 1);            % numero punti dati
    N = hspace.ndof;              % numero dof correnti
    
    % Controllo dimensioni
    if N * M > 1e9
        error('Matrice troppo grande: N=%d, M=%d', N, M);
    end
    
    % Costruzione matrice di valutazione Phi [M x N]
    % Ogni riga i contiene i valori delle N funzioni di base nel punto data(i,:)
    Phi = zeros(M, N);
    for i = 1:M
        % Valuta tutte le funzioni di base nel punto data(i,:)
        basis_vals = sp_eval_alt(speye(N), hspace, data(i,:));
        Phi(i, :) = basis_vals(:)';
    end
    
    % Matrice di regolarizzazione (operatore Laplaciano)
    if lambda_s > 0
        Mass = op_gradgradu_gradgradv_hier(hspace, hspace, hmsh);
    else
        Mass = sparse(N, N);
    end
    
    % --- Inizializzazione con pesi uniformi ---
    w = ones(M, 1) / M;              % pesi uniformi normalizzati
    
    % Sistema lineare iniziale: (Phi^T * W * Phi + lambda_s * Mass) * coeff = Phi^T * W * f
    W = spdiags(w, 0, M, M);
    S = Phi' * W * Phi + lambda_s * Mass;
    rhs = Phi' * (w .* f);
    
    % Controlla condizionamento
    condS = condest(S);
    if condS > 1e10
        warning('Sistema mal condizionato: cond(S) = %.2e', condS);
    end
    
    % Risolvi sistema iniziale
    coeff_col = S \ rhs;
    
    % Calcola previsioni e residui iniziali
    pred = Phi * coeff_col;
    res2 = (f - pred).^2;
    
    % Target E2: errore medio pesato iniziale diviso per fattore di riduzione
    E2_target = sum(w .* res2) / opts.r;
    
    % Inizializza parametri iterativi
    lambda2 = 1e-6;
    converged = false;
    
    % --- Ciclo iterativo MEWLS ---
    for iter = 1:opts.max_iter
        % Salva stato precedente
        w_old = w;
        lambda2_old = lambda2;
        
        % Risolvi sistema corrente
        W = spdiags(w, 0, M, M);
        S = Phi' * W * Phi + lambda_s * Mass;
        rhs = Phi' * (w .* f);
        
        try
            coeff_col = S \ rhs;
        catch ME
            warning(ME.identifier,'Errore nella risoluzione del sistema: %s', ME.message);
            converged = false;
            return;
        end
        
        % Calcola previsioni e residui
        pred = Phi * coeff_col;
        res2 = (f - pred).^2;
        
        % Trova lambda2 usando Newton per soddisfare il vincolo di entropia
        [lambda2, newton_converged] = solve_entropy_constraint(...
            res2, E2_target, lambda2, opts.tol, 500);
        
        if ~newton_converged
            if opts.verbose
                fprintf('   Warning: Newton non convergente all''iterazione %d\n', iter);
            end
        end
        
        % Aggiorna pesi secondo la formula di entropia
        w_new = exp(-lambda2 * res2);
        
        % Applica flooring per stabilità numerica
        w_new = max(w_new, w_floor);
        
        % Normalizza pesi
        w_new = w_new / sum(w_new);
        
        % Calcola metriche di convergenza
        weight_change = norm(w_new - w, inf);
        lambda_change = abs(lambda2 - lambda2_old);
        
        if opts.verbose
            E2_current = sum(w_new .* res2);
            fprintf('   iter %2d: λ₂=%.2e, maxΔw=%.2e, E2=%.2e\n', ...
                    iter, lambda2, weight_change, E2_current);
        end
        
        % Aggiorna pesi
        w = w_new;
        
        % Test di convergenza
        if iter >= opts.min_iter
            if weight_change < opts.tol && lambda_change < opts.tol
                converged = true;
                if opts.verbose
                    fprintf('   Convergenza raggiunta all''iterazione %d\n', iter);
                end
                break;
            end
        end
        
        % Controllo stabilità
        if any(isnan(w)) || any(isinf(w))
            warning('Pesi numericamente instabili all''iterazione %d', iter);
            w = w_old;
            lambda2 = lambda2_old;
            break;
        end
    end
    
    % Risolvi sistema finale
    W = spdiags(w, 0, M, M);
    S = Phi' * W * Phi + lambda_s * Mass;
    rhs = Phi' * (w .* f);
    coeff_col = S \ rhs;
    
    % Aggiorna numero di condizione finale
    condS = condest(S);
    
    if ~converged && opts.verbose
        fprintf('   MEWLS non convergente dopo %d iterazioni\n', opts.max_iter);
    end
end
function [lambda2, converged] = solve_entropy_constraint(res2, E2_target, lambda2_init, tol, max_iter)
% Risolve il vincolo di entropia usando il metodo di Newton
% Trova lambda2 tale che: sum(res2 .* exp(-lambda2 * res2)) = E2_target * sum(exp(-lambda2 * res2))
    lambda2 = max(lambda2_init, 1e-12);  % Assicura che lambda2 sia positivo
    converged = false;
    
    for k = 1:max_iter
        % Calcola esponenziali (con controllo per overflow)
        exp_terms = exp(-lambda2 * res2);
        
        % Verifica overflow/underflow
        if any(~isfinite(exp_terms))
            lambda2 = lambda2 * 0.5;  % Riduce lambda2 se c'è overflow
            continue;
        end
        
        % Funzione obiettivo: g(lambda2) = sum(res2 .* exp_terms) - E2_target * sum(exp_terms)
        sum_weighted = sum(res2 .* exp_terms);
        sum_exp = sum(exp_terms);
        g = sum_weighted - E2_target * sum_exp;
        
        % Test di convergenza
        if abs(g) < tol
            converged = true;
            break;
        end
        
        % Derivata: g'(lambda2) = -sum(res2.^2 .* exp_terms) + E2_target * sum(res2 .* exp_terms)
        gp = -sum(res2.^2 .* exp_terms) + E2_target * sum_weighted;
        
        % Controllo derivata
        if abs(gp) < 1e-15
            break;  % Derivata troppo piccola
        end
        
        % Aggiornamento Newton con controllo per mantenere lambda2 positivo
        lambda2_new = lambda2 - g / gp;
        lambda2 = max(lambda2_new, 1e-12);
        
        % Controllo per step troppo grandi
        if lambda2 > 1e6
            lambda2 = 1e6;
        end
    end
end

function vals = evaluate_spline_at_points(coeff, hspace, points)
% Valuta la spline definita dai coefficienti nei punti specificati
% 
% Input:
%   coeff: coefficienti della spline [N x 1]
%   hspace: struttura dello spazio gerarchico
%   points: punti di valutazione [M x 2]
%
% Output:
%   vals: valori della spline nei punti [M x 1]
    M = size(points, 1);
    vals = zeros(M, 1);
    
    % Valuta la spline in ogni punto
    for i = 1:M
        % Valuta tutte le funzioni di base nel punto points(i,:)
        basis_vals = sp_eval_alt(speye(length(coeff)), hspace, points(i,:));
        vals(i) = coeff' * basis_vals(:);
    end
end
%% -------------------- RISULTATI FINALI ------------------------------
%create_formatted_table_image(results, 'MEWLS', 'black_forest_MEWLS.png');
