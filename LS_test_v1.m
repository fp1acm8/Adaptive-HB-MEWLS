%% main_script_weighted_least_squares.m
% Script MATLAB per HB-splines + Weighted Least Squares con termine di penalizzazione
clearvars; close all; clc; format long;
DEBUG = true;  % Flag per output dettagliato

%% ---------------------- DEFINIZIONE DEL PROBLEMA -------------------------------
% Dati geografici e condizioni al contorno del problema di Laplace
problem_data.geo_name = 'geo_square.txt';
problem_data.nmnn_sides = [];
problem_data.drchlt_sides = [1 2 3 4];
problem_data.c_diff = @(x, y) ones(size(x));

% Definizione della soluzione esatta u(x,y) = exp(-C*r^2) con C=100
C = 100;
normax2 = @(x,y) ((x-.5).^2+(y-.5).^2);
problem_data.uex = @(x,y) exp(-C*normax2(x,y));
problem_data.f = @(x,y) 4*C*(1-C*normax2(x,y)).*problem_data.uex(x,y);
problem_data.g = @(x, y, ind) zeros(size(x));
problem_data.h = @(x, y, ind) problem_data.uex(x,y);
problem_data.graduex = @(x,y) -2*C*cat (1, ...
            reshape (problem_data.uex(x,y).*(x-.5), [1, size(x)]), ...
            reshape (problem_data.uex(x,y).*(y-.5), [1, size(x)]));

%% ------------------- PARAMETRI DELLO SPAZIO DISCRETO ------------------------
% Grado e regolarità degli HB-splines, griglia iniziale e regole di quadratura
method_data = struct(...
    'degree',      [2 2], ...     % B-spline biquadratic
    'regularity',  [1 1], ...     % Continuità C^1
    'nsub_coarse', [16 8], ...     % Suddivisione iniziale n×n elementi
    'nsub_refine', [2 2], ...     % Ogni raffinamento dimezza ogni spigolo
    'nquad',       [3 3], ...     % 3×3 punti di quadratura
    'space_type',  'standard', ...
    'truncated',   1);            % Truncated THB-splines

%% -------------------- PARAMETRI ADATTIVITÀ -------------------------
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

%% --------------------- CARICAMENTO DATI SPARSE -------------------------
% Selezione dataset
filename = '/MATLAB Drive/hierarchical_Bsplines_least_squares_with_GeoPDEs/data/example11_dataset.txt';

% Caricamento dati
fid = fopen(filename,'r');
if fid == -1
    error('File %s non trovato', filename);
end
data_3 = fscanf(fid,'%f %f %f\n',[3 inf]);
fclose(fid);

% Coordinate originali e valori f
data0 = data_3(1:2,:)';
f0 = data_3(3,:)';

% Mappa le coordinate in [0,1]^2
a = min(data0(:,1)); b = max(data0(:,1));
c = min(data0(:,2)); d = max(data0(:,2));
data = [(data0(:,1)-a)/(b-a), (data0(:,2)-c)/(d-c)];
f = f0;
M = size(data,1);

%% ---------- INIZIALIZZAZIONE SPAZIO GERARCHICO ----------
[hmsh, hspace, geometry] = adaptivity_initialize_laplace(problem_data, method_data);

%% -------- PARAMETRI WEIGHTED LEAST SQUARES -------------------------
lambda_pen = 1e-12;  % Parametro di penalizzazione
weight = ones(1,M) / M;  % Inizializzazione pesi uniformi

%% ================== CICLO ADATTATIVO ===========================
results = [];  % Matrice per salvare [livello, ndof, max_err, L2_err, cond_num]
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
    
    % Calcola errori
    error = abs(f - QI');
    max_err = max(error);
    L2_err = norm(error,2)/sqrt(length(error));  % RMSE
    
    % Salva risultati
    results = [results; lev, ndof, max_err, L2_err, cond_num];
    
    if DEBUG
        fprintf('cond(S)=%.2e | max_err=%.2e | L2_err=%.2e\n', ...
                cond_num, max_err, L2_err);
    end
    
    % Condizioni di arresto
    if lev > 1
        prev_err = results(end-1, 4);  % errore L2 precedente
        improvement = (prev_err - L2_err) / prev_err;
        if improvement < 0.01  % miglioramento < 1%
            fprintf('> Miglioramento insufficiente (%.3f%%) -> uscita\n', improvement*100);
            break;
        end
    end
    
    if max_err <= adaptivity_data.tol
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
    
    % Trova punti con errore superiore alla tolleranza
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
fprintf('\n%s\n', repmat('=', 1, 70));
fprintf('RISULTATI FINALI - WEIGHTED LEAST SQUARES\n');
fprintf('%s\n', repmat('=', 1, 70));
fprintf('Liv  ndof   max_err    L2_err     cond_num\n');
fprintf('%s\n', repmat('-', 1, 45));
fprintf('%3d %5d %10.4e %10.4e %10.2e\n', results');


fprintf('\n%s\n', repmat('=', 1, 50));
fprintf('ELABORAZIONE COMPLETATA\n');
fprintf('Dataset: %s\n', filename);
fprintf('Livelli finali: %d\n', size(results,1));
fprintf('DOF finali: %d\n', results(end,2));
fprintf('Errore finale: %.4e\n', results(end,4));
fprintf('%s\n', repmat('=', 1, 50));

%% -------------------- PLOT FINALI ------------------------------
% Mesh adattata
figure('Name', 'Mesh Adattata');
hmsh_plot_cells(hmsh, 10, 1);
view(0, 90); 
axis off; 
%title('Mesh finale adattata');

% Superficie della soluzione
figure('Name', 'Superficie WLS');
xplot = linspace(a, b, 151);
yplot = linspace(c, d, 151);
[Xplot, Yplot] = meshgrid(xplot, yplot);
try
    QIplot = sp_eval(QI_coeff', hspace_temp, geometry, [151 151]);
catch
    QIplot = sp_eval(QI_coeff', hspace, geometry, [151 151]);
end
QIplot = QIplot';
mesh(Xplot, Yplot, QIplot, 'FaceLighting', 'phong', 'FaceColor', [.5 .5 .5], ...
     'EdgeColor', 'none', 'AmbientStrength', 1.0, 'SpecularExponent', 15, ...
     'SpecularStrength', 0);
camlight(-37.5, 30+20);
title('Superficie Weighted Least Squares');
xlabel('X'); ylabel('Y'); zlabel('U');
axis off;

% Errore puntuale 3D
figure('Name', 'Errore 3D');
scatter3(data0(:,1), data0(:,2), error, 20, error, 'filled');
colorbar;
title('Errore |f-u| in 3D');
xlabel('X'); ylabel('Y'); zlabel('Errore');

% Dati originali 3D
figure('Name', 'Dati Originali');
scatter3(data0(:,1), data0(:,2), f0, 20, f0, 'filled');
colorbar;
%title('Dati originali');
xlabel('X'); ylabel('Y'); zlabel('f');

% Punti originali (proiezione sul piano xy)
figure('Name', 'Punti Dati');
plot(data0(:,1), data0(:,2), '.');
axis([a b c d]);
%title('Distribuzione punti dati');
xlabel('X'); ylabel('Y');

% Distribuzione dei pesi
figure('Name', 'Distribuzione Pesi');
histogram(weight, 50);
%title('Distribuzione finale dei pesi');
xlabel('Peso'); ylabel('Frequenza');

% Grafico di convergenza
if size(results, 1) > 1
    figure('Name', 'Convergenza');
    subplot(2,2,1);
    semilogy(results(:,1), results(:,3), 'o-');
    title('Errore Massimo vs Livello');
    xlabel('Livello'); ylabel('Errore Max');
    grid on;
    
    subplot(2,2,2);
    semilogy(results(:,1), results(:,4), 's-');
    title('Errore L2 vs Livello');
    xlabel('Livello'); ylabel('Errore L2');
    grid on;
    
    subplot(2,2,3);
    semilogy(results(:,2), results(:,3), 'o-');
    title('Errore Massimo vs DOF');
    xlabel('DOF'); ylabel('Errore Max');
    grid on;
    
    subplot(2,2,4);
    semilogy(results(:,2), results(:,4), 's-');
    title('Errore L2 vs DOF');
    xlabel('DOF'); ylabel('Errore L2');
    grid on;
end

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
%create_formatted_table_image(results, 'LS', 'black_forest_LS_pen1e-6.png');

