function plot_txt_3d(filename)
% PLOT_TXT_3D  Carica e visualizza in 3D un dataset Nx3 da file .txt
%    PLOT_TXT_3D          -> apre un selettore file
%    PLOT_TXT_3D('file.txt') -> carica direttamente il file
%
% Il file deve avere 3 colonne numeriche: X Y Z (separatore: spazio o tab).
% Gestisce righe vuote, header testuali e virgole decimali.

    %=== Se non passato, chiedi il file ===%
    if nargin < 1 || isempty(filename)
        [f, p] = uigetfile({'*.txt;*.csv','Text/CSV files (*.txt, *.csv)'; '*.*','Tutti i file'}, ...
                            'Seleziona il file con 3 colonne (X Y Z)');
        if isequal(f,0), disp('Operazione annullata.'); return; end
        filename = fullfile(p, f);
    end

    %=== Leggi il file come testo grezzo per gestire casi "difficili" ===%
    raw = fileread(filename);

    % Sostituisci eventuali virgole decimali con punti (es. "3,14" -> "3.14")
    % Evita però di toccare le virgole che fungono da separatore CSV puro:
    % se il file è chiaramente CSV con virgole come separatore tra colonne,
    % salta questa sostituzione globale.
    looksCSV = contains(raw, newline) && contains(raw, ',') && ~contains(raw, ';');
    if ~looksCSV
        raw = regexprep(raw, '(?<=[0-9]),(?=[0-9])', '.');
    end

    % Rimuovi eventuali righe vuote iniziali/finali
    raw = strtrim(raw);

    %=== Prova il parsing robusto con textscan ===%
    % Accetta spazi, tab o virgole come separatore; salta header non numerici.
    fmt = '%f%f%f';
    data = [];
    try
        % Usa un "fopen virtuale" con stringa
        tmpFile = [tempname,'.txt'];
        fidTmp = fopen(tmpFile,'w');
        fwrite(fidTmp, raw);
        fclose(fidTmp);

        fid = fopen(tmpFile,'r');
        % Salta righe che NON iniziano con cifre, segno, o punto (header)
        cleanedLines = strings(0,1);
        while true
            tline = fgetl(fid);
            if ~ischar(tline), break; end
            if ~isempty(regexp(strtrim(tline), '^[\+\-\.0-9]', 'once'))
                cleanedLines(end+1,1) = string(tline); %#ok<AGROW>
            end
        end
        fclose(fid);

        % Reimpacchetta e rileggi con separatori multipli (spazio/tab/virgola/;).
        if isempty(cleanedLines)
            error('Nessuna riga numerica trovata.');
        end
        cleaned = join(cleanedLines, newline);
        cleaned = regexprep(cleaned, '[;\t,]+', ' ');  % normalizza separatori in spazio singolo
        cleaned = regexprep(cleaned, '\s+', ' ');

        % Parse finale
        C = textscan(cleaned, fmt, 'CollectOutput', true);
        data = C{1};

        delete(tmpFile);
    catch ME
        warning('Parsing con textscan fallito (%s). Provo readmatrix...', ME.message);
        try
            data = readmatrix(filename, 'FileType', 'text');
        catch
            error('Impossibile leggere il file: %s', ME.message);
        end
    end

    %=== Validazioni base ===%
    if isempty(data) || size(data,2) < 3
        error('Il file deve contenere almeno 3 colonne numeriche (X Y Z).');
    end
    if size(data,2) > 3
        % Tieni solo le prime 3 colonne
        data = data(:,1:3);
    end

    % Rimuovi righe con NaN
    data = data(all(isfinite(data),2), :);

    if isempty(data)
        error('Dopo la pulizia non restano righe numeriche valide.');
    end

    X = data(:,1); Y = data(:,2); Z = data(:,3);

    %=== Plot 3D a punti, colorati per Z ===%
    figure('Name', sprintf('3D: %s', filename), 'Color','w');
    s = scatter3(X, Y, Z, 16, Z, 'filled'); %#ok<NASGU>
    grid on; box on;
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title(sprintf('Punti 3D (%d campioni): %s', size(data,1), shortName(filename)), 'Interpreter','none');
    axis vis3d; view(45, 30);
    colorbar; c = colormap; %#ok<NASGU> % usa colormap di default
    rotate3d on;

    %=== Opzionale: prova una superficie se i punti sono densi ===%
    try
        if size(data,1) >= 50
            % Delaunay triangulation per avere un "tappeto" indicativo
            DT = delaunayTriangulation(X, Y);
            tri = DT.ConnectivityList;
            figure('Name','Superficie (triangolazione Delaunay)','Color','w');
            trisurf(tri, X, Y, Z, Z, 'EdgeColor','none', 'FaceAlpha', 0.9);
            grid on; box on;
            xlabel('X'); ylabel('Y'); zlabel('Z');
            title('Superficie stimata (Delaunay)'); colorbar;
            axis tight vis3d; view(45,30);
        end
    catch
        % Se fallisce (punti degeneri, duplicati, ecc.), ignora silenziosamente
    end

    fprintf('Caricati %d punti da: %s\n', size(data,1), filename);
end

function s = shortName(p)
% Restituisce solo il nome file + estensione
    [~, name, ext] = fileparts(p);
    s = [name, ext];
end
