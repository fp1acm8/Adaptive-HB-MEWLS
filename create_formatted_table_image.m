function create_formatted_table_image(results, method_name, filename)
    fig = figure('Name', 'Tabella Risultati', 'Position', [100, 100, 800, 600]);
    axis off;
    
    % Titolo
    text(0.5, 0.95, ['RISULTATI FINALI - ' upper(method_name)], ...
         'FontSize', 16, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    
    % Riga separatrice
    line([0.1, 0.9], [0.90, 0.90], 'Color', 'k', 'LineWidth', 2);
    
    % Header della tabella
    headers = {'Liv', 'DOF', 'Max Err', 'L2 Err'};
    if size(results, 2) == 5
        headers = [headers, {'Cond Num'}];
    end
    
    x_positions = linspace(0.15, 0.85, length(headers));
    
    % Stampa header
    for i = 1:length(headers)
        text(x_positions(i), 0.85, headers{i}, ...
             'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    end
    
    % Riga separatrice header
    line([0.1, 0.9], [0.82, 0.82], 'Color', 'k', 'LineWidth', 1);
    
    % Dati
    y_start = 0.78;
    y_step = 0.06;
    
    for i = 1:size(results, 1)
        y_pos = y_start - (i-1) * y_step;
        
        % Livello
        text(x_positions(1), y_pos, sprintf('%d', results(i, 1)), ...
             'FontSize', 10, 'HorizontalAlignment', 'center');
        
        % DOF
        text(x_positions(2), y_pos, sprintf('%d', results(i, 2)), ...
             'FontSize', 10, 'HorizontalAlignment', 'center');
        
        % Max Error
        text(x_positions(3), y_pos, sprintf('%.2e', results(i, 3)), ...
             'FontSize', 10, 'HorizontalAlignment', 'center');
        
        % L2 Error
        text(x_positions(4), y_pos, sprintf('%.2e', results(i, 4)), ...
             'FontSize', 10, 'HorizontalAlignment', 'center');
        
        % Condition Number (se presente)
        if size(results, 2) == 5
            text(x_positions(5), y_pos, sprintf('%.1e', results(i, 5)), ...
                 'FontSize', 10, 'HorizontalAlignment', 'center');
        end
    end
    
    % Riga separatrice finale
    y_final = y_start - size(results, 1) * y_step;
    line([0.1, 0.9], [y_final, y_final], 'Color', 'k', 'LineWidth', 1);
    
    % Salva immagine
    if nargin < 3
        filename = 'tabella_risultati_formatted.png';
    end
    
    exportgraphics(fig, filename, 'Resolution', 300);
    fprintf('Tabella formattata salvata come: %s\n', filename);
end