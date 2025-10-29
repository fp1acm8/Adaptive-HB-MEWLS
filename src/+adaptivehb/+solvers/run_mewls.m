function results = run_mewls(config)
%RUN_MEWLS Execute the MEWLS solver using a configuration struct.
%   RESULTS = RUN_MEWLS(CONFIG) performs the multi-entropy weighted least
%   squares solve driven by CONFIG. The configuration struct is typically
%   created through adaptivehb.config.load_config.

arguments
    config (1, 1) struct
end

mewls_defaults = struct(...
    'lambda_smooth', 1e-12, ...
    'w_floor', 1e-12, ...
    'opts_me', struct('r', 2, 'tol', 1e-4, 'max_iter', 100, ...
                      'min_iter', 3, 'verbose', false));

if isfield(config, 'solvers') && isfield(config.solvers, 'mewls')
    mewls_config = merge_structs(mewls_defaults, config.solvers.mewls);
else
    mewls_config = mewls_defaults;
end

if ~isfield(mewls_config, 'opts_me') || ~isstruct(mewls_config.opts_me)
    mewls_config.opts_me = mewls_defaults.opts_me;
else
    mewls_config.opts_me = merge_structs(mewls_defaults.opts_me, mewls_config.opts_me);
end

if ~isfield(config, 'dataset')
    error('adaptivehb:solvers:MissingDataset', ...
        'The configuration does not define a dataset section.');
end

dataset = adaptivehb.solvers.load_dataset(config.dataset);
M = dataset.num_points;
data = dataset.data;
f_true = dataset.f_true;
f = dataset.f_noise;
is_outlier = dataset.is_outlier;

fprintf('Dataset caricato: %s (formato: %d colonne)\n', ...
    dataset.relative_path, dataset.expected_cols);
fprintf('Numero di punti: %d\n', M);
if dataset.expected_cols == 5
    fprintf('Outliers marcati: %d (%.2f%%)\n', sum(is_outlier), ...
        100 * sum(is_outlier) / M);
end

problem_data = build_problem_data();
method_data = build_method_data(config);
adaptivity_data = build_adaptivity_data(config);

[hmsh, hspace, geometry] = adaptivity_initialize_laplace(problem_data, method_data);

results_matrix = [];
tol_satisfied = false;
space_stuck = false;
DEBUG = logical(mewls_config.opts_me.verbose);

while true
    lev = hspace.nlevels;
    ndof = hspace.ndof;

    if DEBUG
        fprintf('\n=== Livello %d | ndof = %d ===\n', lev, ndof);
    end

    [coeff_col, w, lambda2, condS, converged] = getcoeff_MEWLS(...
        hspace, hmsh, data, f, mewls_config.lambda_smooth, ...
        mewls_config.opts_me, mewls_config.w_floor);

    if ~converged
        warning('MEWLS non è convergente al livello %d', lev);
    end

    vals = evaluate_spline_at_points(coeff_col, hspace, data);
    err = abs(f_true - vals);

    max_err_all = max(err);
    L2_err_all = sqrt(mean(err.^2));

    no_out_mask = (is_outlier == 0);
    if any(no_out_mask)
        err_no_out = err(no_out_mask);
        max_err_no_out = max(err_no_out);
        L2_err_no_out = sqrt(mean(err_no_out.^2));
    else
        max_err_no_out = NaN;
        L2_err_no_out = NaN;
    end

    results_matrix = [results_matrix; lev, ndof, max_err_all, L2_err_all, ...
        max_err_no_out, L2_err_no_out];

    if DEBUG
        fprintf('cond(S)=%.2e | max_err_all=%.2e | L2_err_all=%.2e | conv=%d\n', ...
            condS, max_err_all, L2_err_all, converged);
        fprintf('max_err_no_out=%.2e | L2_err_no_out=%.2e\n', ...
            max_err_no_out, L2_err_no_out);
        fprintf('Peso min=%.2e | Peso max=%.2e | λ₂=%.2e\n', ...
            min(w), max(w), lambda2);
    end

    if size(results_matrix, 1) > 1
        prev_err = results_matrix(end-1, 4);
        improvement = (prev_err - L2_err_all) / prev_err;
        if improvement < 0
            fprintf('> Miglioramento insufficiente (%.3f%%) -> uscita\n', ...
                improvement * 100);
            break;
        end
    end

    if max_err_all <= adaptivity_data.tol
        disp('> Tolleranza globale soddisfatta');
        tol_satisfied = true;
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

    ex_idx = find(err > adaptivity_data.tol);
    if DEBUG
        fprintf('-> marcati %d punti su %d (%.2f%%) con err > %.2f\n', ...
            numel(ex_idx), M, 100 * numel(ex_idx) / M, adaptivity_data.tol);
    end

    if isempty(ex_idx)
        disp('> Nessun supporto da raffinare -> uscita');
        break;
    end

    marked = support_containing_point(hspace, hmsh, data(ex_idx, :));
    hmsh_prev = hmsh;
    hspace_prev = hspace;

    [hmsh, hspace] = adaptivity_refine(hmsh, hspace, marked, adaptivity_data);

    if DEBUG
        fprintf('-> raffinati %d supporti | nuovo ndof = %d\n', ...
            numel(marked{lev}), hspace.ndof);
    end

    if hspace.ndof == hspace_prev.ndof
        space_stuck = true;
        disp('> Spazio invariato -> uscita');
        break;
    end
end

print_summary(dataset, results_matrix, adaptivity_data, tol_satisfied, space_stuck);

if exist('w', 'var')
    weights = w;
else
    weights = [];
end

if exist('coeff_col', 'var')
    coeffs = coeff_col;
else
    coeffs = [];
end

results = struct('history', results_matrix, ...
                 'tol_satisfied', tol_satisfied, ...
                 'space_stuck', space_stuck, ...
                 'weights', weights, ...
                 'coefficients', coeffs, ...
                 'geometry', geometry, ...
                 'hspace', hspace, ...
                 'hmsh', hmsh);
end

function [coeff_col, w, lambda2, condS, converged] = getcoeff_MEWLS(...
    hspace, hmsh, data, f, lambda_s, opts, w_floor)
%GETCOEFF_MEWLS Run the MEWLS iteration for the current hierarchy level.
    M = size(data, 1);
    N = hspace.ndof;

    if N * M > 1e9
        error('Matrice troppo grande: N=%d, M=%d', N, M);
    end

    Phi = zeros(M, N);
    for i = 1:M
        basis_vals = sp_eval_alt(speye(N), hspace, data(i, :));
        Phi(i, :) = basis_vals(:)';
    end

    if lambda_s > 0
        Mass = op_gradgradu_gradgradv_hier(hspace, hspace, hmsh);
    else
        Mass = sparse(N, N);
    end

    w = ones(M, 1) / M;
    W = spdiags(w, 0, M, M);
    S = Phi' * W * Phi + lambda_s * Mass;
    rhs = Phi' * (w .* f);

    condS = condest(S);
    if condS > 1e10
        warning('Sistema mal condizionato: cond(S) = %.2e', condS);
    end

    coeff_col = S \ rhs;
    pred = Phi * coeff_col;
    res2 = (f - pred).^2;

    E2_target = sum(w .* res2) / opts.r;

    lambda2 = 1e-6;
    converged = false;

    for iter = 1:opts.max_iter
        w_old = w;
        lambda2_old = lambda2;

        W = spdiags(w, 0, M, M);
        S = Phi' * W * Phi + lambda_s * Mass;
        rhs = Phi' * (w .* f);

        try
            coeff_col = S \ rhs;
        catch ME
            warning(ME.identifier, 'Errore nella risoluzione del sistema: %s', ME.message);
            converged = false;
            return;
        end

        pred = Phi * coeff_col;
        res2 = (f - pred).^2;

        [lambda2, newton_converged] = solve_entropy_constraint(...
            res2, E2_target, lambda2, opts.tol, 500);

        if ~newton_converged && opts.verbose
            fprintf('   Warning: Newton non convergente all''iterazione %d\n', iter);
        end

        w_new = exp(-lambda2 * res2);
        w_new = max(w_new, w_floor);
        w_new = w_new / sum(w_new);

        weight_change = norm(w_new - w, inf);
        lambda_change = abs(lambda2 - lambda2_old);

        if opts.verbose
            E2_current = sum(w_new .* res2);
            fprintf('   iter %2d: λ₂=%.2e, maxΔw=%.2e, E2=%.2e\n', ...
                iter, lambda2, weight_change, E2_current);
        end

        w = w_new;

        if iter >= opts.min_iter
            if weight_change < opts.tol && lambda_change < opts.tol
                converged = true;
                if opts.verbose
                    fprintf('   Convergenza raggiunta all''iterazione %d\n', iter);
                end
                break;
            end
        end

        if any(isnan(w)) || any(isinf(w))
            warning('Pesi numericamente instabili all''iterazione %d', iter);
            w = w_old;
            lambda2 = lambda2_old;
            break;
        end
    end

    W = spdiags(w, 0, M, M);
    S = Phi' * W * Phi + lambda_s * Mass;
    rhs = Phi' * (w .* f);
    coeff_col = S \ rhs;

    condS = condest(S);

    if ~converged && opts.verbose
        fprintf('   MEWLS non convergente dopo %d iterazioni\n', opts.max_iter);
    end
end

function [lambda2, converged] = solve_entropy_constraint(res2, E2_target, lambda2_init, tol, max_iter)
    lambda2 = max(lambda2_init, 1e-12);
    converged = false;

    for k = 1:max_iter
        exp_terms = exp(-lambda2 * res2);

        if any(~isfinite(exp_terms))
            lambda2 = lambda2 * 0.5;
            continue;
        end

        sum_weighted = sum(res2 .* exp_terms);
        sum_exp = sum(exp_terms);
        g = sum_weighted - E2_target * sum_exp;

        if abs(g) < tol
            converged = true;
            break;
        end

        gp = -sum(res2.^2 .* exp_terms) + E2_target * sum_weighted;

        if abs(gp) < 1e-15
            break;
        end

        lambda2_new = lambda2 - g / gp;
        lambda2 = max(lambda2_new, 1e-12);

        if lambda2 > 1e6
            lambda2 = 1e6;
        end
    end
end

function vals = evaluate_spline_at_points(coeff, hspace, points)
    M = size(points, 1);
    vals = zeros(M, 1);

    for i = 1:M
        basis_vals = sp_eval_alt(speye(length(coeff)), hspace, points(i, :));
        vals(i) = coeff' * basis_vals(:);
    end
end

function print_summary(dataset, results_matrix, adaptivity_data, tol_satisfied, space_stuck)
    if isempty(results_matrix)
        warning('Nessun risultato disponibile per la stampa del sommario.');
        return;
    end

    fprintf('\n%s\n', repmat('=', 1, 80));
    fprintf('RISULTATI FINALI - MEWLS\n');
    fprintf('%s\n', repmat('=', 1, 80));

    if dataset.expected_cols == 5
        fprintf('Liv  DOF   max_err_all    L2_err_all   max_err_no_out    L2_err_no_out\n');
        fprintf('%s\n', repmat('-', 1, 70));
        fprintf('%3d %5d %14.4e %14.4e %14.4e %14.4e\n', results_matrix');
    else
        fprintf('Liv  DOF   Max Err    L2 Err\n');
        fprintf('%s\n', repmat('-', 1, 35));
        fprintf('%3d %5d %10.4e %10.4e\n', results_matrix(:, [1, 2, 3, 4])');
    end

    fprintf('\n%s\n', repmat('=', 1, 50));
    fprintf('ELABORAZIONE COMPLETATA\n');
    fprintf('Dataset: %s\n', dataset.relative_path);
    fprintf('Livelli finali: %d\n', size(results_matrix, 1));
    fprintf('DOF finali: %d\n', results_matrix(end, 2));
    fprintf('Errore finale (all): %.4e\n', results_matrix(end, 4));
    fprintf('Tolleranza soddisfatta: %d | Spazio invariato: %d\n', ...
        tol_satisfied, space_stuck);
    fprintf('Tol target: %.2e\n', adaptivity_data.tol);
    fprintf('%s\n', repmat('=', 1, 50));
end

function problem_data = build_problem_data()
    C = 100;
    r2 = @(x, y) (x - 0.5).^2 + (y - 0.5).^2;

    problem_data = struct(...
        'geo_name',    'geo_square.txt', ...
        'nmnn_sides',  [], ...
        'drchlt_sides', [1 2 3 4]);

    problem_data.uex = @(x, y) exp(-C * r2(x, y));
    problem_data.c_diff = @(x, y) ones(size(x));
    problem_data.f = @(x, y) 4 * C * (1 - C * r2(x, y)) .* problem_data.uex(x, y);
    problem_data.g = @(x, y, ind) zeros(size(x));
    problem_data.h = @(x, y, ind) problem_data.uex(x, y);
end

function method_data = build_method_data(config)
    defaults = struct(...
        'degree',      [2 2], ...
        'regularity',  [1 1], ...
        'nsub_coarse', [16 8], ...
        'nsub_refine', [2 2], ...
        'nquad',       [3 3], ...
        'space_type',  'standard', ...
        'truncated',   1);

    if isfield(config, 'method')
        method_data = merge_structs(defaults, config.method);
    else
        method_data = defaults;
    end

    method_data.degree = reshape(method_data.degree, 1, []);
    method_data.regularity = reshape(method_data.regularity, 1, []);
    method_data.nsub_coarse = reshape(method_data.nsub_coarse, 1, []);
    method_data.nsub_refine = reshape(method_data.nsub_refine, 1, []);
    method_data.nquad = reshape(method_data.nquad, 1, []);
    method_data.space_type = char(method_data.space_type);
end

function adaptivity_data = build_adaptivity_data(config)
    defaults = struct(...
        'flag', 'functions', ...
        'mark_param', 0.5, ...
        'mark_strategy', 'MS', ...
        'max_level', 6, ...
        'max_ndof', 50000, ...
        'num_max_iter', 11, ...
        'max_nel', 50000, ...
        'tol', 5e-3, ...
        'adm', 0);

    if isfield(config, 'adaptivity')
        adaptivity_data = merge_structs(defaults, config.adaptivity);
    else
        adaptivity_data = defaults;
    end

    adaptivity_data.flag = char(adaptivity_data.flag);
    adaptivity_data.mark_strategy = char(adaptivity_data.mark_strategy);
end

function merged = merge_structs(defaults, overrides)
    merged = defaults;
    fields = fieldnames(overrides);
    for idx = 1:numel(fields)
        merged.(fields{idx}) = overrides.(fields{idx});
    end
end
