function result = mewls(params)
%MEWLS Adaptive HB-spline solver based on Maximum Entropy WLS.
%   RESULT = MEWLS(PARAMS) runs the entropy-based weighted solver. PARAMS
%   shares the same structure used by WEIGHTED_LEAST_SQUARES with solver
%   specific fields:
%       params.solver.lambda_smooth
%       params.solver.weight_floor
%       params.solver.me_options (struct with r, tol, max_iter, min_iter)

validate_params(params);
options = ensure_options(params);

% Dataset
config = struct('filename', params.dataset.filename);
if isfield(params.dataset, 'outlierFlag')
    config.outlierFlag = params.dataset.outlierFlag;
end
config.verbose = options.verbose;
dataset = adaptivehb.core.prepare_dataset(config);

% Hierarchical space
[hmsh, hspace, geometry] = adaptivity_initialize_laplace(params.problem, params.method);

state = struct( ...
    'hmsh', hmsh, ...
    'hspace', hspace, ...
    'geometry', geometry, ...
    'dataset', dataset, ...
    'solverParams', params.solver);

callbacks.solve = @(s) solve_step(s, options);
callbacks.refine = @(s, summary, adaptCfg) refine_step(s, summary, adaptCfg, options);

core = adaptivehb.core.run_adaptivity(state, callbacks, params.adaptivity, options);
result = assemble_result(core, dataset);
end

% =========================================================================
function validate_params(params)
required = {'dataset', 'problem', 'method', 'adaptivity', 'solver'};
for i = 1:numel(required)
    if ~isfield(params, required{i})
        error('mewls:MissingField', 'PARAMS must contain ''%s''.', required{i});
    end
end

if ~isfield(params.dataset, 'filename')
    error('mewls:MissingDataset', 'PARAMS.dataset must define filename.');
end

solverFields = {'lambda_smooth', 'weight_floor', 'me_options'};
for i = 1:numel(solverFields)
    if ~isfield(params.solver, solverFields{i})
        error('mewls:MissingSolverField', 'PARAMS.solver must define %s.', solverFields{i});
    end
end
end

function options = ensure_options(params)
if isfield(params, 'options') && isstruct(params.options)
    options = params.options;
else
    options = struct();
end

if ~isfield(options, 'verbose')
    options.verbose = false;
end
end

function step = solve_step(state, options)
solverParams = state.solverParams;
opts = solverParams.me_options;
opts.verbose = options.verbose;

[coefficients, weights, lambda2, condNumber, converged] = getcoeff_MEWLS( ...
    state.hspace, state.hmsh, state.dataset.points, state.dataset.fObserved, ...
    solverParams.lambda_smooth, opts, solverParams.weight_floor);

predictions = evaluate_spline_at_points(coefficients, state.hspace, state.dataset.points);

step = struct( ...
    'coefficients', coefficients, ...
    'weights', weights, ...
    'conditionNumber', condNumber, ...
    'predictions', predictions(:), ...
    'diagnostics', struct('lambda2', lambda2, 'converged', converged));

if options.verbose
    fprintf('    cond(S)=%.2e | lambda2=%.2e | converged=%d\n', condNumber, lambda2, converged);
end
end

function [state, refined] = refine_step(state, summary, adaptivityConfig, options)
errorMask = summary.pointwise > adaptivityConfig.tol;

if ~any(errorMask)
    if options.verbose
        fprintf('    > No supports selected for refinement.\n');
    end
    refined = false;
    return;
end

marked = support_containing_point(state.hspace, state.hmsh, state.dataset.points(errorMask, :));
oldNdof = state.hspace.ndof;
[state.hmsh, state.hspace] = adaptivity_refine(state.hmsh, state.hspace, marked, adaptivityConfig);
refined = state.hspace.ndof > oldNdof;

if options.verbose
    totalMarked = 0;
    for lv = 1:numel(marked)
        totalMarked = totalMarked + numel(marked{lv});
    end
    fprintf('    > Refined %d supports | new ndof %d\n', totalMarked, state.hspace.ndof);
end
end

function result = assemble_result(core, dataset)
result = struct();
result.dataset = dataset;
result.iterations = core.history;
result.terminationReason = core.terminationReason;
result.finalState = core.finalState;

if isempty(core.history)
    result.finalCoefficients = [];
    result.finalWeights = [];
    result.conditionNumbers = [];
    result.errorMetrics = struct([]);
    result.table = struct('headers', {{}}, 'data', []);
    return;
end

last = core.history(end);
result.finalCoefficients = last.coefficients;
result.finalWeights = last.weights;
result.conditionNumbers = [core.history.conditionNumber];
result.errorMetrics = arrayfun(@(h) h.errors.metrics, core.history);
[result.table.headers, result.table.data] = build_table(core.history);
end

function [headers, data] = build_table(history)
firstSummary = history(1).errors;
headers = [{'iter', 'level', 'ndof'}, firstSummary.headers, {'condition'}];
numRows = numel(history);
numCols = numel(headers);
data = zeros(numRows, numCols);
for i = 1:numRows
    row = [history(i).iteration, history(i).level, history(i).ndof, ...
        history(i).errors.values, history(i).conditionNumber];
    data(i, :) = row;
end
end

% =========================================================================
% Local MEWLS implementation (adapted from legacy script)
% =========================================================================
function [coeff_col, w, lambda2, condS, converged] = getcoeff_MEWLS(...
    hspace, hmsh, data, f, lambda_s, opts, w_floor)
M = size(data, 1);
N = hspace.ndof;

if N * M > 1e9
    error('getcoeff_MEWLS:TooLarge', 'Matrix too large: N=%d, M=%d', N, M);
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
        warning(ME.identifier, 'Failed to solve MEWLS system: %s', ME.message);
        converged = false;
        return;
    end

    pred = Phi * coeff_col;
    res2 = (f - pred).^2;

    [lambda2, newton_converged] = solve_entropy_constraint(res2, E2_target, lambda2, opts.tol, 500);
    if ~newton_converged && opts.verbose
        fprintf('      > Newton did not converge at iteration %d\n', iter);
    end

    w_new = exp(-lambda2 * res2);
    w_new = max(w_new, w_floor);
    w_new = w_new / sum(w_new);

    weight_change = norm(w_new - w, inf);
    lambda_change = abs(lambda2 - lambda2_old);

    if opts.verbose
        E2_current = sum(w_new .* res2);
        fprintf('      iter %2d: lambda2=%.2e | maxDw=%.2e | E2=%.2e\n', ...
            iter, lambda2, weight_change, E2_current);
    end

    w = w_new;

    if iter >= opts.min_iter && weight_change < opts.tol && lambda_change < opts.tol
        converged = true;
        break;
    end

    if any(~isfinite(w))
        warning('getcoeff_MEWLS:UnstableWeights', 'Weights became unstable.');
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

    lambda2 = max(lambda2 - g / gp, 1e-12);
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
