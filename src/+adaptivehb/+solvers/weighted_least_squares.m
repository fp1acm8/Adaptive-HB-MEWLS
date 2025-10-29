function result = weighted_least_squares(params)
%WEIGHTED_LEAST_SQUARES Adaptive HB-spline solver using penalised WLS.
%   RESULT = WEIGHTED_LEAST_SQUARES(PARAMS) runs the adaptive solver using
%   the configuration contained in PARAMS. Mandatory fields:
%       params.dataset.filename
%       params.problem      (passed to adaptivity_initialize_laplace)
%       params.method       (space discretisation parameters)
%       params.adaptivity   (tolerance and limits)
%       params.solver.lambda_penalty
%
%   Optional:
%       params.dataset.outlierFlag
%       params.options.verbose
%
%   The returned RESULT struct exposes reusable artefacts (history,
%   coefficients, weights, condition numbers, formatted table data).

validate_params(params);
options = ensure_options(params);

% -------------------------------------------------------------------------
% Dataset preparation
% -------------------------------------------------------------------------
datasetConfig = struct('filename', params.dataset.filename);
if isfield(params.dataset, 'outlierFlag')
    datasetConfig.outlierFlag = params.dataset.outlierFlag;
end
if isfield(options, 'verbose')
    datasetConfig.verbose = options.verbose;
end

dataset = adaptivehb.core.prepare_dataset(datasetConfig);

% -------------------------------------------------------------------------
% Initial hierarchical space
% -------------------------------------------------------------------------
[hmsh, hspace, geometry] = adaptivity_initialize_laplace(params.problem, params.method);

state = struct( ...
    'hmsh', hmsh, ...
    'hspace', hspace, ...
    'geometry', geometry, ...
    'dataset', dataset, ...
    'weights', ones(1, dataset.info.numPoints) / dataset.info.numPoints, ...
    'solverParams', params.solver);

callbacks.solve = @(s) solve_step(s, options);
callbacks.refine = @(s, summary, adaptCfg) refine_step(s, summary, adaptCfg, options);

core = adaptivehb.core.run_adaptivity(state, callbacks, params.adaptivity, options);

result = assemble_result(core, dataset);
end

% =========================================================================
function validate_params(params)
requiredTop = {'dataset', 'problem', 'method', 'adaptivity', 'solver'};
for i = 1:numel(requiredTop)
    if ~isfield(params, requiredTop{i})
        error('weighted_least_squares:MissingField', ...
            'PARAMS must contain the field ''%s''.', requiredTop{i});
    end
end

if ~isfield(params.dataset, 'filename')
    error('weighted_least_squares:MissingDataset', ...
        'PARAMS.dataset must contain a filename field.');
end

if ~isfield(params.solver, 'lambda_penalty')
    error('weighted_least_squares:MissingSolverField', ...
        'PARAMS.solver must define lambda_penalty.');
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

function [step, state] = solve_step(state, options)
weights = state.weights;
solverParams = state.solverParams;

[coefficients, condNumber] = getcoeff_weighted_least_squares_pen( ...
    state.hspace, state.hmsh, state.dataset.points, state.dataset.fObserved, ...
    solverParams.lambda_penalty, weights);

predictions = sp_eval_alt(coefficients, state.hspace, state.dataset.points);

step = struct( ...
    'coefficients', coefficients, ...
    'weights', weights, ...
    'conditionNumber', condNumber, ...
    'predictions', predictions(:), ...
    'diagnostics', struct('lambda_penalty', solverParams.lambda_penalty));

if isfield(solverParams, 'update_weights') && isa(solverParams.update_weights, 'function_handle')
    newWeights = solverParams.update_weights(weights, step);
    if ~isempty(newWeights)
        state.weights = newWeights;
        step.weights = newWeights;
    end
end

if options.verbose
    fprintf('    cond(A)=%.2e\n', condNumber);
end
end

function [state, refined] = refine_step(state, summary, adaptivityConfig, options)
errorMask = summary.pointwise > adaptivityConfig.tol;

if ~any(errorMask)
    if options.verbose
        fprintf('    > No points above tolerance.\n');
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
    fprintf('    > Refined %d supports | new ndof %d\n', ...
        totalMarked, state.hspace.ndof);
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
