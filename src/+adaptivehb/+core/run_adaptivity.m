function result = run_adaptivity(state, callbacks, adaptivityConfig, options)
%RUN_ADAPTIVITY Core loop shared by adaptive HB-spline solvers.
%   RESULT = RUN_ADAPTIVITY(STATE, CALLBACKS, ADAPTIVITYCONFIG) executes the
%   adaptive refinement loop. STATE must contain the fields `hmsh`,
%   `hspace` and `dataset`. CALLBACKS is a struct with two function handles:
%       solve(state) -> stepStruct [,(updatedState)]
%       refine(state, summary, adaptivityConfig) -> [newState, refined]
%   ADAPTIVITYCONFIG stores stopping criteria (tol, max_ndof, max_level,
%   num_max_iter, min_relative_improvement).
%
%   OPTIONS.verbose enables textual output.
%
%   The returned RESULT struct contains:
%       history             Struct array with per-iteration diagnostics
%       terminationReason   Textual reason for exiting the loop
%       finalState          STATE after the last iteration
%
%   Each HISTORY entry holds fields: level, ndof, coefficients, weights,
%   conditionNumber, errors (see summarize_errors), diagnostics.

if nargin < 4
    options = struct();
end

if ~isfield(options, 'verbose')
    options.verbose = false;
end

requiredFields = {'hmsh', 'hspace', 'dataset'};
for k = 1:numel(requiredFields)
    if ~isfield(state, requiredFields{k})
        error('run_adaptivity:MissingField', ...
            'STATE must contain the field ''%s''.', requiredFields{k});
    end
end

if ~isfield(adaptivityConfig, 'tol')
    error('run_adaptivity:MissingTolerance', 'ADAPTIVITYCONFIG must define a tolerance (tol).');
end

maxNdof = get_field_with_default(adaptivityConfig, 'max_ndof', inf);
maxLevel = get_field_with_default(adaptivityConfig, 'max_level', inf);
maxIter = get_field_with_default(adaptivityConfig, 'num_max_iter', inf);
minImprovement = get_field_with_default(adaptivityConfig, 'min_relative_improvement', 0.0);

history = struct([]);
terminationReason = '';

iter = 0;
while iter < maxIter
    iter = iter + 1;
    level = state.hspace.nlevels;
    ndof = state.hspace.ndof;

    if options.verbose
        fprintf('\n=== Iteration %d | level %d | ndof %d ===\n', iter, level, ndof);
    end

    if nargout(callbacks.solve) >= 2
        [step, state] = callbacks.solve(state);
    else
        step = callbacks.solve(state);
    end

    if ~isstruct(step) || ~isfield(step, 'predictions')
        error('run_adaptivity:InvalidStep', 'solve callback must return a struct with a ''predictions'' field.');
    end

    summary = adaptivehb.core.summarize_errors(step.predictions, ...
        state.dataset.fTrue, state.dataset.isOutlier);

    entry = struct( ...
        'iteration', iter, ...
        'level', level, ...
        'ndof', ndof, ...
        'coefficients', get_field_with_default(step, 'coefficients', []), ...
        'weights', get_field_with_default(step, 'weights', []), ...
        'conditionNumber', get_field_with_default(step, 'conditionNumber', NaN), ...
        'errors', summary, ...
        'diagnostics', get_field_with_default(step, 'diagnostics', struct()));

    history = [history; entry]; %#ok<AGROW>

    if options.verbose
        fprintf('    max_err=%.3e | rmse=%.3e\n', ...
            summary.metrics.maxAll, summary.metrics.rmseAll);
    end

    % Check tolerance
    if summary.metrics.maxAll <= adaptivityConfig.tol
        terminationReason = 'tolerance';
        break;
    end

    % Check ndof and level limits
    if ndof >= maxNdof
        terminationReason = 'max_ndof';
        break;
    end

    if level >= maxLevel
        terminationReason = 'max_level';
        break;
    end

    % Check improvement (skip first iteration)
    if numel(history) > 1
        prevRmse = history(end-1).errors.metrics.rmseAll;
        currRmse = summary.metrics.rmseAll;
        if prevRmse > 0
            improvement = (prevRmse - currRmse) / prevRmse;
        else
            improvement = 0;
        end
        if improvement < minImprovement
            terminationReason = 'insufficient_improvement';
            break;
        end
    end

    % Perform refinement
    [newState, refined] = callbacks.refine(state, summary, adaptivityConfig);
    state = newState;

    if ~refined
        terminationReason = 'no_refinement';
        break;
    end
end

if isempty(terminationReason)
    if iter >= maxIter
        terminationReason = 'max_iterations';
    else
        terminationReason = 'completed';
    end
end

result = struct( ...
    'history', {history}, ...
    'terminationReason', terminationReason, ...
    'finalState', state);
end

function value = get_field_with_default(s, fieldName, defaultValue)
if isstruct(s) && isfield(s, fieldName)
    value = s.(fieldName);
else
    value = defaultValue;
end
end
