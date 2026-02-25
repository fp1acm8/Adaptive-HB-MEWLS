function dataset = apply_noise(dataset, noiseCfg)
%APPLY_NOISE Apply synthetic perturbations to a dataset structure.
%   DATASET = APPLY_NOISE(DATASET, NOISECFG) takes the struct returned by
%   adaptivehb.io.load_dataset and injects reproducible perturbations as
%   described by the noise configuration. Supported noise types are:
%
%       - "none"     : return the dataset unchanged
%       - "gaussian" : add zero-mean Gaussian noise with standard deviation
%                      noiseCfg.standardDeviation. A fraction of samples can
%                      be inflated to simulate outliers.
%
%   This function operates on the dataset struct used by run_comparison and
%   is driven by the JSON noise configuration.  For a lower-level interface
%   that accepts raw (x, y, f) vectors and supports additional noise types
%   (e.g. "spike"), see adaptivehb.data.apply_noise.
%
%   The configuration fields are optional and default to sensible values
%   when omitted.
%
%   See also adaptivehb.io.load_dataset, adaptivehb.data.apply_noise.

% --- Input validation -------------------------------------------------
% Both arguments must be scalar structs.  noiseCfg defaults to empty so
% the caller can omit it (resulting in type="none", i.e. no noise).
arguments
    dataset (1, 1) struct
    noiseCfg (1, 1) struct = struct()
end

% --- Default configuration values -------------------------------------
% Any field missing from noiseCfg is filled with a safe default.
defaults = struct(...
    'type', "none", ...              % noise type: "none" or "gaussian"
    'standardDeviation', 0.0, ...    % sigma of the Gaussian perturbation
    'outlierFraction', 0.0, ...      % fraction of points marked as outliers
    'outlierInflation', 0.0, ...     % multiplicative inflation for outliers
    'seed', NaN);                    % RNG seed (NaN = no explicit seeding)

% Merge user-supplied fields with defaults (user values take precedence).
noiseCfg = fill_defaults(noiseCfg, defaults);

% --- Early exit if noise is disabled ----------------------------------
if strcmpi(noiseCfg.type, "none")
    return;  % dataset is returned unchanged
end

% --- Seed the random number generator for reproducibility -------------
% If a numeric seed is provided, fix the RNG state so that running the
% same config twice produces identical noise patterns.
if ~isnan(noiseCfg.seed)
    rng(noiseCfg.seed);
end

% --- Apply noise depending on the configured type ---------------------
switch lower(noiseCfg.type)
    case "gaussian"
        % Generate zero-mean Gaussian noise with standard deviation sigma.
        % Clamp sigma to eps to avoid a degenerate (zero-variance) case.
        sigma = max(noiseCfg.standardDeviation, eps);
        perturb = sigma * randn(size(dataset.f));

        % Add the noise to the observed values.
        f = dataset.f + perturb;

        % --- Outlier simulation ---
        % Select a random subset of k points and inflate their observed
        % values by a multiplicative factor (1 + |outlierInflation|).
        % This simulates gross errors / outlier measurements.
        n = numel(f);
        k = round(noiseCfg.outlierFraction * n);  % number of outliers

        if k > 0
            % Randomly pick k indices without replacement.
            idx = randperm(n, k);

            % Inflation factor: e.g. outlierInflation=0.2 -> multiply by 1.2.
            inflation = 1 + abs(noiseCfg.outlierInflation);
            f(idx) = f(idx) * inflation;

            % Update the outlier mask: reset all to false, then flag the
            % selected indices as outliers.
            dataset.is_outlier(:) = false;
            dataset.is_outlier(idx) = true;
        end

        % Store the perturbed observations back into the dataset struct.
        dataset.f = f;

    otherwise
        % Any unrecognised noise type triggers an informative error.
        error('adaptivehb:io:UnsupportedNoise', ...
            'Unsupported noise type "%s".', noiseCfg.type);
end

end


% =====================================================================
%  Local function: fill_defaults
%  Merges user-supplied config fields with default values.
% =====================================================================
function cfg = fill_defaults(cfg, defaults)
%FILL_DEFAULTS Copy missing or empty fields from defaults into cfg.
fields = fieldnames(defaults);
for i = 1:numel(fields)
    name = fields{i};
    % Only overwrite if the field is missing or empty in the user config.
    if ~isfield(cfg, name) || isempty(cfg.(name))
        cfg.(name) = defaults.(name);
    end
end
end
