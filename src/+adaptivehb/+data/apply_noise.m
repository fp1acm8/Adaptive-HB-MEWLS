function [augmented, metadata] = apply_noise(x, y, f, settings)
%APPLY_NOISE Injects artificial noise and outliers into a dataset.
%   [AUGMENTED, METADATA] = adaptivehb.data.apply_noise(X, Y, F, SETTINGS)
%   adds noise to the scalar field values F defined on the coordinates X
%   and Y. SETTINGS is a struct with the following optional fields:
%       addNoise          - logical flag to enable corruption (default: true)
%       outlierFraction   - fraction of samples to corrupt (default: 0.10)
%       outlierIntensity  - relative amplitude of the injected noise (default: 0.20)
%       seed              - RNG seed for reproducibility (default: 42)
%       noiseType         - 'gauss' or 'spike' (default: 'gauss')
%
%   The function returns AUGMENTED, a struct containing the vectors x, y,
%   f_true, f_noisy and is_outlier (logical mask), and METADATA which stores
%   bookkeeping information about the corruption that was applied.
%
%   The helper adaptivehb.data.format_noise_filename can be used together
%   with METADATA to create consistent filenames when persisting the data.
%
%   Note: for the struct-level interface used by run_comparison (driven by
%   JSON configuration), see adaptivehb.io.apply_noise instead.
%
%   See also adaptivehb.data.format_noise_filename,
%            adaptivehb.io.apply_noise.

    % --- Default settings when the caller omits the argument ----------
    if nargin < 4 || isempty(settings)
        settings = struct();
    end

    % Fill missing fields with defaults and sanitise value ranges.
    settings = local_normalise_settings(settings);

    % --- Force inputs to column vectors for consistent indexing -------
    x = x(:);
    y = y(:);
    f = f(:);

    % All three vectors must have the same number of elements.
    if ~isequal(numel(x), numel(y), numel(f))
        error('adaptivehb:data:apply_noise:DimensionMismatch', ...
              'X, Y and F must contain the same number of samples.');
    end

    % --- Initialise output variables ----------------------------------
    M = numel(f);                   % total number of samples
    rangeF = max(f) - min(f);       % range of the signal (used to scale noise)

    is_outlier = false(M, 1);      % logical mask: true for corrupted points
    f_noisy = f;                    % start with a copy of the clean signal
    outlierIdx = [];                % indices of corrupted samples
    noiseVector = zeros(0, 1);      % actual noise values added

    % --- Noise injection (only if enabled) ----------------------------
    if settings.addNoise
        % Fix the RNG state for reproducibility: same seed -> same noise.
        rng(settings.seed);

        % How many points to corrupt: round(fraction * N), clamped to [0, N].
        numOutliers = round(settings.outlierFraction * M);
        numOutliers = min(max(numOutliers, 0), M);

        if numOutliers > 0
            % Randomly select numOutliers indices without replacement.
            outlierIdx = randperm(M, numOutliers).';
            is_outlier(outlierIdx) = true;  % flag them in the mask

            % Generate the noise vector whose amplitude is proportional to
            % outlierIntensity * range(f).  This ensures the noise scale
            % is meaningful relative to the data.
            switch settings.noiseType
                case 'gauss'
                    % Gaussian noise: N(0, intensity*rangeF).
                    noiseVector = settings.outlierIntensity * rangeF * randn(numOutliers, 1);
                case 'spike'
                    % Uniform spike noise in [-intensity*rangeF, +intensity*rangeF].
                    noiseVector = settings.outlierIntensity * rangeF * (2 * rand(numOutliers, 1) - 1);
                otherwise
                    error('adaptivehb:data:apply_noise:UnsupportedNoiseType', ...
                          'Noise type "%s" is not supported.', settings.noiseType);
            end

            % Add the noise to the selected samples only.
            f_noisy(outlierIdx) = f_noisy(outlierIdx) + noiseVector;
        end
    end

    % --- Pack the augmented dataset into a struct ---------------------
    augmented = struct(...
        'x', x, ...               % x coordinates (unchanged)
        'y', y, ...               % y coordinates (unchanged)
        'f_true', f, ...          % original clean values
        'f_noisy', f_noisy, ...   % corrupted values (= f_true when no noise)
        'is_outlier', is_outlier  % logical mask of corrupted points
    );

    % --- Build metadata for traceability ------------------------------
    % This struct records everything needed to reproduce or describe the
    % noise that was applied, and is consumed by format_noise_filename.
    metadata = struct();
    metadata.noiseApplied = settings.addNoise;           % was noise enabled?
    metadata.numSamples = M;                             % total sample count
    metadata.rangeF = rangeF;                            % range of f (for reference)
    metadata.outlierIndices = outlierIdx;                 % which points were corrupted
    metadata.numOutliers = numel(outlierIdx);             % how many
    if M == 0
        metadata.actualOutlierFraction = 0;
    else
        metadata.actualOutlierFraction = metadata.numOutliers / M;  % actual fraction
    end
    metadata.noiseVector = noiseVector;                  % exact noise values added
    metadata.settings = settings;                        % full settings used
    metadata.noiseLabel = local_build_noise_label(settings);  % string label for filenames
end


% =====================================================================
%  Local function: local_normalise_settings
%  Fills missing fields with defaults and sanitises value ranges.
% =====================================================================
function settings = local_normalise_settings(settings)
    % Default values for all noise parameters.
    defaults = struct( ...
        'addNoise', true, ...           % enable noise injection
        'outlierFraction', 0.10, ...    % 10% of points corrupted
        'outlierIntensity', 0.20, ...   % noise amplitude = 20% of range(f)
        'seed', 42, ...                 % RNG seed for reproducibility
        'noiseType', 'gauss' ...        % Gaussian noise by default
    );

    % Merge: for each field, use the user's value if present and non-empty,
    % otherwise fall back to the default.
    fields = fieldnames(defaults);
    for iField = 1:numel(fields)
        name = fields{iField};
        if ~isfield(settings, name) || isempty(settings.(name))
            settings.(name) = defaults.(name);
        end
    end

    % --- Type coercion and range clamping -----------------------------
    settings.addNoise = logical(settings.addNoise);                 % ensure logical
    settings.outlierFraction = max(0, min(1, settings.outlierFraction));  % clamp to [0,1]
    settings.outlierIntensity = max(0, settings.outlierIntensity);  % non-negative
    settings.seed = double(settings.seed);                          % ensure double
    settings.noiseType = lower(string(settings.noiseType));         % normalise case
    settings.noiseType = char(settings.noiseType);                  % convert to char
end


% =====================================================================
%  Local function: local_build_noise_label
%  Generates a compact string label encoding the noise parameters.
%  Used by format_noise_filename for constructing filenames.
% =====================================================================
function label = local_build_noise_label(settings)
    % If noise was not applied, return an empty label (no suffix in filename).
    if ~settings.addNoise
        label = '';
        return;
    end

    % Format: "noise_gauss_frac0.10_int0.20_seed42"
    label = sprintf('noise_%s_frac%.2f_int%.2f_seed%d', ...
                    settings.noiseType, ...
                    settings.outlierFraction, ...
                    settings.outlierIntensity, ...
                    settings.seed);
end
