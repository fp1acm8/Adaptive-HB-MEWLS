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
%   fTrue, fNoisy and isOutlier (logical mask), and METADATA which stores
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

    if nargin < 4 || isempty(settings)
        settings = struct();
    end

    settings = local_normalise_settings(settings);

    x = x(:);
    y = y(:);
    f = f(:);

    if ~isequal(numel(x), numel(y), numel(f))
        error('adaptivehb:data:apply_noise:DimensionMismatch', ...
              'X, Y and F must contain the same number of samples.');
    end

    M = numel(f);
    rangeF = max(f) - min(f);

    isOutlier = false(M, 1);
    fNoisy = f;
    outlierIdx = [];
    noiseVector = zeros(0, 1);

    if settings.addNoise
        rng(settings.seed);

        numOutliers = round(settings.outlierFraction * M);
        numOutliers = min(max(numOutliers, 0), M);

        if numOutliers > 0
            outlierIdx = randperm(M, numOutliers).';
            isOutlier(outlierIdx) = true;

            switch settings.noiseType
                case 'gauss'
                    noiseVector = settings.outlierIntensity * rangeF * randn(numOutliers, 1);
                case 'spike'
                    noiseVector = settings.outlierIntensity * rangeF * (2 * rand(numOutliers, 1) - 1);
                otherwise
                    error('adaptivehb:data:apply_noise:UnsupportedNoiseType', ...
                          'Noise type "%s" is not supported.', settings.noiseType);
            end

            fNoisy(outlierIdx) = fNoisy(outlierIdx) + noiseVector;
        end
    end

    augmented = struct(
        'x', x,
        'y', y,
        'fTrue', f,
        'fNoisy', fNoisy,
        'isOutlier', isOutlier
    );

    metadata = struct();
    metadata.noiseApplied = settings.addNoise;
    metadata.numSamples = M;
    metadata.rangeF = rangeF;
    metadata.outlierIndices = outlierIdx;
    metadata.numOutliers = numel(outlierIdx);
    if M == 0
        metadata.actualOutlierFraction = 0;
    else
        metadata.actualOutlierFraction = metadata.numOutliers / M;
    end
    metadata.noiseVector = noiseVector;
    metadata.settings = settings;
    metadata.noiseLabel = local_build_noise_label(settings);
end

function settings = local_normalise_settings(settings)
    defaults = struct( ...
        'addNoise', true, ...
        'outlierFraction', 0.10, ...
        'outlierIntensity', 0.20, ...
        'seed', 42, ...
        'noiseType', 'gauss' ...
    );

    fields = fieldnames(defaults);
    for iField = 1:numel(fields)
        name = fields{iField};
        if ~isfield(settings, name) || isempty(settings.(name))
            settings.(name) = defaults.(name);
        end
    end

    settings.addNoise = logical(settings.addNoise);
    settings.outlierFraction = max(0, min(1, settings.outlierFraction));
    settings.outlierIntensity = max(0, settings.outlierIntensity);
    settings.seed = double(settings.seed);
    settings.noiseType = lower(string(settings.noiseType));
    settings.noiseType = char(settings.noiseType);
end

function label = local_build_noise_label(settings)
    if ~settings.addNoise
        label = '';
        return;
    end

    label = sprintf('noise_%s_frac%.2f_int%.2f_seed%d', ...
                    settings.noiseType, ...
                    settings.outlierFraction, ...
                    settings.outlierIntensity, ...
                    settings.seed);
end
