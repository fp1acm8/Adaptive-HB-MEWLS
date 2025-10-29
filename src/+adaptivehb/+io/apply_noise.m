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
%   The configuration fields are optional and default to sensible values
%   when omitted.

arguments
    dataset (1, 1) struct
    noiseCfg (1, 1) struct = struct()
end

defaults = struct(...
    'type', "none", ...
    'standardDeviation', 0.0, ...
    'outlierFraction', 0.0, ...
    'outlierInflation', 0.0, ...
    'seed', NaN);

noiseCfg = fill_defaults(noiseCfg, defaults);

if strcmpi(noiseCfg.type, "none")
    return;
end

if ~isnan(noiseCfg.seed)
    rng(noiseCfg.seed);
end

switch lower(noiseCfg.type)
    case "gaussian"
        sigma = max(noiseCfg.standardDeviation, eps);
        perturb = sigma * randn(size(dataset.fObserved));
        fObserved = dataset.fObserved + perturb;

        n = numel(fObserved);
        k = round(noiseCfg.outlierFraction * n);
        if k > 0
            % Inflate a subset of samples by a multiplicative factor.
            idx = randperm(n, k);
            inflation = 1 + abs(noiseCfg.outlierInflation);
            fObserved(idx) = fObserved(idx) * inflation;
            dataset.isOutlier(:) = false;
            dataset.isOutlier(idx) = true;
        end

        dataset.fObserved = fObserved;
    otherwise
        error('adaptivehb:io:UnsupportedNoise', ...
            'Unsupported noise type "%s".', noiseCfg.type);
end

end

function cfg = fill_defaults(cfg, defaults)
fields = fieldnames(defaults);
for i = 1:numel(fields)
    name = fields{i};
    if ~isfield(cfg, name) || isempty(cfg.(name))
        cfg.(name) = defaults.(name);
    end
end
end
