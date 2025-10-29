function filename = format_noise_filename(baseName, metadata, options)
%FORMAT_NOISE_FILENAME Build a filename that encodes the noise settings.
%   NAME = adaptivehb.data.format_noise_filename(BASENAME, METADATA) returns
%   a filename in the style BASENAME_noise_<type>_frac##_int##_seed##.txt when
%   METADATA.noiseApplied is true. When no noise was applied, NAME reduces to
%   BASENAME.txt.
%
%   Additional options can be provided in the struct OPTIONS:
%       Prefix     - prepended to the basename (default: '')
%       Extension  - file extension, including the dot (default: '.txt')
%
%   METADATA is the struct returned by adaptivehb.data.apply_noise.

    if nargin < 3 || isempty(options)
        options = struct();
    end
    if ~isfield(options, 'Prefix') || isempty(options.Prefix)
        options.Prefix = '';
    end
    if ~isfield(options, 'Extension') || isempty(options.Extension)
        options.Extension = '.txt';
    end

    if ~(ischar(baseName) || isstring(baseName)) || strlength(string(baseName)) == 0
        error('adaptivehb:data:format_noise_filename:InvalidBaseName', ...
              'Basename must be a non-empty character vector or string.');
    end

    if nargin < 2 || ~isstruct(metadata)
        error('adaptivehb:data:format_noise_filename:InvalidMetadata', ...
              'Metadata must be the struct returned by apply_noise.');
    end

    baseName = char(baseName);

    suffix = '';
    if isfield(metadata, 'noiseApplied') && metadata.noiseApplied
        if isfield(metadata, 'noiseLabel') && ~isempty(metadata.noiseLabel)
            suffix = ['_' metadata.noiseLabel];
        elseif isfield(metadata, 'settings')
            suffix = ['_' local_build_noise_label(metadata.settings)];
        end
    end

    filename = sprintf('%s%s%s%s', options.Prefix, baseName, suffix, options.Extension);
end

function label = local_build_noise_label(settings)
    if ~isfield(settings, 'addNoise') || ~settings.addNoise
        label = '';
        return;
    end
    label = sprintf('noise_%s_frac%.2f_int%.2f_seed%d', ...
        settings.noiseType, settings.outlierFraction, settings.outlierIntensity, settings.seed);
end
