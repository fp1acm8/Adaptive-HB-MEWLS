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

    % --- Default options ----------------------------------------------
    % If options is missing or empty, create an empty struct.
    if nargin < 3 || isempty(options)
        options = struct();
    end

    % Prefix: an optional string prepended before the basename.
    % Example: prefix = "outl_" -> "outl_glacier_noise_...txt"
    if ~isfield(options, 'Prefix') || isempty(options.Prefix)
        options.Prefix = '';
    end

    % Extension: the file extension including the dot.
    % Default is '.txt' (plain text, space-separated).
    if ~isfield(options, 'Extension') || isempty(options.Extension)
        options.Extension = '.txt';
    end

    % --- Input validation ---------------------------------------------
    % baseName must be a non-empty character vector or string.
    if ~(ischar(baseName) || isstring(baseName)) || strlength(string(baseName)) == 0
        error('adaptivehb:data:format_noise_filename:InvalidBaseName', ...
              'Basename must be a non-empty character vector or string.');
    end

    % metadata must be a struct (typically from adaptivehb.data.apply_noise).
    if nargin < 2 || ~isstruct(metadata)
        error('adaptivehb:data:format_noise_filename:InvalidMetadata', ...
              'Metadata must be the struct returned by apply_noise.');
    end

    % Ensure baseName is a char vector for consistent concatenation.
    baseName = char(baseName);

    % --- Build the noise suffix ---------------------------------------
    % If noise was applied, append a descriptive suffix encoding the noise
    % parameters. This makes filenames self-documenting.
    suffix = '';
    if isfield(metadata, 'noiseApplied') && metadata.noiseApplied
        % Prefer the pre-built noiseLabel if available in metadata.
        if isfield(metadata, 'noiseLabel') && ~isempty(metadata.noiseLabel)
            suffix = ['_' metadata.noiseLabel];
        % Otherwise, reconstruct the label from the settings sub-struct.
        elseif isfield(metadata, 'settings')
            suffix = ['_' local_build_noise_label(metadata.settings)];
        end
    end
    % If noise was NOT applied, suffix remains '' -> filename = "basename.txt".

    % --- Assemble the final filename ----------------------------------
    % Pattern: <prefix><basename><suffix><extension>
    % Example: "outl_glacier_noise_gauss_frac0.10_int0.20_seed42.txt"
    filename = sprintf('%s%s%s%s', options.Prefix, baseName, suffix, options.Extension);
end


% =====================================================================
%  Local function: local_build_noise_label
%  Fallback label construction when metadata.noiseLabel is not available.
% =====================================================================
function label = local_build_noise_label(settings)
    % If noise was not enabled, return empty (no suffix).
    if ~isfield(settings, 'addNoise') || ~settings.addNoise
        label = '';
        return;
    end
    % Format: "noise_gauss_frac0.10_int0.20_seed42"
    label = sprintf('noise_%s_frac%.2f_int%.2f_seed%d', ...
        settings.noiseType, settings.outlierFraction, settings.outlierIntensity, settings.seed);
end
