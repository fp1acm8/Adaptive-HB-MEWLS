function config = load_config(config_path)
%LOAD_CONFIG Read a JSON configuration file for Adaptive HB solvers.
%   CONFIG = LOAD_CONFIG(CONFIG_PATH) parses the JSON file located at
%   CONFIG_PATH and returns a MATLAB struct containing configuration
%   information. Relative paths inside the JSON are expanded to absolute
%   paths using the repository root. The repository root is assumed to be
%   the parent directory of the folder containing the configuration file.
%
%   The JSON schema is intentionally permissive; missing sections fall back
%   to sensible defaults where possible.
%
%   Example:
%       config = adaptivehb.config.load_config('config/comparison.json');
%
%   See also JSONDECODE, FULLFILE.

arguments
    config_path (1, :) char
end

if ~isfile(config_path)
    error('adaptivehb:config:MissingFile', ...
        'Configuration file not found: %s', config_path);
end

raw_text = fileread(config_path);
config = jsondecode(raw_text);

% Determine repository root (parent of configuration directory)
config_dir = fileparts(config_path);
repo_root = fileparts(config_dir);

config.paths.repo_root = repo_root;
config.paths.config_dir = config_dir;

% Normalise dataset paths
if isfield(config, 'dataset')
    config.dataset = normalise_dataset(config.dataset, repo_root);
end

% Normalise output directories
if isfield(config, 'output')
    config.output = normalise_output(config.output, repo_root);
end

end

function dataset = normalise_dataset(dataset, repo_root)
    if ~isfield(dataset, 'filename')
        error('adaptivehb:config:MissingDataset', ...
            'The configuration must include dataset.filename');
    end

    if isfield(dataset, 'root') && ~isempty(dataset.root)
        dataset.root_relative = char(dataset.root);
        dataset.root = fullfile(repo_root, dataset.root_relative);
    else
        dataset.root_relative = '';
        dataset.root = repo_root;
    end

    dataset.filename = char(dataset.filename);
    dataset.relative_path = fullfile(dataset.root_relative, dataset.filename);
    dataset.full_path = fullfile(repo_root, dataset.relative_path);
end

function output = normalise_output(output, repo_root)
    output_fields = fieldnames(output);
    for idx = 1:numel(output_fields)
        field = output_fields{idx};
        value = output.(field);
        if ischar(value) || (isstring(value) && isscalar(value))
            value_char = char(value);
            output.([field '_relative']) = value_char;
            absolute_path = fullfile(repo_root, value_char);
            output.(field) = absolute_path;
            if ~isfolder(absolute_path)
                mkdir(absolute_path);
            end
        end
    end
end
