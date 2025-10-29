%% LS_v2.m (wrapper)
% This script is kept for backward compatibility. The full implementation
% now lives in adaptivehb.solvers.run_ls. Adjust the configuration JSON to
% change datasets or solver parameters.
%
% Optional: define CONFIG_PATH before running to override the default.
if ~exist('CONFIG_PATH', 'var') || isempty(CONFIG_PATH)
    CONFIG_PATH = fullfile('config', 'comparison.json');
end

config = adaptivehb.config.load_config(CONFIG_PATH);
results = adaptivehb.solvers.run_ls(config); %#ok<NASGU>
