%% MEWLS_v3.m (wrapper)
% This script delegates to adaptivehb.solvers.run_mewls so legacy entry
% points continue to work while configuration is moved to JSON.
%
% Optional: define CONFIG_PATH before running to override the default.
if ~exist('CONFIG_PATH', 'var') || isempty(CONFIG_PATH)
    CONFIG_PATH = fullfile('config', 'comparison.json');
end

config = adaptivehb.config.load_config(CONFIG_PATH);
results = adaptivehb.solvers.run_mewls(config); %#ok<NASGU>
