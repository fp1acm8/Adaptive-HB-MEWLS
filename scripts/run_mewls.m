%% scripts/run_mewls.m
% Thin wrapper around adaptivehb.solvers.run_mewls using the default
% configuration file.
config = adaptivehb.config.load_config(fullfile('config', 'comparison.json'));
adaptivehb.solvers.run_mewls(config);
