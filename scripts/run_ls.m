%% scripts/run_ls.m
% Thin wrapper around adaptivehb.solvers.run_ls using the default
% configuration file.
config = adaptivehb.config.load_config(fullfile('config', 'comparison.json'));
adaptivehb.solvers.run_ls(config);
