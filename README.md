# Adaptive-HB-MEWLS
MATLAB implementation of Hierarchical B-splines (HB-splines) with adaptive refinement driven by Maximum Entropy Weighted Least Squares (MEWLS).  
The project explores entropy-based weighting strategies for local approximation and adaptive mesh refinement, aiming to improve numerical stability and accuracy in spline-based geometric modeling and numerical analysis.

## Features
- Lightweight comparison harness for LS vs MEWLS-inspired polynomial
  solvers.
- Configurable experiments with JSON-driven dataset, noise, and solver
  definitions.
- Reusable utilities for loading datasets, injecting noise, and exporting
  rich visual summaries.
- Clean MATLAB package layout that makes it straightforward to introduce
  additional solver variants.

## Quick start

The repository is organised around a single comparison harness and a
handful of reusable utilities. A typical MATLAB session looks like:

```matlab
cd /path/to/Adaptive-HB-MEWLS
addpath('scripts');
run_comparison();
```

`run_comparison` automatically adds `src/` to the MATLAB path and, when
called without arguments, loads the default configuration found at
`config/comparison_config.json`.

To compare a different set of solvers or datasets, copy the JSON file,
adjust the fields you need, and pass the new path:

```matlab
run_comparison('config/my_experiment.json');
```

The helper [`scripts/generate_noisy_dataset.m`](scripts/generate_noisy_dataset.m)
is available for interactively corrupting an existing clean dataset using
the same noise utilities employed by the harness.

## Repository layout

- `config/` – JSON templates describing datasets, solver lists, and
  metrics.
- `data/` – Sample datasets used in the published experiments.
- `scripts/` – Entry points such as `run_comparison.m` and dataset
  preparation utilities.
- `src/` – MATLAB package (`+adaptivehb`) containing solvers, I/O helpers,
  and visualisation routines consumed by the harness.
- `reports/` – Created on demand; each run of the harness writes its
  outputs to a timestamped subfolder.

## Comparison harness

[`scripts/run_comparison.m`](scripts/run_comparison.m) orchestrates the
entire experiment:

1. loads a JSON configuration describing the dataset, optional noise
   perturbations, and solver definitions;
2. evaluates each solver with identical inputs;
3. aggregates the requested metrics (RMSE, max absolute error, MAE by
   default) and computes per-metric deltas against a reference solver;
4. emits figures and tabular summaries inside a timestamped directory
   under `reports/`.

All figures and tables are written to `reports/comparison_<timestamp>/` so
multiple experiments can coexist without manual bookkeeping.

### Generated artefacts

Each execution creates `reports/comparison_<timestamp>/` containing:

- `metrics.csv` / `metrics.json` / `metrics_table.mat`: the aggregated
  metrics table with delta columns showing the difference between each
  solver and the first solver listed in the configuration;
- `solution_surfaces.png`: ground truth vs solver predictions plotted
  side by side;
- `error_distributions.png`: residual histograms for each solver;
- `convergence_curves.png`: RMSE decay across the polynomial degrees
  explored by the simplified solvers;
- `dataset_snapshot.mat`: a copy of the dataset struct actually consumed
  during the run (after normalisation/noise injection).

### Configuration reference

The configuration file is a plain JSON document with four top-level
sections:

| Section     | Required keys | Notes |
|-------------|---------------|-------|
| `dataset`   | `path`        | Absolute or relative path to a text file with either three columns (clean data) or five columns (true/observed values plus an outlier flag). Optional `normalize` (default `true`) rescales coordinates to `[0, 1]^2`. |
| `noise`     | `type`        | Use `"none"` to keep the dataset untouched or `"gaussian"` together with `standardDeviation`, `outlierFraction`, `outlierInflation`, and `seed` to apply reproducible perturbations. |
| `solvers`   | –             | Array of solver specs, each with `name`, the fully qualified MATLAB function handle (for example `"adaptivehb.solvers.mewlsSolver"`), and an optional `parameters` struct passed through verbatim. |
| `metrics`   | –             | List of metric names to include in the output tables. The stock solvers expose `rmse`, `maxAbsError`, and `mae`. |

Any additional fields are preserved and handed to the respective
functions, making it easy to extend the harness with custom behaviour.

## Author
Francesco Pucci  
PhD Candidate — University of Florence / INdAM 
