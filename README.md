# Adaptive-HB-MEWLS
MATLAB implementation of Hierarchical B-splines (HB-splines) with adaptive refinement driven by Maximum Entropy Weighted Least Squares (MEWLS).  
The project explores entropy-based weighting strategies for local approximation and adaptive mesh refinement, aiming to improve numerical stability and accuracy in spline-based geometric modeling and numerical analysis.

## Features
- Construction of HB-splines basis with hierarchical refinement
- Weighted Least Squares fitting using Maximum Entropy weights
- Adaptive refinement guided by local error and weight distribution
- Visualization tools for datasets, hierarchical meshes, and fitting results
- Modular structure for future extensions (THB-splines, multi-dimensional fitting, etc.)

## Comparison harness

The repository now includes a lightweight harness for comparing solver
variants under common conditions. The entry point is
[`scripts/run_comparison.m`](scripts/run_comparison.m), which

1. loads a JSON configuration describing the dataset, optional noise
   perturbations, and solver definitions;
2. evaluates each solver with identical inputs;
3. aggregates the requested metrics (RMSE, max absolute error, MAE by
   default) and computes per-metric deltas against a reference solver;
4. emits figures and tabular summaries inside a timestamped directory
   under `reports/`.

Run the default experiment from the repository root with:

```matlab
addpath('scripts');
run_comparison();
```

Both commands rely on the default configuration stored in
`config/comparison_config.json`. To launch a custom experiment, copy the
file, adjust the fields, and pass the new path to `run_comparison`.

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

### Customising datasets and noise scenarios

The JSON configuration exposes three key sections:

- `dataset`: provide a path to any supported dataset (three-column clean
  files or five-column files with explicit noise/outlier annotations) and
  opt into coordinate normalisation via `"normalize": true`.
- `noise`: select `"type": "none"` to keep the dataset untouched or use
  `"gaussian"` with additional parameters (`standardDeviation`,
  `outlierFraction`, `outlierInflation`, `seed`) to generate reproducible
  perturbations.
- `solvers`: list solver entries with a user-facing `name`, the fully
  qualified MATLAB function to call, and optional parameter structs. New
  solvers can be integrated by implementing a function that accepts the
  dataset struct and a parameter struct and returns a result struct with
  `name`, `prediction`, `metrics`, and `convergence` fields.

Additions or modifications to these sections are automatically picked up
by the harness—no changes to `run_comparison.m` are required.

## Author
Francesco Pucci  
PhD Candidate — University of Florence / INdAM 
