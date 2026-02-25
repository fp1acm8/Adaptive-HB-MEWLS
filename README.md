# Adaptive-HB-MEWLS

MATLAB framework for comparing polynomial surface-fitting solvers — in particular Ordinary Least Squares (OLS) versus Maximum Entropy Weighted Least Squares (MEWLS) — on 2-D point-cloud datasets with optional synthetic noise and outlier injection.

The project explores entropy-based weighting strategies for local approximation and adaptive mesh refinement, aiming to improve numerical stability and accuracy in spline-based geometric modelling and numerical analysis. The polynomial component implemented here is designed to be integrated with the adaptive hierarchical B-spline (HBS) code based on GeoPDEs.

## Requirements

| Requirement | Details |
|-------------|---------|
| **MATLAB**  | R2020b or later (required for `arguments` blocks and `readmatrix`) |
| **Toolboxes** | None — the polynomial solver framework uses only base MATLAB functions |
| **OS** | Any platform supported by MATLAB (Windows, macOS, Linux) |
| **GeoPDEs** | Required for the hierarchical B-spline component (HBS). Not needed for the standalone polynomial solver. Official site: https://rafavzqz.github.io/geopdes/ |

## GeoPDEs Dependency (HBS component)

The complete extension of the framework to adaptive hierarchical B-splines (the HBS code, with adaptive mesh refinement) requires the **GeoPDEs** library installed on the MATLAB path. The GeoPDEs functions used in the HBS code are:

- `adaptivity_initialize_laplace` — initialises `hmsh`/`hspace` on the square domain (`geo_square.txt`)
- `adaptivity_refine` — refines the hierarchical mesh and B-spline space
- `hspace_subdivision_matrix` — subdivision matrix for hierarchical level change
- `op_gradgradu_gradgradv_hier` — assembles the penalisation matrix (biharmonic)
- `sp_eval` — evaluates the spline on a regular grid
- `hmsh_plot_cells` — plots the hierarchical mesh cells
- `sp_get_cells` — returns the elements in the support of a B-spline basis function

The corresponding HBS code files are: `getcoeff_weighted_least_squares_pen.m`, `sp_eval_alt.m`, `basisfun_multi.m`, `support_containing_point.m`.

GeoPDEs is available at: **https://rafavzqz.github.io/geopdes/**

## Quick start

```matlab
cd /path/to/Adaptive-HB-MEWLS

% Option A: run the full comparison pipeline (generates reports/)
addpath('scripts');
run_comparison();

% Option B: see a quick metrics printout without file output
run_quick_demo();

% Option C: verify everything works
run_tests();
```

`run_comparison` automatically adds `src/` to the MATLAB path and loads the default configuration at `config/comparison_config.json`.

## Repository layout

```
Adaptive-HB-MEWLS/
├── config/
│   └── comparison_config.json      % Default experiment configuration
├── data/                           % Sample datasets (clean and noisy)
│   ├── glacier.txt
│   ├── black_forest.txt
│   ├── example11_dataset.txt
│   ├── example13_dataset.txt
│   └── ...                         % Pre-generated noisy variants
├── scripts/
│   ├── run_comparison.m            % Main orchestrator (entry point)
│   ├── run_quick_demo.m            % Non-interactive demo
│   ├── run_tests.m                 % Verification suite
│   └── generate_noisy_dataset.m    % Interactive noise generation wizard
├── src/+adaptivehb/                % MATLAB package
│   ├── +solvers/
│   │   ├── leastSquaresSolver.m    % OLS polynomial fitter
│   │   ├── mewlsSolver.m          % Entropy-weighted WLS fitter (MEWLS)
│   │   ├── polynomialDesignMatrix.m% 2-D polynomial design matrix builder
│   │   └── compute_metrics.m      % Shared RMSE / MaxAbs / MAE calculator
│   ├── +io/
│   │   ├── load_dataset.m         % Dataset loader (3- or 5-column files)
│   │   └── apply_noise.m          % Struct-level noise injector (JSON-driven)
│   ├── +data/
│   │   ├── apply_noise.m          % Vector-level noise injector (gauss/spike)
│   │   └── format_noise_filename.m% Filename builder encoding noise params
│   └── +viz/
│       ├── plot_solution_surface.m % 3-D scatter comparison
│       ├── plot_error_distribution.m% Error histograms
│       └── plot_convergence_curve.m % RMSE vs polynomial degree
├── reports/                        % Auto-generated (one subfolder per run)
├── .gitignore
├── LICENSE                         % Apache 2.0
└── README.md
```

## How it works

The comparison pipeline follows four steps:

```
 ┌──────────────┐     ┌──────────────┐     ┌──────────────┐     ┌──────────────┐
 │ 1. Load data │ ──> │ 2. Add noise │ ──> │ 3. Run       │ ──> │ 4. Aggregate │
 │    from disk │     │ (optional)   │     │    solvers   │     │    & report  │
 └──────────────┘     └──────────────┘     └──────────────┘     └──────────────┘
```

1. **Load** — `adaptivehb.io.load_dataset` reads a text file (3 or 5 columns) and optionally normalises coordinates to [0,1]². Data are in the form `(u,v,f(u,v))`.
2. **Noise** — `adaptivehb.io.apply_noise` injects Gaussian noise and outliers according to the JSON configuration.
3. **Solve** — Each solver listed in the config is called with identical inputs. The built-in solvers sweep polynomial degrees 1 to `degree`.
4. **Report** — Metrics (RMSE, Max Absolute Error, MAE) are aggregated into CSV/JSON/MAT tables; diagnostic plots are saved as PNGs.

## Configuration reference

The JSON configuration has four sections:

### `dataset`

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `path` | string | *required* | Path to the text file (relative to repo root) |
| `normalize` | boolean | `true` | Rescale coordinates to [0,1]² |

### `noise`

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `type` | string | `"none"` | `"none"` or `"gaussian"` |
| `standardDeviation` | number | `0.0` | Std. dev. of the Gaussian perturbation |
| `outlierFraction` | number | `0.0` | Fraction of samples to mark as outliers |
| `outlierInflation` | number | `0.0` | Multiplicative inflation factor for outlier values |
| `seed` | integer | `NaN` | RNG seed for reproducibility (`NaN` = no seeding) |

### `solvers`

An array of solver specifications:

| Field | Type | Description |
|-------|------|-------------|
| `name` | string | Display name used in output tables and legends |
| `function` | string | Fully qualified MATLAB function (e.g. `"adaptivehb.solvers.mewlsSolver"`) |
| `parameters` | object | Passed verbatim to the solver as `method_data` struct |

### `metrics`

An array of metric names to include in the output: `"rmse"`, `"maxAbsError"`, `"mae"`.

## Dataset format

| Columns | Format | Description |
|---------|--------|-------------|
| 3 | `x  y  f` | Clean dataset — observed value equals true value |
| 5 | `x  y  f_true  f  is_outlier` | Noisy dataset with outlier flags (0/1) |

All columns are space-separated. The loader auto-detects the format.

## Generated output

Each call to `run_comparison` creates `reports/comparison_<timestamp>/` containing:

| File | Description |
|------|-------------|
| `metrics.csv` | Aggregated metrics table with delta columns |
| `metrics.json` | Same data in JSON format |
| `metrics_table.mat` | MATLAB table saved for programmatic use |
| `solution_surfaces.png` | Ground truth vs solver QI approximations (3-D scatter) |
| `error_distributions.png` | Per-solver error histograms (PDF normalised) |
| `convergence_curves.png` | RMSE decay over polynomial degree |
| `dataset_snapshot.mat` | Copy of the dataset struct used during the run |

## How to add a custom solver

1. Create a function file in `src/+adaptivehb/+solvers/` (or anywhere on the MATLAB path):

```matlab
function result = mySolver(dataset, method_data)
    % dataset.data       — N-by-2 coordinates (normalised to [0,1]^2)
    % dataset.f          — N-by-1 observed values
    % dataset.f_true     — N-by-1 true values (for metrics)
    % dataset.is_outlier — N-by-1 logical mask
    % method_data        — struct with your custom parameters

    % ... your fitting logic ...

    QI = ...;  % N-by-1 approximation

    m = adaptivehb.solvers.compute_metrics(dataset.f_true, QI);
    result.name = "MySolver";
    result.QI = QI;
    result.QI_coeff = [];
    result.metrics = struct('rmse', m(1), 'maxAbsError', m(2), 'mae', m(3));
    result.convergence = table();  % optional
end
```

2. Add it to the configuration JSON:

```json
{
  "name": "MySolver",
  "function": "adaptivehb.solvers.mySolver",
  "parameters": { "degree": 3, "myParam": 42 }
}
```

3. Run `run_comparison('config/my_config.json')`.

## Mathematical background

### Ordinary Least Squares (OLS)

Given N data points (x_i, y_i, f_i) and a polynomial design matrix **A** (Phi) whose columns are the monomials x^p · y^q with p + q ≤ d, OLS finds the coefficient vector **QI_coeff** that minimises:

```
||A·QI_coeff - f||²
```

Solution: **QI_coeff** = **A** \ **f** (MATLAB backslash operator).

### Maximum Entropy Weighted Least Squares (MEWLS)

MEWLS introduces a diagonal weight matrix **W** = diag(w₁, …, wₙ) and solves:

```
||W^(1/2) (A·QI_coeff - f)||²
```

Solution: **QI_coeff** = (A'WA)⁻¹ A'W**f**

Weights are computed as:

```
w_i = exp(-lambda · ||p_i - centroid||²)
```

where `centroid` is the mean of the point cloud and `lambda` controls the decay rate. Samples identified as outliers have their weight multiplied by `alpha_out` (default 0.5).

**Parameter guidance:**
- `lambda` = 1–5: mild weighting, close to OLS behaviour
- `lambda` = 10–20: moderate weighting (recommended starting range)
- `lambda` > 30: aggressive weighting, strong centroid bias
- `alpha_out` = 0.0: outliers are completely ignored
- `alpha_out` = 0.5: outliers contribute at half weight (default)
- `alpha_out` = 1.0: no penalty (equivalent to ignoring outlier flags)

### Error metrics

| Metric | Formula | Interpretation |
|--------|---------|----------------|
| RMSE | sqrt(mean((f_true - QI)²)) | Average error magnitude, penalises large errors |
| MaxAbsError | max(\|f_true - QI\|) | Worst-case point error |
| MAE | mean(\|f_true - QI\|) | Average absolute deviation |

## API reference

| Function | Signature | Description |
|----------|-----------|-------------|
| `adaptivehb.solvers.leastSquaresSolver` | `result = leastSquaresSolver(dataset, method_data)` | OLS polynomial surface fitter |
| `adaptivehb.solvers.mewlsSolver` | `result = mewlsSolver(dataset, method_data)` | MEWLS entropy-weighted fitter |
| `adaptivehb.solvers.polynomialDesignMatrix` | `[A, terms] = polynomialDesignMatrix(data, degree)` | 2-D polynomial design matrix |
| `adaptivehb.solvers.compute_metrics` | `metrics = compute_metrics(f_true, QI)` | Compute [RMSE, MaxAbs, MAE] |
| `adaptivehb.io.load_dataset` | `dataset = load_dataset(filePath, options)` | Load point-cloud data from text file |
| `adaptivehb.io.apply_noise` | `dataset = apply_noise(dataset, noiseCfg)` | Struct-level noise (JSON-driven) |
| `adaptivehb.data.apply_noise` | `[augmented, metadata] = apply_noise(x, y, f, settings)` | Vector-level noise (gauss/spike) |
| `adaptivehb.data.format_noise_filename` | `filename = format_noise_filename(baseName, metadata, options)` | Build noise-encoding filename |
| `adaptivehb.viz.plot_solution_surface` | `fig = plot_solution_surface(dataset, solverResults)` | 3-D surface comparison figure |
| `adaptivehb.viz.plot_error_distribution` | `fig = plot_error_distribution(dataset, solverResults)` | Error histogram figure |
| `adaptivehb.viz.plot_convergence_curve` | `fig = plot_convergence_curve(solverResults)` | RMSE convergence figure |

## Running the tests

```matlab
addpath('scripts');
run_tests();
```

The test suite validates all modules (loaders, noise injectors, solvers, visualisation, and the end-to-end pipeline) and prints a pass/fail summary. No toolboxes are required.

## Author

Francesco Pucci
PhD Candidate — University of Florence / INdAM

## License

This project is licensed under the Apache License 2.0. See [LICENSE](LICENSE) for details.
