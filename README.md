# ot-policy — Replication package for
Distributional Policy Evaluation via Optimal Transport (preprint)

This repository contains the runnable replication code and materials for the preprint
"Distributional Policy Evaluation via Optimal Transport" (submitted to The Econometric Journal).
The codebase is written in Julia and reproduces the paper's Monte Carlo tables and the representative
quantile-shift figure.

Julia version
- Use Julia 1.11.2 (recommended for reproducibility).

Required Julia packages
Install the packages listed below before running the scripts:

- CSV
- DataFrames
- Distributions
- StatsBase
- CairoMakie

Install example (run inside a Julia REPL):

```julia
using Pkg
Pkg.add.("CSV","DataFrames","Distributions","StatsBase","CairoMakie")
Pkg.precompile()
```

Repository files
- `Distributional_Policy_Evaluation_via_Optimal_Transport.pdf` — the article (PDF).
- `replicate_all.jl` — replication driver (runs tables + representative figure; supports CLI flags and environment-variable overrides).
- `csv_mc_grid.jl` — Monte Carlo grid runner; provides `export_tables(...)`.
- `make_shift_figure_one_rep.jl` — figure builder; provides `make_shift_figure_one_rep(...)`.
- `ot_policy_simulation.jl` — core DGP, OT estimation, bootstrap, Monte Carlo routines.
- `sim_outputs.jl` — quantile/Wasserstein utilities and plotting helpers.
- `run_mc_delta.jl` — auxiliary/placeholder driver (may be empty).
- `tables/` — output directory for table CSVs (created by scripts).
- `figs/` — output directory for figure PDFs (created by scripts).

How to run

1. Clone the repository and change to its root:

```bash
git clone https://github.com/profsms/ot-policy.git
cd ot-policy
```

2. Install the required packages (see "Required Julia packages" above).

3. Run the full replication

A — From the shell (repository root):

```bash
julia replicate_all.jl
```

B — From the Julia REPL (interactive):

```julia
include("replicate_all.jl")
main()
```

CLI flags
- `--skip-tables` — skip table export.
- `--skip-figure` — skip figure creation.
- `--only-tables` — run only tables.
- `--only-figure` — run only the figure.

Environment-variable overrides
You may override defaults by setting environment variables when launching Julia.

- `REPLICATE_NS` — comma-separated sample sizes for the grid (e.g. "600,1200,2400").
- `REPLICATE_R` — Monte Carlo repetitions (integer R).
- `REPLICATE_B` — bootstrap draws (integer B).
- `REPLICATE_XVAL` — covariate value used by the experiments (numeric).
- `REPLICATE_C` — threshold parameter for tail functional (numeric).
- `REPLICATE_OUTDIR` — output directory for tables (default `tables`).
- `SHIFT_N` — sample size for the representative shift figure (integer).
- `SHIFT_XVAL` — covariate value used by the figure (numeric).
- `SHIFT_SEED` — RNG seed for reproducibility of the figure (integer).
- `SHIFT_OUTPATH` — output path for the figure (default `figs/shift_curve_x0.pdf`).

Example with environment overrides:

```bash
REPLICATE_NS="600,1200,2400" REPLICATE_R=200 REPLICATE_B=300 REPLICATE_XVAL=0 REPLICATE_C=0.0 \
  julia replicate_all.jl

SHIFT_N=2000 SHIFT_SEED=2026 SHIFT_OUTPATH="figs/shift_curve_x0.pdf" \
  julia replicate_all.jl --skip-tables
```

Quick smoke tests
To check the pipeline quickly, reduce Monte Carlo and bootstrap workload:

```bash
REPLICATE_R=20 REPLICATE_B=20 julia replicate_all.jl --only-tables
SHIFT_N=500 SHIFT_SEED=42 julia replicate_all.jl --only-figure
```

Expected outputs
- `tables/` — CSV files produced by `export_tables`, typically `sim_mean.csv` and `sim_tail.csv`. Each row corresponds to a sample-size cell and contains metrics such as bias, RMSE, and coverage.
- `figs/shift_curve_x0.pdf` — representative quantile-shift figure produced by `make_shift_figure_one_rep` (path configurable via `SHIFT_OUTPATH`).

Notes and troubleshooting
- If a script fails with a missing-package error, install the named package using `Pkg.add("PackageName")` in the Julia REPL.
- If figure rendering fails, ensure `CairoMakie` is installed and that a suitable graphical environment is available for PDF rendering.
- The Monte Carlo runtime is controlled by `REPLICATE_R` and `REPLICATE_B`; raising these yields more accurate results but increases runtime.

Reproducibility
- Use Julia 1.11.2 and the package list above for reproducible results.
- The driver and scripts set RNG seeds in example runs; to reproduce identical outputs, keep the default seeds and parameters.

Citation
If you use this code or reproduce results, cite the accompanying preprint:
"Distributional Policy Evaluation via Optimal Transport" — preprint submitted to The Econometric Journal (see included PDF).

Contact / issues
Open an issue in this repository with:
- Julia version (`julia -v`),
- the exact command you ran,
- the error message or backtrace (if any).

License
- See repository for any license file included. If no license is present, request permission before reusing substantial portions of the code or results.