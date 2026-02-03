# ot-policy — Replication package for
"Distributional Policy Evaluation via Optimal Transport" (preprint)

This repository contains the replication code and materials for the paper "Distributional Policy Evaluation via Optimal Transport" (preprint, submitted to The Econometric Journal). The codebase is written in Julia and provides scripts to reproduce the paper's Monte Carlo tables and the representative shift figure.

This README documents the exact files in the bundle and shows how to run the replication using `replicate_all.jl`. It does NOT assume a Project/Manifest file — install packages manually as described below.

---

## Repository files (root)
- `Distributional_Policy_Evaluation_via_Optimal_Transport.pdf` — the paper (PDF).
- `replicate_all.jl` — replication driver that runs the Monte Carlo table export and the representative shift figure. (Supports CLI flags and environment overrides.)
- `csv_mc_grid.jl` — grid routine that calls the Monte Carlo runner and exports tables (defines `export_tables(...)`).
- `make_shift_figure_one_rep.jl` — builds and saves the representative shift figure (defines `make_shift_figure_one_rep(...)`).
- `ot_policy_simulation.jl` — core DGP, OT estimation, bootstrap, and Monte Carlo routines (defines `simulate_data`, `ot_map_rankmatch`, `estimate_delta_phi_at_x`, `bootstrap_ci_delta`, `run_mc`, etc.).
- `sim_outputs.jl` — quantile utilities, W2 and plotting helpers used by the figure script.
- `run_mc_delta.jl` — (placeholder / optional) additional driver (may be empty).
- `tables/` — output directory where table CSVs are written by the grid (created by the scripts).
- `figs/` — output directory where figure PDFs are written by the figure script (created by the scripts).

(If you list other files locally, they belong to auxiliary workflows; the files above are the main scripts used by the replication driver.)

---

## Requirements
- Julia 1.11.2 (please use this exact version for reproducibility).
- The codebase does not include a Project.toml/Manifest.toml; install packages manually in your Julia environment.

Suggested packages to install before running (install only those the scripts require):
```julia
using Pkg
Pkg.add.("CSV", "DataFrames", "Distributions", "StatsBase", "Random", "CairoMakie")
```
- If a script fails with a missing-package error, add the named package with `Pkg.add("Name")`.
- Standard library modules (e.g. `Statistics`, `LinearAlgebra`) do not require `Pkg.add`.

---

## How to run the replication (replicate_all.jl)
The driver `replicate_all.jl` orchestrates table generation and figure creation. It supports CLI flags and environment-variable overrides. Example ways to run:

1) Run both tables and figure with defaults (from repo root):
```bash
julia replicate_all.jl
# or inside Julia REPL:
# include("replicate_all.jl"); main()
```

2) Run only tables or only figure using CLI flags:
```bash
julia replicate_all.jl --only-tables
julia replicate_all.jl --only-figure
julia replicate_all.jl --skip-figure
julia replicate_all.jl --skip-tables
```

3) Override defaults with environment variables (examples):
```bash
REPLICATE_NS="600,1200,2400" REPLICATE_R=200 REPLICATE_B=300 REPLICATE_XVAL=0 REPLICATE_C=0.0 \
julia replicate_all.jl

SHIFT_N=2000 SHIFT_SEED=2026 SHIFT_OUTPATH="figs/shift_curve_x0.pdf" julia replicate_all.jl --skip-tables
```

### Environment variables recognized by `replicate_all.jl`
- `REPLICATE_NS` — comma-separated ns for table grid (e.g. `"600,1200,2400"`).
- `REPLICATE_R` — Monte Carlo repetitions (Int).
- `REPLICATE_B` — bootstrap draws (Int).
- `REPLICATE_XVAL` — numeric covariate value.
- `REPLICATE_C` — numeric threshold parameter.
- `REPLICATE_OUTDIR` — tables output directory (default `tables`).
- `SHIFT_N` — sample size for the shift figure (Int).
- `SHIFT_XVAL` — xval for the shift figure (numeric).
- `SHIFT_SEED` — RNG seed for the shift figure (Int).
- `SHIFT_OUTPATH` — output path for the figure (default `figs/shift_curve_x0.pdf`).

### CLI flags
- `--skip-tables` — skip table export.
- `--skip-figure` — skip figure creation.
- `--only-tables` — run tables only.
- `--only-figure` — run figure only.

---

## Typical outputs
- `tables/sim_mean.csv` and `tables/sim_tail.csv` (or similar filenames) — CSVs written by `export_tables`. Each row corresponds to an `n` in the grid and columns contain `bias`, `rmse`, `coverage_normal`, `coverage_percentile`, etc.
- `figs/shift_curve_x0.pdf` — the representative quantile-shift figure produced by the figure script.

You can run the driver interactively (as you did):
```julia
include("replicate_all.jl")
main()

# Example console output (your run):
# replicate_all.jl — replication driver
# Options summary:
#   Tables: true
#   Shift figure: true
# ...
# Done.
```

This confirms the driver calls `export_tables(...)` and `make_shift_figure_one_rep(...)` and writes outputs into `tables/` and `figs/`.

---

## Quick testing tips
- Reduce run time by lowering `REPLICATE_R` and `REPLICATE_B` (e.g. `R=20,B=20`) for smoke tests.
- If plotting fails, ensure `CairoMakie` (or plotting backend used) is installed.
- If a script raises a `MethodError` or missing symbol, confirm you included the correct script (e.g. `include("ot_policy_simulation.jl")`) or run the `replicate_all.jl` driver which includes the necessary scripts.

---

## If you want changes
I can: add more detailed README sections, make `replicate_all.jl` accept command-line args, or pin package versions by providing a Project/Manifest — tell me which you prefer.
