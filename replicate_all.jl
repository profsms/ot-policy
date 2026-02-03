# replicate_all.jl
# Simple driver to produce tables and the shift figure (example).
# This script creates output directories, includes the necessary helper scripts,
# and runs the grid/table generation and an example shift figure.

# Create output directories if they don't exist
mkpath("tables")
mkpath("figs")

# Include supporting scripts (must exist in the same repo)
include("csv_mc_grid.jl")
include("make_shift_figure_one_rep.jl")

# Default parameters (can be edited before running)
s = [600, 1200, 2400]   # sample sizes for the table grid
R = 200                 # number of Monte Carlo repetitions per cell (controls runtime)
B = 300                 # bootstrap replications (controls runtime)
xval = 0               # covariate value used in the examples
c = 0.0                 # cost/regularization parameter used in examples

# Run grid to produce tables (warning: R and B determine runtime)
println("Starting export_tables with ns=", ns, " R=", R, " B=", B)
export_tables(outdir="tables", ns=ns, R=R, B=B, xval=xval, c=c)
println("Tables written to tables/")

# Make one representative shift figure
fig_n = 2000
fig_seed = 2026
fig_out = "figs/shift_curve_x0.pdf"
println("Creating representative shift figure: ", fig_out)
make_shift_figure_one_rep(n=fig_n, xval=xval, seed=fig_seed, outpath=fig_out)
println("Figure written to ", fig_out)