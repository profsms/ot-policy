# replicate_all.jl
# Simple driver to produce tables and the shift figure (example).
#
# Usage:
#   julia replicate_all.jl            # runs both tables and figure with defaults
#   julia replicate_all.jl --skip-tables
#   julia replicate_all.jl --skip-figure
#   julia replicate_all.jl --only-tables
#   julia replicate_all.jl --only-figure
#
# You can also set a few environment variables to override defaults:
#   REPLICATE_NS     -> comma-separated ns for tables (e.g. "600,1200,2400")
#   REPLICATE_R      -> R (Int)
#   REPLICATE_B      -> B (Int)
#   REPLICATE_XVAL   -> xval (numeric)
#   REPLICATE_C      -> c (numeric)
#   SHIFT_N          -> n for shift figure (Int)
#   SHIFT_SEED       -> seed for shift figure (Int)
#   SHIFT_OUTPATH    -> output path for shift figure (string)
#
# The functions exported by the package's scripts are expected to be:
#   export_tables(; outdir=..., ns=..., R=..., B=..., xval=..., c=...)
#   make_shift_figure_one_rep(; n=..., xval=..., seed=..., outpath=...)
#
# This script loads the repository files:
include("csv_mc_grid.jl")
include("make_shift_figure_one_rep.jl")

# -----------------------
# Default parameters
# -----------------------
default_ns = [600, 1200, 2400]
default_R = 200
default_B = 300
default_xval = 0
default_c = 0.0
default_tables_outdir = "tables"

default_shift_n = 2000
default_shift_xval = 0
default_shift_seed = 2026
default_shift_outpath = "figs/shift_curve_x0.pdf"

# -----------------------
# Helpers: parse env / ARGS
# -----------------------
function parse_ns(s::AbstractString)
    parts = split(strip(s), ',')
    try
        return [parse(Int, strip(p)) for p in parts if !isempty(strip(p))]
    catch
        @warn "Could not parse REPLICATE_NS, falling back to defaults" s
        return default_ns
    end
end

# Read environment overrides (if provided)
ns_env = haskey(ENV, "REPLICATE_NS") ? parse_ns(ENV["REPLICATE_NS"]) : default_ns
R_env = haskey(ENV, "REPLICATE_R") ? parse(Int, ENV["REPLICATE_R"]) : default_R
B_env = haskey(ENV, "REPLICATE_B") ? parse(Int, ENV["REPLICATE_B"]) : default_B
xval_env = haskey(ENV, "REPLICATE_XVAL") ? parse(Float64, ENV["REPLICATE_XVAL"]) : default_xval
c_env = haskey(ENV, "REPLICATE_C") ? parse(Float64, ENV["REPLICATE_C"]) : default_c
tables_outdir_env = haskey(ENV, "REPLICATE_OUTDIR") ? ENV["REPLICATE_OUTDIR"] : default_tables_outdir

shift_n_env = haskey(ENV, "SHIFT_N") ? parse(Int, ENV["SHIFT_N"]) : default_shift_n
shift_xval_env = haskey(ENV, "SHIFT_XVAL") ? parse(Float64, ENV["SHIFT_XVAL"]) : default_shift_xval
shift_seed_env = haskey(ENV, "SHIFT_SEED") ? parse(Int, ENV["SHIFT_SEED"]) : default_shift_seed
shift_outpath_env = haskey(ENV, "SHIFT_OUTPATH") ? ENV["SHIFT_OUTPATH"] : default_shift_outpath

# CLI flags
skip_tables = "--skip-tables" in ARGS
skip_figure = "--skip-figure" in ARGS
only_tables = "--only-tables" in ARGS
only_figure = "--only-figure" in ARGS

if only_tables
    skip_figure = true
end
if only_figure
    skip_tables = true
end

# -----------------------
# Run
# -----------------------
function ensure_parent_dir(path::AbstractString)
    dir = dirname(path)
    mkpath(dir)
end

function main()
    println("replicate_all.jl — replication driver")
    println("Options summary:")
    println("  Tables: ", !skip_tables)
    println("  Shift figure: ", !skip_figure)
    println()

    # Run tables
    if !skip_tables
        println("Preparing to export tables...")
        mkpath(tables_outdir_env)
        println("  outdir = ", tables_outdir_env)
        println("  ns = ", ns_env)
        println("  R = ", R_env, ", B = ", B_env)
        println("  xval = ", xval_env, ", c = ", c_env)
        try
            # export_tables is expected to be defined in csv_mc_grid.jl
            export_tables(; outdir = tables_outdir_env,
                             ns = ns_env,
                             R = R_env,
                             B = B_env,
                             xval = xval_env,
                             c = c_env)
            println("Tables exported to: ", tables_outdir_env)
        catch err
            @error "export_tables failed" exception=(err, catch_backtrace())
            rethrow(err)
        end
    else
        println("Skipping tables step.")
    end

    # Make representative shift figure
    if !skip_figure
        println("Preparing to make shift figure (one representative replication)...")
        ensure_parent_dir(shift_outpath_env)
        println("  n = ", shift_n_env, ", xval = ", shift_xval_env, ", seed = ", shift_seed_env)
        println("  outpath = ", shift_outpath_env)
        try
            # make_shift_figure_one_rep is expected to be defined in make_shift_figure_one_rep.jl
            make_shift_figure_one_rep(n = shift_n_env,
                                     xval = shift_xval_env,
                                     seed = shift_seed_env,
                                     outpath = shift_outpath_env)
            println("Shift figure written to: ", shift_outpath_env)
        catch err
            @error "make_shift_figure_one_rep failed" exception=(err, catch_backtrace())
            rethrow(err)
        end
    else
        println("Skipping shift figure step.")
    end

    println("Done.")
end

# If this script is run directly, execute main()
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end