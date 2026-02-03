"""
Runs a grid over sample sizes for a given functional and returns a DataFrame:
n, bias, rmse, coverage_normal, coverage_percentile.
"""
function run_grid(; ns=[600,1200,2400], R=200, B=300, xval=0, phi=:mean, c=0.0, base_seed=2026)
    rows = DataFrame(n=Int[], bias=Float64[], rmse=Float64[],
                     cov_normal=Float64[], cov_percentile=Float64[])
    for n in ns
        res = run_mc(R; n=n, xval=xval, phi=phi, c=c, B=B, α=0.05, base_seed=base_seed + 10n)
        push!(rows, (n, res.bias, res.rmse, res.coverage_normal, res.coverage_percentile))
        println("Done: phi=$phi n=$n => ", res)
    end
    return rows
end

function export_tables(; outdir="tables", ns=[600,1200,2400], R=200, B=300, xval=0, c=0.0)
    mkpath(outdir)

    df_mean = run_grid(; ns=ns, R=R, B=B, xval=xval, phi=:mean, c=c, base_seed=2026)
    CSV.write(joinpath(outdir, "sim_mean.csv"), df_mean)

    df_tail = run_grid(; ns=ns, R=R, B=B, xval=xval, phi=:tail, c=c, base_seed=3030)
    CSV.write(joinpath(outdir, "sim_tail.csv"), df_tail)

    # Optional: Wasserstein magnitude (truth & estimate) — see function below
    return (mean=df_mean, tail=df_tail)
end
