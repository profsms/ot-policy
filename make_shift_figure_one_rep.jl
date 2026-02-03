function make_shift_figure_one_rep(; n=2000, xval=0, seed=2026, outpath="figs/shift_curve_x0.pdf")
    dat = simulate_data(n; seed=seed)  # from your original script
    Y = dat.Y
    D = Int.(dat.D)
    X = Int.(dat.X)

    # observed within stratum x: split by D
    idx0 = [i for i in eachindex(Y) if X[i]==xval && D[i]==0]
    idx1 = [i for i in eachindex(Y) if X[i]==xval && D[i]==1]
    y0_obs = Float64[Y[i] for i in idx0]
    y1_obs = Float64[Y[i] for i in idx1]

    ps, _, _, δ_hat = quantile_shift_curve(y0_obs, y1_obs; grid_n=200)

    # optional: "truth" within the same replication using potential outcomes
    idx = findall(x -> x == xval, X)
    y0_true = Float64[dat.Y0[i] for i in idx]
    y1_true = Float64[dat.Y1[i] for i in idx]
    _, _, _, δ_true = quantile_shift_curve(y0_true, y1_true; grid_n=200)

    # plot
    mkpath(dirname(outpath))
    fig = Figure(size=(700, 420))
    ax = Axis(fig[1,1], xlabel="p (quantile index)", ylabel="shift  q1(p) - q0(p)")
    lines!(ax, ps, δ_hat, label="estimated")
    lines!(ax, ps, δ_true, linestyle=:dash, label="truth (from potential outcomes)")
    axislegend(ax, position=:rb)
    save(outpath, fig)
    return outpath
end

# run it:
# make_shift_figure_one_rep()
