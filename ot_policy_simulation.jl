##############################
# Simulation: OT-based policy effects (1D)
# Author: you
##############################

using Random
using Distributions
using Statistics
using StatsBase

# -----------------------------
# Utilities: empirical quantile with linear interpolation
# -----------------------------
"""
    emp_quantile(sorted_y::Vector{Float64}, p::Float64)

Compute empirical quantile at probability p in [0,1] using linear interpolation
over the sorted sample `sorted_y`.
"""
function emp_quantile(sorted_y::Vector{Float64}, p::Float64)
    n = length(sorted_y)
    @assert n >= 2 "Need at least 2 points for interpolation."
    p = clamp(p, 0.0, 1.0)
    # Use "Type 7" style index: h = 1 + (n-1)*p
    h = 1.0 + (n - 1) * p
    i = floor(Int, h)
    j = min(i + 1, n)
    γ = h - i
    # handle boundaries
    if i <= 1
        return sorted_y[1]
    elseif i >= n
        return sorted_y[n]
    else
        return (1 - γ) * sorted_y[i] + γ * sorted_y[j]
    end
end

"""
    emp_quantiles(sorted_y::Vector{Float64}, ps::Vector{Float64})

Vectorized version of `emp_quantile`.
"""
function emp_quantiles(sorted_y::Vector{Float64}, ps::Vector{Float64})
    return [emp_quantile(sorted_y, p) for p in ps]
end

# -----------------------------
# 1D empirical OT map via monotone rearrangement (rank matching)
# -----------------------------
"""
    ot_map_rankmatch(y0::Vector{Float64}, y1::Vector{Float64}; grid_n=200)

Constructs an empirical 1D OT map from distribution of y0 to y1 using
a common probability grid and empirical quantile interpolation.

Returns:
- ps: grid probabilities
- q0: quantiles of y0 at ps
- q1: quantiles of y1 at ps
- function T_hat(y): map y -> q1(F0_hat(y)) approximated via interpolation on (q0,q1)
"""
function ot_map_rankmatch(y0::Vector{Float64}, y1::Vector{Float64}; grid_n::Int=200)
    @assert length(y0) >= 5 "Too few controls for OT estimation."
    @assert length(y1) >= 5 "Too few treated for OT estimation."

    s0 = sort(y0)
    s1 = sort(y1)

    # Common grid away from 0 and 1 to avoid extreme interpolation noise
    ps = range(0.005, 0.995; length=grid_n) |> collect

    q0 = emp_quantiles(s0, ps)
    q1 = emp_quantiles(s1, ps)

    # Build an approximate inverse map y -> p using interpolation on q0:
    # If y lies between q0[k] and q0[k+1], interpolate p between ps[k] and ps[k+1].
    function F0_invprob(y::Float64)
        if y <= q0[1]
            return ps[1]
        elseif y >= q0[end]
            return ps[end]
        else
            k = searchsortedlast(q0, y)  # index k s.t. q0[k] <= y < q0[k+1]
            k = clamp(k, 1, length(q0) - 1)
            yL, yR = q0[k], q0[k+1]
            pL, pR = ps[k], ps[k+1]
            # linear interpolation in y
            t = (y - yL) / (yR - yL + eps())
            return (1 - t) * pL + t * pR
        end
    end

    # Transport map y -> q1(p) with p approx from control-quantile grid
    function T_hat(y::Float64)
        p = F0_invprob(y)
        # map to treated quantile at probability p (again via interpolation on (ps,q1))
        # Since q1 is already on ps grid, interpolate directly:
        if p <= ps[1]
            return q1[1]
        elseif p >= ps[end]
            return q1[end]
        else
            k = searchsortedlast(ps, p)
            k = clamp(k, 1, length(ps) - 1)
            pL, pR = ps[k], ps[k+1]
            yL, yR = q1[k], q1[k+1]
            t = (p - pL) / (pR - pL + eps())
            return (1 - t) * yL + t * yR
        end
    end

    return ps, q0, q1, T_hat
end

# -----------------------------
# Policy functionals Δ_phi
# -----------------------------
phi_mean(u, y=nothing) = u
phi_tail(u, y; c=0.0) = (y + u <= c) - (y <= c)   # needs y baseline
phi_loss(u, y; ℓ = y -> y^2) = ℓ(y + u) - ℓ(y)    # generic welfare/loss

# -----------------------------
# DGP: tail-only policy shift
# -----------------------------
"""
    simulate_data(n; kappa=0.2, tau=1.0, sigma=0.1, seed=1)

Creates i.i.d. sample (Y, D, X) where:
- X is binary (0/1) with P(X=1)=0.5
- D ~ Bernoulli(p(X)) independent of (Y0,Y1) given X  (unconfoundedness)
- Y0|X is skewed (lognormal + shift)
- Policy affects only lower tail: Y1 = Y0 - tau * 1{Y0 <= q_kappa(X)} + noise

Returns a NamedTuple with vectors: Y, D, X, Y0, Y1
"""
function simulate_data(n::Int; kappa=0.2, tau=1.0, sigma=0.1, seed=1)
    rng = MersenneTwister(seed)

    # covariate X ∈ {0,1}
    X = rand(rng, Bernoulli(0.5), n)

    # propensity depends on X (still unconfounded)
    p0, p1 = 0.35, 0.65
    p = [xi == 1 ? p1 : p0 for xi in X]
    D = [rand(rng) < p[i] ? 1 : 0 for i in 1:n]

    # baseline outcome Y0 | X is skewed
    # lognormal with different location by X
    μ = [xi == 1 ? 0.2 : 0.0 for xi in X]
    σ0 = 0.6
    Y0 = [rand(rng, LogNormal(μ[i], σ0)) - 1.0 for i in 1:n]  # shift so can be near 0

    # compute within-X tail cutoff q_kappa(X=x) on the population draw
    # (in a real application this is unknown; in simulation it defines the policy)
    qk = Dict{Int,Float64}()
    for xval in (0, 1)
        idx = findall(x -> x == xval, X)
        ys = Y0[idx]
        qk[xval] = quantile(ys, kappa)
    end

    # policy outcome: tail-only shift + additive noise
    ε = rand(rng, Normal(0, sigma), n)
    Y1 = similar(Y0)
    for i in 1:n
        cutoff = qk[X[i]]
        Y1[i] = Y0[i] - tau * (Y0[i] <= cutoff ? 1.0 : 0.0) + ε[i]
    end

    # observed
    Y = [D[i] == 1 ? Y1[i] : Y0[i] for i in 1:n]

    return (Y=Y, D=D, X=X, Y0=Y0, Y1=Y1, qk=qk)
end

# -----------------------------
# Estimation at fixed x: Δ_phi(x)
# -----------------------------
"""
    estimate_delta_phi_at_x(Y, D, X; xval=0, phi=:mean, c=0.0, grid_n=200)

Computes Δ̂_phi(x) using:
- μ̂0 from controls with X=x
- μ̂1 from treated with X=x
- 1D OT map via ot_map_rankmatch
- plug-in Δ̂_phi = mean over controls of φ(T̂(y)-y, y)

phi options: :mean, :tail, :loss
"""
function estimate_delta_phi_at_x(Y::Vector{Float64}, D::Vector{Int}, X::Vector{Int};
                                 xval::Int=0, phi::Symbol=:mean, c::Float64=0.0,
                                 grid_n::Int=200)
    idx0 = [i for i in eachindex(Y) if (X[i] == xval && D[i] == 0)]
    idx1 = [i for i in eachindex(Y) if (X[i] == xval && D[i] == 1)]

    y0 = Float64[Y[i] for i in idx0]
    y1 = Float64[Y[i] for i in idx1]

    # OT map
    _, _, _, T̂ = ot_map_rankmatch(y0, y1; grid_n=grid_n)

    # Compute shifts and plug-in functional
    shifts = [T̂(y) - y for y in y0]

    if phi == :mean
        return mean(shifts)
    elseif phi == :tail
        # tail probability shift at threshold c
        vals = [phi_tail(shifts[i], y0[i]; c=c) for i in eachindex(y0)]
        return mean(vals)
    elseif phi == :loss
        # quadratic loss by default
        vals = [phi_loss(shifts[i], y0[i]; ℓ = y -> y^2) for i in eachindex(y0)]
        return mean(vals)
    else
        error("Unknown phi option: $phi")
    end
end

# -----------------------------
# "True" Δ_phi(x) on the realized potential outcomes (Monte Carlo truth proxy)
# -----------------------------
"""
    true_delta_phi_at_x(Y0, Y1, X; xval=0, phi=:mean, c=0.0, grid_n=4000)

Computes a high-precision approximation to Δ_phi(x) using the realized potential outcomes
(Y0,Y1) in the simulated population draw. This is NOT available in practice; it's for coverage checks.

We treat μ0|x as empirical distribution of Y0 among X=x, μ1|x as empirical distribution of Y1 among X=x,
then compute OT map via rankmatch and plug-in functional.
"""
function true_delta_phi_at_x(Y0::Vector{Float64}, Y1::Vector{Float64}, X::Vector{Int};
                             xval::Int=0, phi::Symbol=:mean, c::Float64=0.0,
                             grid_n::Int=4000)
    idx = findall(x -> x == xval, X)
    y0 = sort(Float64[Y0[i] for i in idx])
    y1 = sort(Float64[Y1[i] for i in idx])

    # OT via quantile grid
    ps = range(0.001, 0.999; length=grid_n) |> collect
    q0 = emp_quantiles(y0, ps)
    q1 = emp_quantiles(y1, ps)
    shifts = [q1[k] - q0[k] for k in eachindex(ps)]

    if phi == :mean
        return mean(shifts)
    elseif phi == :tail
        # Approx with baseline y at q0 grid
        vals = [(q0[k] + shifts[k] <= c) - (q0[k] <= c) for k in eachindex(ps)]
        return mean(vals)
    elseif phi == :loss
        vals = [(q0[k] + shifts[k])^2 - (q0[k])^2 for k in eachindex(ps)]
        return mean(vals)
    else
        error("Unknown phi option: $phi")
    end
end

# -----------------------------
# Bootstrap CI for Δ̂_phi(x)
# -----------------------------
"""
    bootstrap_ci_delta(Y, D, X; xval=0, phi=:mean, c=0.0, B=200, α=0.05, seed=1)

Returns:
- point estimate
- bootstrap SE
- normal CI
- percentile CI
"""
function bootstrap_ci_delta(Y::Vector{Float64}, D::Vector{Int}, X::Vector{Int};
                            xval::Int=0, phi::Symbol=:mean, c::Float64=0.0,
                            B::Int=200, α::Float64=0.05, seed::Int=1)
    rng = MersenneTwister(seed)
    n = length(Y)

    θ̂ = estimate_delta_phi_at_x(Y, D, X; xval=xval, phi=phi, c=c)

    boots = Vector{Float64}(undef, B)
    for b in 1:B
        idx = rand(rng, 1:n, n)  # resample indices
        Yb = Y[idx]
        Db = D[idx]
        Xb = X[idx]
        boots[b] = estimate_delta_phi_at_x(Yb, Db, Xb; xval=xval, phi=phi, c=c)
    end

    se = std(boots)
    z = quantile(Normal(), 1 - α/2)

    ci_normal = (θ̂ - z*se, θ̂ + z*se)
    ci_perc   = (quantile(boots, α/2), quantile(boots, 1 - α/2))

    return (theta_hat=θ̂, se=se, ci_normal=ci_normal, ci_percentile=ci_perc, boots=boots)
end

# -----------------------------
# Monte Carlo study
# -----------------------------
"""
    run_mc(R; n=1000, xval=0, phi=:mean, c=0.0, B=200, α=0.05, base_seed=123)

Runs R Monte Carlo replications.
Reports bias, RMSE, and coverage of bootstrap normal CI + percentile CI.
Truth is approximated using potential outcomes from each replication.
"""
function run_mc(R::Int; n::Int=1000, xval::Int=0, phi::Symbol=:mean, c::Float64=0.0,
                B::Int=200, α::Float64=0.05, base_seed::Int=123)

    θhats = Float64[]
    θtrues = Float64[]
    cover_normal = 0
    cover_perc = 0

    for r in 1:R
        dat = simulate_data(n; seed=base_seed + r)
        Y = dat.Y
        D = Int.(dat.D)
        X = Int.(dat.X)

        θtrue = true_delta_phi_at_x(dat.Y0, dat.Y1, X; xval=xval, phi=phi, c=c)
        out = bootstrap_ci_delta(Y, D, X; xval=xval, phi=phi, c=c, B=B, α=α, seed=10_000 + r)

        push!(θhats, out.theta_hat)
        push!(θtrues, θtrue)

        if out.ci_normal[1] <= θtrue <= out.ci_normal[2]
            cover_normal += 1
        end
        if out.ci_percentile[1] <= θtrue <= out.ci_percentile[2]
            cover_perc += 1
        end
    end

    errs = θhats .- θtrues
    bias = mean(errs)
    rmse = sqrt(mean(errs .^ 2))
    covN = cover_normal / R
    covP = cover_perc / R

    return (bias=bias, rmse=rmse, coverage_normal=covN, coverage_percentile=covP,
            thetahat_mean=mean(θhats), thetatrue_mean=mean(θtrues))
end

# -----------------------------
# Example run
# -----------------------------
if abspath(PROGRAM_FILE) == @__FILE__
    # Example: evaluate mean shift for X=0
    R = 50
    res = run_mc(R; n=1200, xval=0, phi=:mean, c=0.0, B=150, α=0.05, base_seed=2026)
    println("Monte Carlo results (R=$R):")
    println(res)

    # Example: tail probability shift at threshold c=0.0
    res_tail = run_mc(R; n=1200, xval=0, phi=:tail, c=0.0, B=150, α=0.05, base_seed=3030)
    println("\nTail-shift results (R=$R, c=0.0):")
    println(res_tail)
end
