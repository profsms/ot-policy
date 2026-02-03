using Random, Distributions, Statistics, StatsBase
using CSV, DataFrames
using CairoMakie

# -----------------------------
# Quantile grid utilities (truth & estimated)
# -----------------------------
function quantile_shift_curve(y0::Vector{Float64}, y1::Vector{Float64}; grid_n=200)
    s0 = sort(y0); s1 = sort(y1)
    ps = collect(range(0.005, 0.995; length=grid_n))
    q0 = emp_quantiles(s0, ps)
    q1 = emp_quantiles(s1, ps)
    δ  = q1 .- q0
    return ps, q0, q1, δ
end

# Wasserstein-2 distance in 1D via quantiles: W2^2 = ∫ (q1-q0)^2 dp
function w2_1d(y0::Vector{Float64}, y1::Vector{Float64}; grid_n=500)
    ps, _, _, δ = quantile_shift_curve(y0, y1; grid_n=grid_n)
    # approximate integral by average over grid
    return sqrt(mean(δ .^ 2))
end

# QTE at specific probabilities
function qte_at_ps(y0::Vector{Float64}, y1::Vector{Float64}, ps::Vector{Float64})
    s0 = sort(y0); s1 = sort(y1)
    q0 = emp_quantiles(s0, ps)
    q1 = emp_quantiles(s1, ps)
    return q1 .- q0
end
