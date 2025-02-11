# fooo.jl

using OptimalControl
using NLPModelsIpopt
using Plots

const τ = 1.0
const σ = 49.0
const fmax = 10.0
const T = 20.0
const e0 = 100.0

ocp = @def begin
    t ∈ [0, T], time
    x = (v, e) ∈ R², state
    f ∈ R, control
    ẋ(t) == [-v(t) / τ + f(t), σ - f(t) * v(t)]
    v(0) == 0
    e(0) == e0
    0 ≤ f(t) ≤ fmax
    e(t) ≥ 0
    ∫( -v(t) ) → min
end

N = 1000
tol = 1e-4

sol = solve(ocp; grid_size=N, tol=tol)
plot(sol)

ocp_r(ε) = @def begin
    t ∈ [0, T], time
    x = (v, e) ∈ R², state
    f ∈ R, control
    ẋ(t) == [-v(t) / τ + f(t), σ - f(t) * v(t)]
    v(0) == 0
    e(0) == e0
    0 ≤ f(t) ≤ fmax
    e(t) ≥ 0
    ∫( -v(t) + ε * f(t)^2 ) → min
end

ε = 1e-4 # L^2 regularisation

sol_r = solve(ocp_r(ε); grid_size=N, tol=tol)
plot(sol_r)