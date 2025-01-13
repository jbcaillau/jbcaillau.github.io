# Solve Goddard problem using ADNLPModels.jl

using ADNLPModels
using NLPModelsIpopt

# Parameters

const n = 3 # State dim
const m = 1 # Control dim
const Cd = 310 # Drag (1/2)
const β = 500 # Drag (2/2)
const Tmax = 3.5 # Max thrust
const b = 2 # Fuel consumption

const r0 = 1 # Initial altitude
const v0 = 0 # Initial speed
const m0 = 1 # Initial const mass
const vmax = 0.1 # Maximal authorized speed
const mf = 0.6 # Final mass to target

N = 100 
tol = 1e-7 # 1e-5

# ADNLPModel
# z = [tf, x[:], u[:]], tf : 1, x : n x (N + 1), u : m x (N + 1)

@inline function get(z, i, j, d; offset=0)
    k = offset + i + (j - 1) * d # j is assumed to start from 1 
    return z[k]
end

@inline function get(z, j, d; offset=0)
    k = offset + 1 + (j - 1) * d
    return @view z[k:k + d - 1]
end

c_dim = 5 + n + 6N # constraints size

@inline _tf(z) = z[1]
@inline _x(z, i, j) = get(z, i, j, n; offset=1)
@inline _x(z, j) = get(z, j, n; offset=1)
@inline _u(z, i, j) = get(z, i, j, m; offset=1 + (N + 1) * n)
@inline _u(z, j) = get(z, j, m; offset=1 + (N + 1) * n)

## objective: -r[N+1]

f(z) = -_x(z, 1, N+1)

## constraints

dr(r, v, m, u) = v
dv(r, v, m, u) = -Cd * v^2 * exp(-β * (r - 1)) / m - 1 / r^2 + u * Tmax / m
dm(r, v, m, u) = -b * Tmax * u
rk2(x1, x2, rhs1, rhs2, dt) = x2 - x1 - dt / 2 * (rhs1 + rhs2)

function con!(c, z) # todo: add N as a parameter of the model
    k = 1
    dt = _tf(z) / N

    # 0 ≤ tf
    c[k] = _tf(z)
    k = k + 1

    # x[:, 1] - [r0, v0, m0] == 0
    c[k:k + n - 1] .= _x(z, 1) - [r0, v0, m0]
    k = k + n

    # x[:, N + 1] - mf == 0
    c[k] = _x(z, 3, N + 1) - mf
    k = k + 1 

    # 0 ≤ u[1, :] ≤ 1 
    for j ∈ 1:N + 1
        c[k + j - 1] = _u(z, 1, j)
    end
    k = k + N + 1

    # r0 ≤ x[1, :]
    for j ∈ 1:N + 1
        c[k + j - 1] = _x(z, 1, j)
    end
    k = k + N + 1

    # 0 ≤ x[2, :] ≤ vmax
    for j ∈ 1:N + 1
        c[k + j - 1] = _x(z, 2, j)
    end
    k = k + N + 1

    # rk2 on r
    dj = dr(_x(z, 1, 1), _x(z, 2, 1), _x(z, 3, 1), _u(z, 1, 1)) 
    for j ∈ 1:N
        dj1 = dr(_x(z, 1, j + 1), _x(z, 2, j + 1), _x(z, 3, j + 1), _u(z, 1, j + 1)) 
        c[k + j - 1] = rk2(_x(z, 1, j), _x(z, 1, j + 1), dj, dj1, dt)
        dj = dj1
    end
    k = k + N

    # rk2 on v
    dj = dv(_x(z, 1, 1), _x(z, 2, 1), _x(z, 3, 1), _u(z, 1, 1)) 
    for j ∈ 1:n
        dj1 = dv(_x(z, 1, j + 1), _x(z, 2, j + 1), _x(z, 3, j + 1), _u(z, 1, j + 1)) 
        c[k + j - 1] = rk2(_x(z, 2, j), _x(z, 2, j + 1), dj, dj1, dt)
        dj = dj1
    end
    k = k + N

    # rk2 on m
    dj = dm(_x(z, 1, 1), _x(z, 2, 1), _x(z, 3, 1), _u(z, 1, 1)) 
    for j ∈ 1:n
        dj1 = dm(_x(z, 1, j + 1), _x(z, 2, j + 1), _x(z, 3, j + 1), _u(z, 1, j + 1)) 
        c[k + j - 1] = rk2(_x(z, 3, j), _x(z, 3, j + 1), dj, dj1, dt)
        dj = dj1
    end
    k = k + N

    @assert(k - 1 == c_dim) # check constraints size

    return nothing

end

## constraints bounds

lcon = zeros(c_dim)
ucon = zeros(c_dim)

k = 1

# 0 ≤ tf
lcon[k] = 0
ucon[k] = Inf
k = k + 1

# x[:, 1] - [r0, v0, m0] == 0
lcon[k:k + n - 1] .= 0
ucon[k:k + n - 1] .= 0
k = k + n

# x[:, N + 1] - mf == 0
lcon[k] = 0
ucon[k] = 0
k = k + 1 

# 0 ≤ u[1, :] ≤ 1 
for j ∈ 1:N + 1
    lcon[k] = 0
    ucon[k] = 1
end
k = k + N + 1
# r0 ≤ x[1, :]
for j ∈ 1:N + 1
    lcon[k] = r0
    ucon[k] = Inf
end
k = k + N + 1

# 0 ≤ x[2, :] ≤ vmax
for j ∈ 1:N + 1
    lcon[k] = 0
    ucon[k] = vmax
end
k = k + N + 1

# rk2 on r
for j ∈ 1:N
    lcon[k + j - 1] = 0
    ucon[k + j - 1] = 0
end
k = k + N

# rk2 on v
for j ∈ 1:n
    lcon[k + j - 1] = 0
    ucon[k + j - 1] = 0
end
k = k + N

# rk2 on m
for j ∈ 1:n
    lcon[k + j - 1] = 0
    ucon[k + j - 1] = 0
end
k = k + N

@assert(k - 1 == c_dim) # check constraints size