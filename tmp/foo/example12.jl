# Example 12 from
# Ansel, Q. et al. Introduction to theoretical and experimental aspects of quantum optimal control,
# J. of Phys. B: Atomic, Molecular and Optical Physics 57.13 (2024), 133001. doi: 10.1088/1361-6455/ad46a5

using OptimalControl
using NLPModelsIpopt 
using OrdinaryDiffEq  
using MINPACK      
using Plots           
using LinearAlgebra
using DifferentiationInterface
import ForwardDiff

# -----------------  Direct solve  --------------------

q0 = [0.0, 0.0, 1.0] 
qf = [0.0, 0.0, -1.0] 
Δ = 0.5
tf = 2 * π / (√( 1 + Δ^2))
p_zero = 0.1
down = [0.0, 0.0, -1.0] 

### ocp, unconstrained

ocp1 = @def begin
    t ∈ [0, tf], time
    q = [x, y, z] ∈ R³, state
    u ∈ R, control
    
    q(0) == q0 
    q̇(t) == [ -Δ * y(t),
             Δ * x(t) - u(t) * z(t),
             u(t) * y(t) ]
    sum((q(tf) - qf).^2) + (p_zero / 2) * ∫(u(t)^2) → min
end

N = 100
direct_sol1 = solve(ocp1; grid_size = N)
plt = plot(direct_sol1; size=(1000,1000),background_color=:transparent)

### ocp with constraits (1/2)

ocp2 = @def begin
    t ∈ [0, tf], time
    q = [x, y, z] ∈ R³, state
    u ∈ R, control
    
    q(0) == q0 
    x(t)^2 + y(t)^2 + z(t)^2 ≤ 1
    q̇(t) == [ -Δ * y(t),
             Δ * x(t) - u(t) * z(t),
             u(t) * y(t) ]
    sum((q(tf) - qf).^2) + (p_zero / 2) * ∫(u(t)^2) → min
end

direct_sol2 = solve(ocp2; grid_size = N)
plot!(plt,direct_sol2)

### ocp with constraits (2/2)

ocp3 = @def begin
    t ∈ [0, tf], time
    q = [x, y, z] ∈ R³, state
    u ∈ R, control
    
    q(0) == q0 
    x(t)^2 + y(t)^2 + z(t)^2 ≤ 1
    -1 ≤ x(t) ≤ 1
    -1 ≤ y(t) ≤ 1
    -1 ≤ z(t) ≤ 1
    q̇(t) == [ -Δ * y(t),
             Δ * x(t) - u(t) * z(t),
             u(t) * y(t) ]
    sum((q(tf) - qf).^2) + (p_zero / 2) * ∫(u(t)^2) → min
end

direct_sol3 = solve(ocp3; grid_size = N)
plot!(plt, direct_sol3)

### Indirect Solve fot the ocp unconstrained

p0 = direct_sol1.costate(0); 

u(q, p) = (p[3] * q[2] - p[2] * q[3]) / p_zero

f = Flow(ocp1, u)

function shoot!(s, p0)
    qqf, pf = f(0, q0, p0, tf)
    s[1:3] .= pf + 2(qqf - qf)
    return nothing
end

s = similar(p0, 3)
shoot!(s, p0)
println("\nNorm of the shooting function: ‖s‖ = ", norm(s), "\n")

#-------------------------USING MINPACK-------------------------------------

backend = AutoForwardDiff();

ξ = p0 # initial guess
nle! = (s, ξ) -> shoot!(s, ξ)
jnle! = (js, ξ) -> jacobian!(nle!, similar(ξ), js, backend, ξ)

indirect_sol = fsolve(nle!, jnle!, ξ; show_trace=true)
p0 = indirect_sol.x

println("")
println("p0 = ", p0)

shoot!(s, p0)
println("\nNorm of the shooting function: ‖s‖ = ", norm(s), "\n")

flow_sol = f((0, tf), q0, p0)
plot!(plt, flow_sol, solution_label="(indirect)")

savefig(plt, "Case2.png")

x = []
y = []
z = []

for i in flow_sol.time_grid
    push!(x,flow_sol.state(i)[1])
    push!(y,flow_sol.state(i)[2])
    push!(z,flow_sol.state(i)[3])
end
plotlyjs()
θ = 0:0.01:π    
φ = 0:0.01:2π  
xs = [sin(t) * cos(p) for t in θ, p in φ]
ys = [sin(t) * sin(p) for t in θ, p in φ]
zs = [cos(t) for t in θ, p in φ]

p = plot(xs, ys, zs, 
    st=:surface, 
    color=:lightblue,
    alpha=0.5, 
    legend=false, 
    axis = nothing, 
    background_color=:transparent,
    grid=false,
)

plot!(x, y, z, lw=2, color=:blue, label="Trajectory")
scatter!([x[1]], [y[1]], [z[1]], markersize=2, color=:green, label="Start")  
scatter!([x[end]], [y[end]], [z[end]], markersize=2, color=:red, label="End")  

savefig(p, "BS2.png")