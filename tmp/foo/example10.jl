# Example 10 from
# Ansel, Q. et al. Introduction to theoretical and experimental aspects of quantum optimal control,
# J. of Phys. B: Atomic, Molecular and Optical Physics 57.13 (2024), 133001. doi: 10.1088/1361-6455/ad46a5

using OptimalControl  
using NLPModelsIpopt  
using OrdinaryDiffEq  
using MINPACK      
using Plots           
using LinearAlgebra: norm
using DifferentiationInterface
import ForwardDiff

# -----------------  Direct solve  --------------------

q0 = [1.0, 0.0, 0.0] 
qf = [0.0, 1.0, 0.0]  

ocp = @def begin
    tf ∈ R, variable
    t ∈ [0, tf], time
    q = [x, y, z] ∈ R³, state
    u = [ux, uy] ∈ R², control

    q(0)  == q0 
    q(tf) == qf 
    ux(t)^2 + uy(t)^2 ≤ 1 
    q̇(t) == [ uy(t) * z(t)
             -ux(t) * z(t)
             -uy(t) * x(t) + ux(t) * y(t) ]
    1 ≤ tf ≤ 25 

    tf → min
end

N = 100
direct_sol = solve(ocp; grid_size = N, linear_solver = "mumps")
plt = plot(direct_sol; size = (1000,1000), background_color=:transparent)
p0 = direct_sol.costate(0); np0 = norm(p0); p0 = p0 / np0
tf = direct_sol.time_grid[end]

# -----------------  Indirect solve  --------------------

function u(q, p) 
    H1 = p[3] * q[2] - p[2] * q[3]
    H2 = p[1] * q[3] - p[3] * q[1]
    n = sqrt(H1^2 + H2^2)
    return [H1, H2] / n
end

f = Flow(ocp, (q, p, tf) -> u(q, p))

function shoot!(s, p0, tf)
    qqf, pf = f(0, q0, p0, tf)
    s[1:3] = qqf - qf
    s[4] = p0' * p0 - 1
    return nothing
end

s = similar(p0, 4);
shoot!(s, p0, tf);
println("\nNorm of the shooting function: ‖s‖ = ", norm(s), "\n")

#-------------------------USING MINPACK-------------------------------------
# exemple reference: https://control-toolbox.org/OptimalControl.jl/stable/tutorial-iss.html

backend = AutoForwardDiff();

ξ = [p0 ; tf] # initial guess
nle! = (s, ξ) -> shoot!(s, ξ[1:3], ξ[4])
jnle! = (js, ξ) -> jacobian!(nle!, similar(ξ), js, backend, ξ)

indirect_sol = fsolve(nle!, jnle!, ξ; show_trace=true)
p0 = indirect_sol.x[1:3];
tf = indirect_sol.x[4];

println("")
println("p0 = ", p0)
println("tf = ", tf)

shoot!(s, p0, tf);
println("\nNorm of the shooting function: ‖s‖ = ", norm(s), "\n")

flow_sol = f((0, tf), q0, np0 * p0);         
plot!(plt, flow_sol, solution_label="(indirect)")

p0_sol = [0.0, √3 / 3, 1.0]; p0_sol = p0_sol / norm(p0_sol)
println("p0 error:", p0 - p0_sol)
println("tf error:", tf - π * √3 / 2)

savefig(plt, "Case1.png")

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

savefig(p, "BS1.png")