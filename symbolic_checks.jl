### A Pluto.jl notebook ###
# v0.19.25

using Markdown
using InteractiveUtils

# ╔═╡ 411f41f8-e676-4bec-b092-61af470bab1b
begin
	using Symbolics
	using LinearAlgebra
	# using PlutoTest
	using BenchmarkTools
	
	import Pipe.@pipe

	# import Symbolics.derivative as sderiv
	sderiv(args...) = Symbolics.derivative(args...; simplify=true)
end

# ╔═╡ 2a1df9be-6c52-456b-bf13-375cf871b5b6
md"#### Benchmarks"

# ╔═╡ 9b5553df-254f-4c26-9e7f-10d2af1c25a6
md"""#### Goals
 1. Check that $F_R, F_θ, F_ϕ$ matches the definition on Wikipedia
 1. Compute $\vec F_\mathrm{res}$
 1. Compute $\overrightarrow{F_\mathrm{res}} ⋅ \partial \overrightarrow{OK} / \partial R$ and $\overrightarrow{F_\mathrm{res}} ⋅ \partial \overrightarrow{OK} / \partial \tau$
 1. Find how to compute the mass matrix and the right-hand side
"""

# ╔═╡ 730092d1-4d6f-491e-bba4-3cb3df8eb91a
let
	md"""
	#### Mapple sheet

	```math
	\begin{align}
	\boldsymbol{F}_\mathrm{res} &= m \frac{\mathrm{d^2} \mathbf{OK({\color{red}R, \tau}/{\color{green}\alpha, \theta_2, \phi_2})}}{\mathrm{d}t^2} \\
	F_\mathrm{tract} &= \frac{I_g \ddot{\alpha} - C_g}{l \sin \theta_2 \sin \phi_2} \\
	\boldsymbol{F}_\mathrm{res} &= \boldsymbol{F}_\mathrm{app} - \boldsymbol{F}_\mathrm{grav} - \boldsymbol{F}_\mathrm{aero} + \boldsymbol{F}_\mathrm{Tract}
	\end{align}
	```
	
	```math
	\begin{cases}
	\boldsymbol{F}_\mathrm{res} \cdot \frac{\mathrm{d} \mathbf{OK({\color{red}R, \tau}/{\color{green}\alpha, \theta_2, \phi_2})}}{\mathrm{d}\tau} &= 0
	\quad (\ddot \alpha, {\color{red}\ddot R, \ddot \tau}/{\color{green}\ddot \theta_2, \ddot \phi_2}) \\
	\boldsymbol{F}_\mathrm{res} \cdot \frac{\mathrm{d} \mathbf{OK({\color{red}R, \tau}/{\color{green}\alpha, \theta_2, \phi_2})}}{\mathrm{d}R} &= 0
	\quad (\ddot \alpha, {\color{red}\ddot R, \ddot \tau}/{\color{green}\ddot \theta_2, \ddot \phi_2}) \\
	\left( \frac{\mathrm{d^2} \mathbf{OK}(R, \tau)}{\mathrm{d}t^2} - \frac{\mathrm{d^2} \mathbf{OK}(\alpha, \theta_2, \phi_2)}{\mathrm{d}t^2} \right) \cdot \hat x &= 0
	\quad (\ddot \alpha, \ddot \theta_2, \ddot \phi_2, \ddot R, \ddot \tau) \\
	\left( \frac{\mathrm{d^2} \mathbf{OK}(R, \tau)}{\mathrm{d}t^2} - \frac{\mathrm{d^2} \mathbf{OK}(\alpha, \theta_2, \phi_2)}{\mathrm{d}t^2} \right) \cdot \hat y &= 0
	\quad (\ddot \alpha, \ddot \theta_2, \ddot \phi_2, \ddot R, \ddot \tau) \\
	\left( \frac{\mathrm{d^2} \mathbf{OK}(R, \tau)}{\mathrm{d}t^2} - \frac{\mathrm{d^2} \mathbf{OK}(\alpha, \theta_2, \phi_2)}{\mathrm{d}t^2} \right) \cdot \hat z &= 0
	\quad (\ddot \theta_2, \ddot R, \ddot \tau) \\
	\end{cases}
	```
	```math
	\boldsymbol{X} = (\ddot \alpha, \ddot \theta_2, \ddot \phi_2, \ddot R, \ddot \tau)
	```

	-> Trouver l'équvialent de `GenerateMatrix` de chez Maple, ou bien résoudre ce système linéaire sans passer par la forme matricielle ?
	"""
end

# ╔═╡ 1fe7f71a-6dc4-4d01-91bb-76605d8173f5
begin
	@variables t
	@variables θ2(t) ϕ2(t) α(t) τ(t) R(t)
	@variables m l r Ig θ0 ϕ0 Δθ Δϕ
	@variables Cg
	# @variables Fres[1:3] Fgrav[1:3] Faero[1:3] Ftract[1:3]
end;

# ╔═╡ f8b2147d-08b1-481d-8787-df33517caef5
begin
	compute_θϕ(τ, θ0, ϕ0, Δθ, Δϕ) = θ0 + Δθ * sin(2τ), ϕ0 + Δϕ * sin(τ)
	θ, ϕ = compute_θϕ(τ, θ0, ϕ0, Δθ, Δϕ)

	# @variables θ(τ, θ0, ϕ0, Δθ, Δϕ) ϕ(τ, θ0, ϕ0, Δθ, Δϕ)

	# @variables θ(τ) ϕ(τ)
end;

# ╔═╡ 8df9c6c6-9bec-43c3-937a-200f9edd9ebc
md"""
$\vec \Omega = \vec x \times \vec v / \lVert \vec x \rVert^2$
"""

# ╔═╡ cac0f793-514f-41b4-86e1-a82d34b8677f
begin
	d(y) = Differential(t)(y)
	d(y, x) = sderiv(y, x)
	dd(y) = d(d(y))
	dd(y, x) = d(d(y, x), x)
	dd(y, x1, x2) = d(d(y, x1), x2)
end

# ╔═╡ 82bd8198-9b8b-4d9c-ac01-6c6420fc6c3a
"""Evaluates an expression with deterministic random values given to each variable.
Already vectorized, calling `reval(M)` where `M` is a matrix works.

`reval.(norm([sin(θ)*cos(ϕ) sin(θ)*sin(ϕ) cos(θ)]))` gives 1 whereas `norm([sin(θ)*cos(ϕ) sin(θ)*sin(ϕ) cos(θ)])` is not able to simplify to 1.

TODO: Give derivatives a numeric value aswell.
"""
function reval(expr; seed=0)
	function reval_(expr; seed)
		random_value(var) = @pipe (var, seed, ℯ*π) .|> hash |> xor(_...) |> float |> _ * 2^-64
		
		generate_derivatives(var) = [var, d(var), dd(var), sderiv(var, τ), sderiv(var, R)]
		subs = Dict([(derivative, random_value(var)) for var in Symbolics.get_variables(expr) for derivative in generate_derivatives(var)])
		
		# subs = Dict([(var, random_value(var)) for var in Symbolics.get_variables(expr)])
		
		substitute(expand.(expand_derivatives.(expr)), subs)
	end
	reval_.(expr; seed=seed)
end

# ╔═╡ 064c857e-82df-44ce-b215-6d540da22a99
begin
	@variables θ_ dθ_ ϕ_ dϕ_
	substitutions = Dict(
		θ0 + Δθ * sin(2τ) => θ_,
		2Δθ * cos(2τ) * d(τ) => dθ_,
		ϕ0 + Δϕ * sin(τ) => ϕ_,
		Δϕ * cos(τ) * d(τ) => dϕ_,
	)
end

# ╔═╡ 301eea73-5cf6-47c8-9c23-13e3ec2c2250
begin
	Fr = m * (r * d(θ2)^2 - r * (d(α) + d(θ2))^2 * sin(θ2)^2 + dd(α) * l * sin(θ2) * sin(ϕ2) - d(α)^2 * l * sin(θ2) * cos(ϕ2))
	Fθ2 = m * (r * dd(θ2) - d(α)^2 * l * cos(θ2) * cos(ϕ2) + dd(α) * l * cos(θ2) * sin(ϕ2) - (d(α) + d(ϕ2))^2 * r * sin(θ2) * cos(θ2))
	Fϕ2 = m * (r * dd(ϕ2) * sin(θ2) + (2 * d(r) * d(ϕ2) + dd(α) * r) * sin(θ2) + 2 * (d(α + ϕ2) * d(θ2) * r * cos(θ2) + dd(α) * l * cos(ϕ2) + d(α)^2 * l * sin(ϕ2)))
end

# ╔═╡ abefb1f4-710b-4c36-82b1-33b0cc165b10
begin
	function compute_pos1(R, θ, ϕ)
	    [R * sin(θ) * cos(ϕ); R * sin(θ) * sin(ϕ); R * cos(θ)]
	end
	function compute_pos2(α, θ2, ϕ2, r, l)
	    [l * cos(α) + r * sin(θ2) * cos(ϕ2 + α); l * sin(α) + r * sin(θ2) * sin(ϕ2 + α); r * cos(θ2)]
	end
end

# ╔═╡ 1f7011aa-a5fa-44ae-ad96-4fff36e18b2b
begin
	OK = compute_pos1(R, θ, ϕ)
	Rhat = sderiv(OK, R)
	θhat = sderiv(Rhat, θ)
	ϕhat = Rhat × θhat
	rot = hcat(Rhat, θhat, ϕhat)'; rot = simplify.(rot)
end;

# ╔═╡ 79900870-0c3c-4437-abca-6021cc1e459a
rot .|> reval |> eachcol .|> norm

# ╔═╡ 82ee007e-9971-4b6f-af1d-605cd6f99b9d
begin
	OK2 = compute_pos2(α, θ2, ϕ2, r, l)
	rhat = sderiv(OK2, r)
	θ2hat = sderiv(rhat, θ2)
	ϕ2hat = rhat × θ2hat
	rot2 = hcat(rhat, θ2hat, ϕ2hat)'; rot2 = simplify.(rot2)
end;

# ╔═╡ 9f713155-2971-4817-af95-fd51042f4878
(sderiv(rhat, t) - sderiv.(rhat, [θ2 ϕ2 α]) * sderiv([θ2; ϕ2; α], t))

# ╔═╡ c8a73f8e-5db8-4c4f-afba-da275534466b
function divvec(v1, v2)
	# return (q, r) ∈ R×R^n such that v1 = q * v2 + r
	q = v1 ⋅ v2 / norm(v2)^2
	r = v1 - q*v2
	return q, r
end

# ╔═╡ e990b110-f51f-42e0-aec7-7a4f056558f6
function sincosvec2(v1, v2)
	# equivalent to
	# ```julia
	# q = v1 ⋅ v2 / norm(v2)^2
	# r = v1 - q*v2
	# sin_v1v2 = norm(r) / norm(v1)
	# cos_v1v2 = q * norm(v2) / norm(v1)
	# ```
	n1, n2 = norm(v1), norm(v2)
	cos_v1v2 = v1 ⋅ v2 / (n1 * n2)
	sin_v1v2 = norm(v1 / n1 - cos_v1v2 * v2 / n2)
	return sin_v1v2, cos_v1v2
end

# ╔═╡ 331602c8-890e-4480-9984-931897b12dec
rhat, θ2hat, ϕ2hat

# ╔═╡ 39ea2f24-d70d-4227-8463-3e8abf45f26c
rot12 = [cos(α) sin(α) 0; -sin(α) cos(α) 0; 0 0 1]

# ╔═╡ 8912c20b-5b78-45ff-bbcf-c32d390f0c28
rot23 = [sin(θ2)*cos(ϕ2) sin(θ2)*sin(ϕ2) cos(θ2); cos(θ2)*cos(ϕ2) cos(θ2)*sin(ϕ2) -sin(θ2); -sin(ϕ2) cos(ϕ2) 0]

# ╔═╡ 3965398d-72fa-456a-839e-581f23b3afa9
rot2 .|> reval |> eachcol .|> norm

# ╔═╡ 2798d6c2-c97d-4806-be42-2327efce59e0
begin
	@variables Fgrav_r Fgrav_θ2
	@variables Faero_r Faero_θ2 Faero_ϕ2
	# @variables Ftrac_r
	Ftrac_r = (Ig*dd(α) - Cg)/(l*sin(θ2)*sin(ϕ2))
	# @variables Fapp[1:3]
	Fapp = [Fr; Fθ2; Fϕ2]
	Fgrav = rot2 * [Fgrav_r; Fgrav_θ2; 0]
	Faero = rot2 * [Faero_r; Faero_θ2; Faero_ϕ2]
	Ftrac = rot2 * [Ftrac_r; 0; 0]
	Fres = Fapp - Fgrav - Faero + Ftrac
	# @variables Fres[1:3]
end;

# ╔═╡ 2ab708ed-3394-4a56-b024-ed3104db32e9
begin
	M = [0 0 m*r*sin(θ2)*((-sin(ϕ2)*cos(α)-cos(ϕ2)*sin(α))*(-Δϕ*cos(τ)*sin(ϕ)*sin(θ)+2*Δθ*cos(2*τ)*cos(ϕ)*cos(θ))+(-sin(ϕ2)*sin(α)+cos(ϕ2)*cos(α))*(Δϕ*cos(τ)*cos(ϕ)*sin(θ)+2*Δθ*cos(2*τ)*sin(ϕ)*cos(θ))) m*r*((cos(θ2)*cos(ϕ2)*cos(α)-cos(θ2)*sin(ϕ2)*sin(α))*(-Δϕ*cos(τ)*sin(ϕ)*sin(θ)+2*Δθ*cos(2*τ)*cos(ϕ)*cos(θ))+(cos(θ2)*cos(ϕ2)*sin(α)+cos(θ2)*sin(ϕ2)*cos(α))*(Δϕ*cos(τ)*cos(ϕ)*sin(θ)+2*Δθ*cos(2*τ)*sin(ϕ)*cos(θ))+2*sin(θ2)*Δθ*cos(2*τ)*sin(θ)) (m*l*sin(θ2)*sin(ϕ2)+Ig/l/sin(θ2)/sin(ϕ2))*((sin(θ2)*cos(ϕ2)*cos(α)-sin(θ2)*sin(ϕ2)*sin(α))*(-Δϕ*cos(τ)*sin(ϕ)*sin(θ)+2*Δθ*cos(2*τ)*cos(ϕ)*cos(θ))+(sin(θ2)*cos(ϕ2)*sin(α)+sin(θ2)*sin(ϕ2)*cos(α))*(Δϕ*cos(τ)*cos(ϕ)*sin(θ)+2*Δθ*cos(2*τ)*sin(ϕ)*cos(θ))-2*cos(θ2)*Δθ*cos(2*τ)*sin(θ))+m*l*cos(θ2)*sin(ϕ2)*((cos(θ2)*cos(ϕ2)*cos(α)-cos(θ2)*sin(ϕ2)*sin(α))*(-Δϕ*cos(τ)*sin(ϕ)*sin(θ)+2*Δθ*cos(2*τ)*cos(ϕ)*cos(θ))+(cos(θ2)*cos(ϕ2)*sin(α)+cos(θ2)*sin(ϕ2)*cos(α))*(Δϕ*cos(τ)*cos(ϕ)*sin(θ)+2*Δθ*cos(2*τ)*sin(ϕ)*cos(θ))+2*sin(θ2)*Δθ*cos(2*τ)*sin(θ))+m*(r*sin(θ2)+l*cos(ϕ2))*((-sin(ϕ2)*cos(α)-cos(ϕ2)*sin(α))*(-Δϕ*cos(τ)*sin(ϕ)*sin(θ)+2*Δθ*cos(2*τ)*cos(ϕ)*cos(θ))+(-sin(ϕ2)*sin(α)+cos(ϕ2)*cos(α))*(Δϕ*cos(τ)*cos(ϕ)*sin(θ)+2*Δθ*cos(2*τ)*sin(ϕ)*cos(θ)));0 0 m*r*sin(θ2)*((-sin(ϕ2)*cos(α)-cos(ϕ2)*sin(α))*cos(ϕ)*sin(θ)+(-sin(ϕ2)*sin(α)+cos(ϕ2)*cos(α))*sin(ϕ)*sin(θ)) m*r*((cos(θ2)*cos(ϕ2)*cos(α)-cos(θ2)*sin(ϕ2)*sin(α))*cos(ϕ)*sin(θ)+(cos(θ2)*cos(ϕ2)*sin(α)+cos(θ2)*sin(ϕ2)*cos(α))*sin(ϕ)*sin(θ)-sin(θ2)*cos(θ)) (m*l*sin(θ2)*sin(ϕ2)+Ig/l/sin(θ2)/sin(ϕ2))*((sin(θ2)*cos(ϕ2)*cos(α)-sin(θ2)*sin(ϕ2)*sin(α))*cos(ϕ)*sin(θ)+(sin(θ2)*cos(ϕ2)*sin(α)+sin(θ2)*sin(ϕ2)*cos(α))*sin(ϕ)*sin(θ)+cos(θ2)*cos(θ))+m*l*cos(θ2)*sin(ϕ2)*((cos(θ2)*cos(ϕ2)*cos(α)-cos(θ2)*sin(ϕ2)*sin(α))*cos(ϕ)*sin(θ)+(cos(θ2)*cos(ϕ2)*sin(α)+cos(θ2)*sin(ϕ2)*cos(α))*sin(ϕ)*sin(θ)-sin(θ2)*cos(θ))+m*(r*sin(θ2)+l*cos(ϕ2))*((-sin(ϕ2)*cos(α)-cos(ϕ2)*sin(α))*cos(ϕ)*sin(θ)+(-sin(ϕ2)*sin(α)+cos(ϕ2)*cos(α))*sin(ϕ)*sin(θ));-2*R*cos(2*τ)*Δθ*sin(θ) cos(θ) 0 r*sin(θ2) 0;-R*cos(τ)*Δϕ*sin(ϕ)*sin(θ)+2*R*cos(2*τ)*Δθ*cos(ϕ)*cos(θ) cos(ϕ)*sin(θ) r*sin(θ2)*sin(ϕ2+α) -r*cos(θ2)*cos(ϕ2+α) r*sin(θ2)*sin(ϕ2+α)+l*sin(α);R*cos(τ)*Δϕ*cos(ϕ)*sin(θ)+2*R*cos(2*τ)*Δθ*sin(ϕ)*cos(θ) sin(ϕ)*sin(θ) -r*sin(θ2)*cos(ϕ2+α) -r*cos(θ2)*sin(ϕ2+α) -r*sin(θ2)*cos(ϕ2+α)-l*cos(α)]

	# d([^\*\+\-\^...]) -> d($1)
	# Faero_kite[([1-3])] - Faero_kite[$1] -> Faero[$1]
	sm = [
		-(m*(-r*d(θ2)^2-r*(d(α)+d(ϕ2))^2*sin(θ2)^2-d(α)^2*l*sin(θ2)*cos(ϕ2))-Fgrav[1]-Faero[1]-Cg/l/sin(θ2)/sin(ϕ2))*((sin(θ2)*cos(ϕ2)*cos(α)-sin(θ2)*sin(ϕ2)*sin(α))*(-Δϕ*cos(τ)*sin(ϕ)*sin(θ)+2*Δθ*cos(2*τ)*cos(ϕ)*cos(θ))+(sin(θ2)*cos(ϕ2)*sin(α)+sin(θ2)*sin(ϕ2)*cos(α))*(Δϕ*cos(τ)*cos(ϕ)*sin(θ)+2*Δθ*cos(2*τ)*sin(ϕ)*cos(θ))-2*cos(θ2)*Δθ*cos(2*τ)*sin(θ))-(m*(-d(α)^2*l*cos(θ2)*cos(ϕ2)-(d(α)+d(ϕ2))^2*r*sin(θ2)*cos(θ2))-Fgrav[2]-Faero[2])*((cos(θ2)*cos(ϕ2)*cos(α)-cos(θ2)*sin(ϕ2)*sin(α))*(-Δϕ*cos(τ)*sin(ϕ)*sin(θ)+2*Δθ*cos(2*τ)*cos(ϕ)*cos(θ))+(cos(θ2)*cos(ϕ2)*sin(α)+cos(θ2)*sin(ϕ2)*cos(α))*(Δϕ*cos(τ)*cos(ϕ)*sin(θ)+2*Δθ*cos(2*τ)*sin(ϕ)*cos(θ))+2*sin(θ2)*Δθ*cos(2*τ)*sin(θ))-(m*(2*(d(α)+d(ϕ2))*d(θ2)*r*cos(θ2)+d(α)^2*l*sin(ϕ2))-Faero[3])*((-sin(ϕ2)*cos(α)-cos(ϕ2)*sin(α))*(-Δϕ*cos(τ)*sin(ϕ)*sin(θ)+2*Δθ*cos(2*τ)*cos(ϕ)*cos(θ))+(-sin(ϕ2)*sin(α)+cos(ϕ2)*cos(α))*(Δϕ*cos(τ)*cos(ϕ)*sin(θ)+2*Δθ*cos(2*τ)*sin(ϕ)*cos(θ)))
		-(m*(-r*d(θ2)^2-r*(d(α)+d(ϕ2))^2*sin(θ2)^2-d(α)^2*l*sin(θ2)*cos(ϕ2))-Fgrav[1]-Faero[1]-Cg/l/sin(θ2)/sin(ϕ2))*((sin(θ2)*cos(ϕ2)*cos(α)-sin(θ2)*sin(ϕ2)*sin(α))*cos(ϕ)*sin(θ)+(sin(θ2)*cos(ϕ2)*sin(α)+sin(θ2)*sin(ϕ2)*cos(α))*sin(ϕ)*sin(θ)+cos(θ2)*cos(θ))-(m*(-d(α)^2*l*cos(θ2)*cos(ϕ2)-(d(α)+d(ϕ2))^2*r*sin(θ2)*cos(θ2))-Fgrav[2]-Faero[2])*((cos(θ2)*cos(ϕ2)*cos(α)-cos(θ2)*sin(ϕ2)*sin(α))*cos(ϕ)*sin(θ)+(cos(θ2)*cos(ϕ2)*sin(α)+cos(θ2)*sin(ϕ2)*cos(α))*sin(ϕ)*sin(θ)-sin(θ2)*cos(θ))-(m*(2*(d(α)+d(ϕ2))*d(θ2)*r*cos(θ2)+d(α)^2*l*sin(ϕ2))-Faero[3])*((-sin(ϕ2)*cos(α)-cos(ϕ2)*sin(α))*cos(ϕ)*sin(θ)+(-sin(ϕ2)*sin(α)+cos(ϕ2)*cos(α))*sin(ϕ)*sin(θ))
		4*R*d(τ)^2*Δθ^2*cos(2*τ)^2*cos(θ)-4*R*d(τ)^2*Δθ*sin(2*τ)*sin(θ)+4*d(R)*d(τ)*Δθ*cos(2*τ)*sin(θ)-r*d(θ2)^2*cos(θ2)
		2*d(R)*d(τ)*Δϕ*cos(τ)*sin(ϕ)*sin(θ)-4*d(R)*d(τ)*Δθ*cos(2*τ)*cos(ϕ)*cos(θ)+4*R*d(τ)^2*Δϕ*cos(τ)*Δθ*cos(2*τ)*sin(ϕ)*cos(θ)+R*(d(τ)^2*Δϕ^2*cos(τ)^2+4*d(τ)^2*Δθ^2*cos(2*τ)^2)*cos(ϕ)*sin(θ)-R*d(τ)^2*Δϕ*sin(τ)*sin(ϕ)*sin(θ)+4*R*d(τ)^2*Δθ*sin(2*τ)*cos(ϕ)*cos(θ)-l*d(α)^2*cos(α)-r*d(θ2)^2*sin(θ2)*cos(ϕ2+α)-r*(d(α)+d(ϕ2))^2*sin(θ2)*cos(ϕ2+α)-2*r*d(θ2)*(d(α)+d(ϕ2))*cos(θ2)*sin(ϕ2+α)
		-2*d(R)*d(τ)*Δϕ*cos(τ)*cos(ϕ)*sin(θ)-4*d(R)*d(τ)*Δθ*cos(2*τ)*sin(ϕ)*cos(θ)-4*R*d(τ)^2*Δϕ*cos(τ)*Δθ*cos(2*τ)*cos(ϕ)*cos(θ)+R*(d(τ)^2*Δϕ^2*cos(τ)^2+4*d(τ)^2*Δθ^2*cos(2*τ)^2)*sin(ϕ)*sin(θ)+R*d(τ)^2*Δϕ*sin(τ)*cos(ϕ)*sin(θ)+4*R*d(τ)^2*Δθ*sin(2*τ)*sin(ϕ)*cos(θ)-l*d(α)^2*sin(α)-r*d(θ2)^2*sin(θ2)*sin(ϕ2+α)-r*(d(α)+d(ϕ2))^2*sin(θ2)*sin(ϕ2+α)+2*r*d(θ2)*(d(α)+d(ϕ2))*cos(θ2)*cos(ϕ2+α)]
end;

# ╔═╡ a4881f7e-0b94-4a7a-8bd4-019d4fe4c503
size(M), size(sm)

# ╔═╡ 4f780255-ee2e-4166-bd76-a8ae70af3555
reval(M)

# ╔═╡ 09873caa-5286-41b4-84fc-27d270864f53
let
	dθ = 2*Δθ*cos(2*τ)
	dϕ = Δϕ*cos(τ)
	e1 = ((cos(θ2)*cos(ϕ2)*cos(α)-cos(θ2)*sin(ϕ2)*sin(α))*(-dϕ*sin(ϕ)*sin(θ)+dθ*cos(ϕ)*cos(θ))+(cos(θ2)*cos(ϕ2)*sin(α)+cos(θ2)*sin(ϕ2)*cos(α))*(dϕ*cos(ϕ)*sin(θ)+dθ*sin(ϕ)*cos(θ))+sin(θ2)*dθ*sin(θ))
	
	[0 0 m*r*sin(θ2)*((-sin(ϕ2)*cos(α)-cos(ϕ2)*sin(α))*(-dϕ*sin(ϕ)*sin(θ)+dθ*cos(ϕ)*cos(θ))+(-sin(ϕ2)*sin(α)+cos(ϕ2)*cos(α))*(dϕ*cos(ϕ)*sin(θ)+dθ*sin(ϕ)*cos(θ))) m*r*e1 (m*l*sin(θ2)*sin(ϕ2)+Ig/l/sin(θ2)/sin(ϕ2))*((sin(θ2)*cos(ϕ2)*cos(α)-sin(θ2)*sin(ϕ2)*sin(α))*(-dϕ*sin(ϕ)*sin(θ)+dθ*cos(ϕ)*cos(θ))+(sin(θ2)*cos(ϕ2)*sin(α)+sin(θ2)*sin(ϕ2)*cos(α))*(dϕ*cos(ϕ)*sin(θ)+dθ*sin(ϕ)*cos(θ))-cos(θ2)*dθ*sin(θ))+m*l*cos(θ2)*sin(ϕ2)*e1+m*(r*sin(θ2)+l*cos(ϕ2))*((-sin(ϕ2)*cos(α)-cos(ϕ2)*sin(α))*(-dϕ*sin(ϕ)*sin(θ)+dθ*cos(ϕ)*cos(θ))+(-sin(ϕ2)*sin(α)+cos(ϕ2)*cos(α))*(dϕ*cos(ϕ)*sin(θ)+dθ*sin(ϕ)*cos(θ)));0 0 m*r*sin(θ2)*((-sin(ϕ2)*cos(α)-cos(ϕ2)*sin(α))*cos(ϕ)*sin(θ)+(-sin(ϕ2)*sin(α)+cos(ϕ2)*cos(α))*sin(ϕ)*sin(θ)) m*r*((cos(θ2)*cos(ϕ2)*cos(α)-cos(θ2)*sin(ϕ2)*sin(α))*cos(ϕ)*sin(θ)+(cos(θ2)*cos(ϕ2)*sin(α)+cos(θ2)*sin(ϕ2)*cos(α))*sin(ϕ)*sin(θ)-sin(θ2)*cos(θ)) (m*l*sin(θ2)*sin(ϕ2)+Ig/l/sin(θ2)/sin(ϕ2))*((sin(θ2)*cos(ϕ2)*cos(α)-sin(θ2)*sin(ϕ2)*sin(α))*cos(ϕ)*sin(θ)+(sin(θ2)*cos(ϕ2)*sin(α)+sin(θ2)*sin(ϕ2)*cos(α))*sin(ϕ)*sin(θ)+cos(θ2)*cos(θ))+m*l*cos(θ2)*sin(ϕ2)*((cos(θ2)*cos(ϕ2)*cos(α)-cos(θ2)*sin(ϕ2)*sin(α))*cos(ϕ)*sin(θ)+(cos(θ2)*cos(ϕ2)*sin(α)+cos(θ2)*sin(ϕ2)*cos(α))*sin(ϕ)*sin(θ)-sin(θ2)*cos(θ))+m*(r*sin(θ2)+l*cos(ϕ2))*((-sin(ϕ2)*cos(α)-cos(ϕ2)*sin(α))*cos(ϕ)*sin(θ)+(-sin(ϕ2)*sin(α)+cos(ϕ2)*cos(α))*sin(ϕ)*sin(θ));-R*dθ*sin(θ) cos(θ) 0 r*sin(θ2) 0;-R*dϕ*sin(ϕ)*sin(θ)+R*dθ*cos(ϕ)*cos(θ) cos(ϕ)*sin(θ) r*sin(θ2)*sin(ϕ2+α) -r*cos(θ2)*cos(ϕ2+α) r*sin(θ2)*sin(ϕ2+α)+l*sin(α);R*dϕ*cos(ϕ)*sin(θ)+dθ*R*sin(ϕ)*cos(θ) sin(ϕ)*sin(θ) -r*sin(θ2)*cos(ϕ2+α) -r*cos(θ2)*sin(ϕ2+α) -r*sin(θ2)*cos(ϕ2+α)-l*cos(α)] - M
end

# ╔═╡ 2ee4357e-72a6-4a4e-8f19-f02a083c3625
Fgrav

# ╔═╡ 4f2495ca-37df-46a0-831f-a3748ae3a30e
begin
	∂OK∂R = sderiv(OK, R)
	∂OK∂τ = sderiv(OK, τ)
	eq1 = Fres ⋅ ∂OK∂R ~ 0
	eq2 = Fres ⋅ ∂OK∂τ ~ 0
end;

# ╔═╡ 368a77e7-ea58-4fe2-be08-a6a4cbb64aa5
begin
	# sderiv(sderiv(OK - OK2, t), t) - expand_derivatives.(d.(d.(OK - OK2))) == [0 0 0]
	eq3, eq4, eq5 = sderiv(sderiv(OK - OK2, t), t) .~ 0
end;

# ╔═╡ d230f7cf-03c4-469f-8cf5-75f4e4b7c4c6
begin
	eqs = [eq1, eq2, eq3, eq4, eq5][[2, 1, 5, 3, 4]]
	vars = [dd(α), dd(θ2), dd(ϕ2), dd(R), dd(τ)][[5, 4, 3, 2, 1]]
end;

# ╔═╡ 8d11ff9b-f241-42f0-b8cc-4b150a8bf408
# Symbolics.solve_for(eqs, vars)

# ╔═╡ 13ff6044-88f8-45aa-9707-9d49d2284886
begin
	var_to_zero = Dict(var => 0 for var in vars)
	M_ = [sderiv(eq.lhs, var) for eq in eqs, var in vars]
	sm_ = [eq.lhs - M_coeffs ⋅ vars for (eq, M_coeffs) in zip(eqs, eachrow(M))]
	sm_ = [substitute(-eq.lhs, var_to_zero) for (eq, M_coeffs) in zip(eqs, eachrow(M))]
	sm_ = convert(Vector{Num}, sm_)
end;

# ╔═╡ ea75501c-f118-4e99-9f9b-6721e608d012
length.(string.(collect(M - M_)))

# ╔═╡ db0105b7-15d2-4a3f-92f4-e1c200ab8b16
print(simplify(substitute(M[2, 3], substitutions)))

# ╔═╡ c0a61504-cbdb-4fc0-80f7-878cd765333e
Fres

# ╔═╡ 1068dc15-fbae-4b40-98ef-cf509bc7701c
Fres ⋅ ∂OK∂τ

# ╔═╡ 668eb2fa-3225-499f-8150-4672ab3b92d8
length.(string.(collect(expand.(expand_derivatives.(sm - sm_)))))

# ╔═╡ 0d1ff03b-faf6-42ae-bf03-adddca84c84a
print(expand.(expand_derivatives.((sm - sm_)[[1, 2]])))

# ╔═╡ ba51358a-fb61-4007-9df5-7d6a436102bf
print(expand.(expand_derivatives.((M - M_)[[1, 2], :])))

# ╔═╡ ecc5ca13-f914-45b3-8c1a-476c2bac8ca1
print(expand.(expand_derivatives.((M - M_)[2, 3])))

# ╔═╡ 906585b1-f8ff-4e98-b704-d32d03ecd296
print(expand.(expand_derivatives.((M_)[2, 3])))

# ╔═╡ 0286b162-c9e2-4012-8d8b-8f1751ce7449
print(expand.(expand_derivatives.((M)[2, 3])))

# ╔═╡ 5e18f1c8-1b13-4246-b6c1-5f1d4da3096b
expand.(expand_derivatives.((sm - sm_)[[3, 4, 5]]))

# ╔═╡ 515e54b5-bf6f-4fbb-bf2d-f47085f2f00e
expand.(expand_derivatives.((M - M_)[[3, 4, 5],:]))

# ╔═╡ 7b6ab57a-0e6e-4308-9035-21d50562cddb
reval(α + ϕ2)

# ╔═╡ 9516f0b7-b111-445d-90e0-335d29d0e135
(d(α + ϕ) - d(α) - d(ϕ))

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
BenchmarkTools = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Pipe = "b98c9c47-44ae-5843-9183-064241ee97a0"
Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"

[compat]
BenchmarkTools = "~1.3.2"
Pipe = "~1.3.0"
Symbolics = "~5.3.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.0"
manifest_format = "2.0"
project_hash = "5b72692d0a9ec784916137766939bcdd3656fcef"

[[deps.AbstractAlgebra]]
deps = ["GroupsCore", "InteractiveUtils", "LinearAlgebra", "MacroTools", "Random", "RandomExtensions", "SparseArrays", "Test"]
git-tree-sha1 = "3ee5c58774f4487a5bf2bb05e39d91ff5022b4cc"
uuid = "c3fe647b-3220-5bb0-a1ea-a7954cac585d"
version = "0.29.4"

[[deps.AbstractTrees]]
git-tree-sha1 = "faa260e4cb5aba097a73fab382dd4b5819d8ec8c"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.4"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "cc37d689f599e8df4f464b2fa3870ff7db7492ef"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.6.1"
weakdeps = ["StaticArrays"]

    [deps.Adapt.extensions]
    AdaptStaticArraysExt = "StaticArrays"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra", "Requires", "SnoopPrecompile", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "38911c7737e123b28182d89027f4216cfc8a9da7"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.4.3"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BenchmarkTools]]
deps = ["JSON", "Logging", "Printf", "Profile", "Statistics", "UUIDs"]
git-tree-sha1 = "d9a9701b899b30332bbcb3e1679c41cce81fb0e8"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.3.2"

[[deps.Bijections]]
git-tree-sha1 = "fe4f8c5ee7f76f2198d5c2a06d3961c249cce7bd"
uuid = "e2ed5e7c-b2de-5872-ae92-c73ca462fb04"
version = "0.1.4"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "e30f2f4e20f7f186dc36529910beaedc60cfa644"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.16.0"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonSolve]]
git-tree-sha1 = "9441451ee712d1aec22edad62db1a9af3dc8d852"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.3"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["UUIDs"]
git-tree-sha1 = "7a60c856b9fa189eb34f5f8a6f6b5529b7942957"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.6.1"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.2+0"

[[deps.CompositeTypes]]
git-tree-sha1 = "02d2316b7ffceff992f3096ae48c7829a8aa0638"
uuid = "b152e2b5-7a66-4b01-a709-34e65c35f657"
version = "0.1.3"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "738fec4d684a9a6ee9598a8bfee305b26831f28c"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.2"
weakdeps = ["IntervalSets", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.DataAPI]]
git-tree-sha1 = "e8119c1a33d267e16108be441a287a6981ba1630"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.14.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "a4ad7ef19d2cdc2eff57abbbe68032b1cd0bd8f8"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.13.0"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "eead66061583b6807652281c0fbf291d7a9dc497"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.90"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.DomainSets]]
deps = ["CompositeTypes", "IntervalSets", "LinearAlgebra", "Random", "StaticArrays", "Statistics"]
git-tree-sha1 = "698124109da77b6914f64edd696be8dccf90229e"
uuid = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
version = "0.6.6"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[deps.DynamicPolynomials]]
deps = ["DataStructures", "Future", "LinearAlgebra", "MultivariatePolynomials", "MutableArithmetics", "Pkg", "Reexport", "Test"]
git-tree-sha1 = "8b84876e31fa39479050e2d3395c4b3b210db8b0"
uuid = "7c1d4256-1411-5781-91ec-d7bc3513ac07"
version = "0.4.6"

[[deps.EnumX]]
git-tree-sha1 = "bdb1942cd4c45e3c678fd11569d5cccd80976237"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.4"

[[deps.ExprTools]]
git-tree-sha1 = "c1d06d129da9f55715c6c212866f5b1bddc5fa00"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.9"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "fc86b4fd3eff76c3ce4f5e96e2fdfa6282722885"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.0.0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "00e252f4d706b3d55a8863432e742bf5717b498d"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.35"
weakdeps = ["StaticArrays"]

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

[[deps.FunctionWrappers]]
git-tree-sha1 = "d62485945ce5ae9c0c48f124a84998d755bae00e"
uuid = "069b7b12-0de2-55c6-9aab-29f3d0a68a2e"
version = "1.1.3"

[[deps.FunctionWrappersWrappers]]
deps = ["FunctionWrappers"]
git-tree-sha1 = "b104d487b34566608f8b4e1c39fb0b10aa279ff8"
uuid = "77dc65aa-8811-40c2-897b-53d922fa7daf"
version = "0.1.3"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "1cd7f0af1aa58abc02ea1d872953a97359cb87fa"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.1.4"

[[deps.Groebner]]
deps = ["AbstractAlgebra", "Combinatorics", "Logging", "MultivariatePolynomials", "Primes", "Random", "SnoopPrecompile"]
git-tree-sha1 = "b6c3e9e1eb8dcc6fd9bc68fe08dcc7ab22710de6"
uuid = "0b43b601-686d-58a3-8a1c-6623616c7cd4"
version = "0.3.4"

[[deps.GroupsCore]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "9e1a5e9f3b81ad6a5c613d181664a0efc6fe6dd7"
uuid = "d5909c97-4eac-4ecc-a3dc-fdd0858a4120"
version = "0.4.0"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "84204eae2dd237500835990bcade263e27674a93"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.16"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.IntegerMathUtils]]
git-tree-sha1 = "f366daebdfb079fd1fe4e3d560f99a0c892e15bc"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IntervalSets]]
deps = ["Dates", "Random", "Statistics"]
git-tree-sha1 = "16c0cc91853084cb5f58a78bd209513900206ce6"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.4"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.LabelledArrays]]
deps = ["ArrayInterface", "ChainRulesCore", "ForwardDiff", "LinearAlgebra", "MacroTools", "PreallocationTools", "RecursiveArrayTools", "StaticArrays"]
git-tree-sha1 = "cd04158424635efd05ff38d5f55843397b7416a9"
uuid = "2ee39098-c373-598a-b85f-a56591580800"
version = "1.14.0"

[[deps.LambertW]]
git-tree-sha1 = "c5ffc834de5d61d00d2b0e18c96267cffc21f648"
uuid = "984bce1d-4616-540c-a9ee-88d1112d94c9"
version = "0.4.6"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Printf", "Requires"]
git-tree-sha1 = "099e356f267354f46ba65087981a77da23a279b7"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.0"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

[[deps.Lazy]]
deps = ["MacroTools"]
git-tree-sha1 = "1370f8202dac30758f3c345f9909b97f53d87d3f"
uuid = "50d2b5c4-7a5e-59d5-8109-a42b560f39c0"
version = "0.15.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "0a1b7c2863e44523180fdb3146534e265a91870b"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.23"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "42324d08725e200c23d4dfb549e0d5d89dede2d2"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.10"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+0"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.10.11"

[[deps.MultivariatePolynomials]]
deps = ["ChainRulesCore", "DataStructures", "LinearAlgebra", "MutableArithmetics"]
git-tree-sha1 = "eaa98afe2033ffc0629f9d0d83961d66a021dfcc"
uuid = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
version = "0.4.7"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "964cb1a7069723727025ae295408747a0b36a854"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.3.0"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.21+4"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "d321bf2de576bf25ec4d3e4360faca399afca282"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.0"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "67eae2738d63117a196f497d7db789821bce61d1"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.17"

[[deps.Parsers]]
deps = ["Dates", "SnoopPrecompile"]
git-tree-sha1 = "478ac6c952fddd4399e71d4779797c538d0ff2bf"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.8"

[[deps.Pipe]]
git-tree-sha1 = "6842804e7867b115ca9de748a0cf6b364523c16d"
uuid = "b98c9c47-44ae-5843-9183-064241ee97a0"
version = "1.3.0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.9.0"

[[deps.PreallocationTools]]
deps = ["Adapt", "ArrayInterface", "ForwardDiff", "Requires"]
git-tree-sha1 = "f739b1b3cc7b9949af3b35089931f2b58c289163"
uuid = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
version = "0.4.12"

    [deps.PreallocationTools.extensions]
    PreallocationToolsReverseDiffExt = "ReverseDiff"

    [deps.PreallocationTools.weakdeps]
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "259e206946c293698122f63e2b513a7c99a244e8"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.1.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "7eb1686b4f04b82f96ed7a4ea5890a4f0c7a09f1"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.0"

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "311a2aa90a64076ea0fac2ad7492e914e6feeb81"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Profile]]
deps = ["Printf"]
uuid = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "6ec7ac8412e83d57e313393220879ede1740f9ee"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.8.2"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RandomExtensions]]
deps = ["Random", "SparseArrays"]
git-tree-sha1 = "062986376ce6d394b23d5d90f01d81426113a3c9"
uuid = "fb686558-2515-59ef-acaa-46db3789a887"
version = "0.4.3"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecursiveArrayTools]]
deps = ["Adapt", "ArrayInterface", "DocStringExtensions", "GPUArraysCore", "IteratorInterfaceExtensions", "LinearAlgebra", "RecipesBase", "Requires", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables"]
git-tree-sha1 = "68078e9fa9130a6a768815c48002d0921a232c11"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "2.38.4"

    [deps.RecursiveArrayTools.extensions]
    RecursiveArrayToolsMeasurementsExt = "Measurements"
    RecursiveArrayToolsTrackerExt = "Tracker"
    RecursiveArrayToolsZygoteExt = "Zygote"

    [deps.RecursiveArrayTools.weakdeps]
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "f65dcb5fa46aee0cf9ed6274ccbd597adc49aa7b"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.1"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6ed52fdd3382cf21947b15e8870ac0ddbff736da"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.4.0+0"

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "d7d9ebe28062161c1e314ed643097b0c6fe657d9"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.7"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SciMLBase]]
deps = ["ArrayInterface", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "EnumX", "FunctionWrappersWrappers", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "Preferences", "RecipesBase", "RecursiveArrayTools", "Reexport", "RuntimeGeneratedFunctions", "SciMLOperators", "SnoopPrecompile", "StaticArraysCore", "Statistics", "SymbolicIndexingInterface", "Tables", "TruncatedStacktraces"]
git-tree-sha1 = "392d3e28b05984496af37100ded94dc46fa6c8de"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "1.91.7"

[[deps.SciMLOperators]]
deps = ["ArrayInterface", "DocStringExtensions", "Lazy", "LinearAlgebra", "Setfield", "SparseArrays", "StaticArraysCore", "Tricks"]
git-tree-sha1 = "d4b254372887d7cbad20ef34021cd69cb0e1d05f"
uuid = "c0aeaf25-5076-4817-a8d5-81caf7dfa961"
version = "0.2.4"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.SnoopPrecompile]]
deps = ["Preferences"]
git-tree-sha1 = "e760a70afdcd461cf01a575947738d359234665c"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "a4ada03f999bd01b3a25dcaa30b2d929fe537e00"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.1.0"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "ef28127915f4229c971eb43f3fc075dd3fe91880"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.2.0"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "c262c8e978048c2b095be1672c9bee55b4619521"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.5.24"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6b7ba252635a5eff6a0b0664a41ee140a1c9e72a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.9.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "45a7769a04a3cf80da1c1c7c60caf932e6f4c9f7"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.6.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "75ebe04c5bed70b91614d684259b661c9e6274a4"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.0"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "f625d686d5a88bcd2b15cd81f18f98186fdc0c9a"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.3.0"

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

    [deps.StatsFuns.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "5.10.1+6"

[[deps.SymbolicIndexingInterface]]
deps = ["DocStringExtensions"]
git-tree-sha1 = "f8ab052bfcbdb9b48fad2c80c873aa0d0344dfe5"
uuid = "2efcf032-c050-4f8e-a9bb-153293bab1f5"
version = "0.2.2"

[[deps.SymbolicUtils]]
deps = ["AbstractTrees", "Bijections", "ChainRulesCore", "Combinatorics", "ConstructionBase", "DataStructures", "DocStringExtensions", "DynamicPolynomials", "IfElse", "LabelledArrays", "LinearAlgebra", "MultivariatePolynomials", "NaNMath", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "TimerOutputs", "Unityper"]
git-tree-sha1 = "5cb1f963f82e7b81305102dd69472fcd3e0e1483"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "1.0.5"

[[deps.Symbolics]]
deps = ["ArrayInterface", "ConstructionBase", "DataStructures", "DiffRules", "Distributions", "DocStringExtensions", "DomainSets", "Groebner", "IfElse", "LaTeXStrings", "LambertW", "Latexify", "Libdl", "LinearAlgebra", "MacroTools", "Markdown", "NaNMath", "RecipesBase", "Reexport", "Requires", "RuntimeGeneratedFunctions", "SciMLBase", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicUtils", "TreeViews"]
git-tree-sha1 = "e23ec62c083ca8f15a4b7174331b3b8d1c511e47"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "5.3.1"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "1544b926975372da01227b382066ab70e574a3ec"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.10.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "f548a9e9c490030e545f72074a41edfd0e5bcdd7"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.23"

[[deps.TreeViews]]
deps = ["Test"]
git-tree-sha1 = "8d0d7a3fe2f30d6a7f833a5f19f7c7a5b396eae6"
uuid = "a2a6695c-b41b-5b7d-aed9-dbfdeacea5d7"
version = "0.3.0"

[[deps.Tricks]]
git-tree-sha1 = "aadb748be58b492045b4f56166b5188aa63ce549"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.7"

[[deps.TruncatedStacktraces]]
deps = ["InteractiveUtils", "MacroTools", "Preferences"]
git-tree-sha1 = "7bc1632a4eafbe9bd94cf1a784a9a4eb5e040a91"
uuid = "781d530d-4396-4725-bb49-402e4bee1e77"
version = "1.3.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Unityper]]
deps = ["ConstructionBase"]
git-tree-sha1 = "d5f4ec8c22db63bd3ccb239f640e895cfde145aa"
uuid = "a7c27f48-0311-42f6-a7f8-2c11e75eb415"
version = "0.1.2"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.7.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"
"""

# ╔═╡ Cell order:
# ╟─2a1df9be-6c52-456b-bf13-375cf871b5b6
# ╟─9b5553df-254f-4c26-9e7f-10d2af1c25a6
# ╟─730092d1-4d6f-491e-bba4-3cb3df8eb91a
# ╠═1fe7f71a-6dc4-4d01-91bb-76605d8173f5
# ╟─82bd8198-9b8b-4d9c-ac01-6c6420fc6c3a
# ╠═f8b2147d-08b1-481d-8787-df33517caef5
# ╟─2ab708ed-3394-4a56-b024-ed3104db32e9
# ╠═a4881f7e-0b94-4a7a-8bd4-019d4fe4c503
# ╠═4f780255-ee2e-4166-bd76-a8ae70af3555
# ╠═09873caa-5286-41b4-84fc-27d270864f53
# ╟─8df9c6c6-9bec-43c3-937a-200f9edd9ebc
# ╠═064c857e-82df-44ce-b215-6d540da22a99
# ╠═cac0f793-514f-41b4-86e1-a82d34b8677f
# ╠═301eea73-5cf6-47c8-9c23-13e3ec2c2250
# ╠═9f713155-2971-4817-af95-fd51042f4878
# ╠═abefb1f4-710b-4c36-82b1-33b0cc165b10
# ╠═1f7011aa-a5fa-44ae-ad96-4fff36e18b2b
# ╠═79900870-0c3c-4437-abca-6021cc1e459a
# ╠═82ee007e-9971-4b6f-af1d-605cd6f99b9d
# ╠═c8a73f8e-5db8-4c4f-afba-da275534466b
# ╠═e990b110-f51f-42e0-aec7-7a4f056558f6
# ╠═331602c8-890e-4480-9984-931897b12dec
# ╠═39ea2f24-d70d-4227-8463-3e8abf45f26c
# ╠═8912c20b-5b78-45ff-bbcf-c32d390f0c28
# ╠═3965398d-72fa-456a-839e-581f23b3afa9
# ╠═2798d6c2-c97d-4806-be42-2327efce59e0
# ╠═2ee4357e-72a6-4a4e-8f19-f02a083c3625
# ╠═4f2495ca-37df-46a0-831f-a3748ae3a30e
# ╠═368a77e7-ea58-4fe2-be08-a6a4cbb64aa5
# ╠═d230f7cf-03c4-469f-8cf5-75f4e4b7c4c6
# ╠═8d11ff9b-f241-42f0-b8cc-4b150a8bf408
# ╠═13ff6044-88f8-45aa-9707-9d49d2284886
# ╠═ea75501c-f118-4e99-9f9b-6721e608d012
# ╠═db0105b7-15d2-4a3f-92f4-e1c200ab8b16
# ╠═c0a61504-cbdb-4fc0-80f7-878cd765333e
# ╠═1068dc15-fbae-4b40-98ef-cf509bc7701c
# ╠═668eb2fa-3225-499f-8150-4672ab3b92d8
# ╠═0d1ff03b-faf6-42ae-bf03-adddca84c84a
# ╠═ba51358a-fb61-4007-9df5-7d6a436102bf
# ╠═ecc5ca13-f914-45b3-8c1a-476c2bac8ca1
# ╠═906585b1-f8ff-4e98-b704-d32d03ecd296
# ╠═0286b162-c9e2-4012-8d8b-8f1751ce7449
# ╠═5e18f1c8-1b13-4246-b6c1-5f1d4da3096b
# ╠═515e54b5-bf6f-4fbb-bf2d-f47085f2f00e
# ╠═7b6ab57a-0e6e-4308-9035-21d50562cddb
# ╠═9516f0b7-b111-445d-90e0-335d29d0e135
# ╠═411f41f8-e676-4bec-b092-61af470bab1b
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
