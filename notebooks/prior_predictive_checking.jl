### A Pluto.jl notebook ###
# v0.12.18

using Markdown
using InteractiveUtils

# ╔═╡ d72ad6f2-5b0a-11eb-2ebd-fb584a8b4696
#Underlying dependencies
using LinearAlgebra,Distributions,Roots,Parameters,DelimitedFiles,Plots

# ╔═╡ efbe0910-5b09-11eb-38b8-854a30129d2e
md"
## Prior predictive checking for KenyaCoV

This notebook describes the **prior predictive checking** done for KenyaCoV. Prior predictive checking is an important part of the workflow of Bayesian data analysis. In short, we interrogate whether the prior assumptions on our parameters align with our prior beliefs about observables from the model. 

Specifically, for KenyaCoV, we are interested in aligning priors over transmission model parameters derived from the literature of transmission modelling outside of Kenya, where parameter uncertainty is usually presented univariately, to a prior belief about the exponential growth rate $r$, which is the most robustly estimated statistic from early-stage epidemic case data. The reasons for performing this analysis is that the prior uncertainty range for underlying transmission parameters (e.g. relative infectiousness of symptomatic cases and age-dependent susceptibility to infection) can give values of $r$ that are unrealistic. This can happen because $r$ is a function of the underlying parameters, and, often, the growth rate can be explained by a higher value of one parameter **and** a lower value of another parameter. Therefore, the prior predictive checking for KenyaCoV aims at taking univariate prior distributions from the literature and reconstructing a multivariate prior which represents our prior beliefs about the possible epidemic growth rates we might observe.

### KenyaCoV model structure for a county
![KenyaCoV3 model structure](https://warwick.ac.uk/fac/cross_fac/zeeman_institute/staffv2/sam_brand/images/model_structure2.jpeg)

The force of infection on a susceptible in age group $a$ depends on the combination of infected types in other age groups:

$\lambda_a = β₀ \frac{\sigma_a}{N_a} \sum_b T_{ab}(\epsilon A_b + P_b + M_b + V_b).$

### Analytic results for $R_0$ and $r$

We define the average (over disease progression) rate of infection from a typical person in age group $b$ to a typical person in age group $a$, at time $\tau$ after infection:

$\beta(a,b,\tau) = {\sigma_a \over N_a}T_{ab}\tilde{\beta}_b(\tau).$
Where $\tilde{\beta}_b(\tau)$ is the average infectiousness of a person in age group $b$ at time $\tau$ after infection. 

Where the conditionality is over whether the infected person in age group $b$ has an asymptomatic disease progression or any other type of disease progression. Standard theory gives:

$K_{ab}(r') := \int_0^\infty e^{-r'\tau} N_a \beta(a,b,\tau) d\tau = \sigma_a T_{ab} \int_0^\infty e^{-r'\tau} \tilde{\beta}_b(\tau)d\tau,$
$R_0 = \rho(K_{ab}(0)),$
$1 = \rho(K_{ab}(r)).$
Where $\rho(K_{ab}(r'))$ is the spectral radius of the matrix. Note that the exponential growth rate $r$ is defined implicitly. The case distribution across ages is the eigenvector associated with the leading eigenvalue of $K_{ab}(0)$.

We define the probablity of a person of age group $b$ being in some disease state $X$ at time $\tau$ after infection as, $P_{X,b}(\tau)$. Therefore,

$\tilde{\beta}_b(\tau) = \beta_0\Big[\epsilon P_{A,b}(\tau) + P_{P,b}(\tau)+P_{M,b}(\tau)+P_{V,b}(\tau)\Big].$

The time evolution of the state probabilities are given by the linear system:

$\partial_t P_{E,b}(\tau) = -\Big((1-\delta_b)\alpha_A + \delta_b \alpha_P \Big)P_{A,b}(\tau),$
$\partial_t P_{A,b}(\tau) = (1-\delta_b)\alpha_A P_{E,b}(\tau) - \gamma_A P_{A,b}(\tau),$
$\partial_t P_{P,b}(\tau) = \delta_b \alpha_P P_{E,b}(\tau) - \alpha_P P_{P,b}(\tau),$
$\partial_t P_{M,b}(\tau) =(1-\upsilon_b) \alpha_P P_{P,b}(\tau) - \gamma_M P_{M,b}(\tau),$
$\partial_t P_{V,b}(\tau) =\upsilon_b \alpha_P P_{P,b}(\tau) - (\gamma_V + \mu_{V,a}) P_{V,b}(\tau).$

With initial condition $P_{E,b}(0) = 1$. Using standard results on the Laplace transforms for dynamical systems we can resolve the integral in the expression of $K_{ab}(r')$:

$\phi_b(r') := \int_0^\infty e^{-r'\tau} \tilde{\beta}_b(\tau)d\tau = {\beta_0\alpha\over \alpha + r'}\Big[{\epsilon (1-\delta_b)\over \gamma_A + r'} + {\delta_b \over \alpha_P+r'}\Big(1 + {(1-\upsilon_b)\alpha_P\over\gamma_M + r'} + {\upsilon_b\alpha_P\over\gamma_V + r'}\Big) \Big].$


"

# ╔═╡ 76b31d52-5b14-11eb-228b-bbcb9cf17fee
md"
### Early growth data that informs our prior

The early growth data that informs our priors is drawn from the early stages of the epidemic in the UK:
* The early exponential growth rate $r$ in the UK corresponded to a doubling time of 3.3 days (r = $(round(log(2)/3.3,digits=3))) (SAGE minutes, 2020).
* The relative attack rate of age groups (>18 year olds) can be estimated from the the REACT2 survey in June (Ward et al, 2020).
"

# ╔═╡ 2507bcba-5b3e-11eb-02f8-df5974baad20
#REACT2 data uniformised by year of life. Missing data is a missing entry
begin
	prop_attack_rate = Union{Float64,Missing}[]
	for a = 1:100
		if a < 18
			push!(prop_attack_rate,missing)
		end
		if a >= 18 && a<25
			push!(prop_attack_rate,0.0786)
		end
		if a >= 25 && a<35
			push!(prop_attack_rate,0.0783)
		end
		if a >= 35 && a<45
			push!(prop_attack_rate,0.0609)
		end
		if a >= 45 && a<55
			push!(prop_attack_rate,0.0641)
		end
		if a >= 55 && a<65
			push!(prop_attack_rate,0.0592)
		end
		if a >= 65 && a<75
			push!(prop_attack_rate,0.0316)
		end
		if a >= 75
			push!(prop_attack_rate,0.0331)
		end		
	end
	prop_attack_rate
end

# ╔═╡ 2a6a1870-5b40-11eb-2a50-83af2c5e752b
begin
	UK_pop_5years_bands = Float64.([3914000,
						3517000,
						3670000,
						3997000,
					4297000,
		4307000,
	4126000,
	4194000,
	4626000,
	4643000,
	4095000,
	3614000,
	3807000,
	3017000,
	2463000,
	2006000,
	1496000,
	918000,
	476000])
	UK_pop_by_year = Float64[]
	for N in UK_pop_5years_bands[1:(end-1)]
		append!(UK_pop_by_year,fill(N/5,5))
	end
	append!(UK_pop_by_year,fill(UK_pop_5years_bands[end]/10,10))
end

# ╔═╡ d4440036-5b40-11eb-17d6-c52e21f9150c
begin
	UK_attack_rate = prop_attack_rate.*UK_pop_by_year
	scatter(UK_attack_rate)
end

# ╔═╡ 16773d86-5b3e-11eb-324b-514d22072f47
md"
Other data we use:

* The pre-measures age-structured mixing matrix for UK (estimated from Prem et al, 2017).
* UK population size by age."

# ╔═╡ 47ee4d2e-5b38-11eb-221e-97edf8e0af66
#We use the transposed from of the mixing matrix from Prem et al 2017
# and normalise the columns
begin
	T=Matrix(readdlm("resources/mixingmatrix_UK_from_to.csv", ',', Float64, '\n';header = true)[1]')
	for j = 1:size(T,2)
		T[:,j] = normalize(T[:,j],1)
	end
end

# ╔═╡ 12c8ed54-5b3e-11eb-256b-871303bc4414


# ╔═╡ e2199394-5b37-11eb-0444-3713fd68d281
begin
	θ = (β₀ = 1.,
		α = 1/3.1,
		ϵ=0.4,
		δ = ones(16),
		υ = fill(0.1,16),
		αP = 1/2,
		γA = 1/9,
		γM = 1/3.1,
		γV = 1/5,
		σ = ones(16))
end

# ╔═╡ 89a91d60-5b36-11eb-0efc-6324bdf17549
#Define the ϕ function described above.
function ϕ(b::Integer,r::Real,params)
	@unpack β₀,α,ϵ,δ,υ,αP,γA,γM,γV = params
	return β₀*(α/(α+r))*((ϵ*(1-δ[b])/(γA + r)) + (δ[b]/(αP + r))*(1 + ((1-υ[b])*αP/(γM + r) ) +( υ[b]*αP /(γV + r) ) ))
end

# ╔═╡ c2c359a8-5b39-11eb-04b2-a9d9ad3b4526
function create_K_matrix(T,params,r::Float64)
	@unpack σ = params
	[σ[a]*T[a,b]*ϕ(b,r,params) for a = 1:length(σ),b = 1:length(σ)]
end

# ╔═╡ f4eb1352-5b3c-11eb-2304-0b9fcf5e7210
K = create_K_matrix(T,θ,0.21)

# ╔═╡ 52f48d46-5b3d-11eb-29d3-21e62f690be8
F= eigen(K)

# ╔═╡ 8c5009b0-5b3a-11eb-056d-cdc747f38701
function find_group_attack_rate(T,params)
	K = create_K_matrix(T,params,0.)
	F = eigen(K)
	return normalize(abs.(F.vectors[:,end]),1)
end

# ╔═╡ c87577f2-5b3c-11eb-0643-01a1cdbcd023
function find_leading_eval_K(T,params,r)
	K = create_K_matrix(T,params,r)
	F = eigen(K)
	leading_eval = Real(F.values[end])
end

# ╔═╡ 6b010218-5b3d-11eb-3e61-57b213348c54
R₀ = find_leading_eval_K(T,θ,0.)

# ╔═╡ 5f401a52-5b3b-11eb-225a-59634b976264
function find_exp_growth_rate(T,params)
	f = r -> find_leading_eval_K(T,params,r) - 1.
	find_zero(f, (0.0001, 1.))
end

# ╔═╡ 8daf5ff6-5b3d-11eb-260b-2f927f507cb8
find_exp_growth_rate(T,θ)

# ╔═╡ Cell order:
# ╠═d72ad6f2-5b0a-11eb-2ebd-fb584a8b4696
# ╟─efbe0910-5b09-11eb-38b8-854a30129d2e
# ╟─76b31d52-5b14-11eb-228b-bbcb9cf17fee
# ╠═2507bcba-5b3e-11eb-02f8-df5974baad20
# ╠═2a6a1870-5b40-11eb-2a50-83af2c5e752b
# ╠═d4440036-5b40-11eb-17d6-c52e21f9150c
# ╟─16773d86-5b3e-11eb-324b-514d22072f47
# ╠═47ee4d2e-5b38-11eb-221e-97edf8e0af66
# ╠═12c8ed54-5b3e-11eb-256b-871303bc4414
# ╠═e2199394-5b37-11eb-0444-3713fd68d281
# ╠═89a91d60-5b36-11eb-0efc-6324bdf17549
# ╠═c2c359a8-5b39-11eb-04b2-a9d9ad3b4526
# ╠═f4eb1352-5b3c-11eb-2304-0b9fcf5e7210
# ╠═52f48d46-5b3d-11eb-29d3-21e62f690be8
# ╠═8c5009b0-5b3a-11eb-056d-cdc747f38701
# ╠═c87577f2-5b3c-11eb-0643-01a1cdbcd023
# ╠═6b010218-5b3d-11eb-3e61-57b213348c54
# ╠═5f401a52-5b3b-11eb-225a-59634b976264
# ╠═8daf5ff6-5b3d-11eb-260b-2f927f507cb8
