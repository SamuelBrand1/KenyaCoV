### A Pluto.jl notebook ###
# v0.12.20

using Markdown
using InteractiveUtils

# ╔═╡ d72ad6f2-5b0a-11eb-2ebd-fb584a8b4696
#Underlying dependencies
using LinearAlgebra,Distributions,Roots,Parameters,DelimitedFiles,StatsPlots,JLD2

# ╔═╡ 1d163474-5c3a-11eb-23f0-53daebc2852c
plotlyjs()

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

First: attach the June REACT2 % prevalence estimate to each year of life.
"

# ╔═╡ 2507bcba-5b3e-11eb-02f8-df5974baad20
#REACT2 data uniformised by year of life. Missing data is a missing entry.
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
	
	scatter(prop_attack_rate,
			xlabel = "Year of life",
			ylabel = "Est. prop. of group infected in UK first wave",
			lab = "",
			title = "REACT2 June Sero-prevalence")
end

# ╔═╡ 83b39838-5bd2-11eb-0ed2-4518b6194b28
md"
Second: Gather the UK population size by year of life. 
"

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
	bar(UK_pop_by_year,
			xlabel = "Year of life",
			ylabel = "N",
			lab = "",
			title = "UK population size")
end

# ╔═╡ 9f84377a-5bd2-11eb-0657-2314b0925317
md"
Third: create an estimated numbers of infected by year of life
"

# ╔═╡ d4440036-5b40-11eb-17d6-c52e21f9150c
begin
	UK_attack_rate = prop_attack_rate.*UK_pop_by_year
	scatter(UK_attack_rate,
			xlabel = "Year of life",
			ylabel = "Est. number of infected in UK first wave",
			lab = "")
end

# ╔═╡ 16773d86-5b3e-11eb-324b-514d22072f47
md"
In order to apply the analytic results outlined above we also use the pre-measures age-structured mixing matrix for UK (estimated from Prem et al, 2017).
"

# ╔═╡ 47ee4d2e-5b38-11eb-221e-97edf8e0af66
#We use the transposed from of the mixing matrix from Prem et al 2017
# and normalise the columns
begin
	T=Matrix(readdlm("resources/mixingmatrix_UK_from_to.csv", ',', Float64, '\n';header = true)[1]')
	heatmap(T,
			xticks = (1:1:16,vcat([string(k)*"-"*string(k+4) for k = 0:5:70],["75+"])),
		yticks = (1:1:16,vcat([string(k)*"-"*string(k+4) for k = 0:5:70],["75+"])),
		title = "Prem et al estimates for age-structured mixing (UK)",
		ylabel = "Age group being contacted",
		xlabel = "Age group of person seeking contact"
	)
end

# ╔═╡ 970828f2-5bd4-11eb-36f5-2b64c6d793df
md"
Its useful to prefactor the `T` matrix so that its leading eigenvalue is one, this gives the $\beta_0$ parameter the interpretation of being $R_0$ for a counter-factual disease that **doesn't have strong age differences in susceptibility/transmissibility**.
"

# ╔═╡ 8bf233da-5bd5-11eb-3e64-31299446486d
begin
	prefactor = Real(eigen(T).values[end])
	rescaled_T = T./prefactor
	heatmap(rescaled_T,
			xticks = (1:1:16,vcat([string(k)*"-"*string(k+4) for k = 0:5:70],["75+"])),
		yticks = (1:1:16,vcat([string(k)*"-"*string(k+4) for k = 0:5:70],["75+"])),
		title = "Rescaled age-structured mixing; prefactor = $(round(prefactor,digits=3))",
		ylabel = "Age group being contacted",
		xlabel = "Age group of person seeking contact"
	)
end

# ╔═╡ 29f311be-5c01-11eb-1036-6953abca540a
md"
It is useful to group the infected number so that they are in 5 year bins and can be put on a scale relative to the 75+ category.
"

# ╔═╡ 43d37af8-5c01-11eb-1b2a-bd4b071568b0
begin	
	UK_rel_attack_rate = [sum(UK_attack_rate[((a-1)*5 + 1):(a*5 )]) for a = 1:15]
	append!(UK_rel_attack_rate,sum(UK_attack_rate[(15*5 + 1):end]))
	UK_rel_attack_rate .= UK_rel_attack_rate./UK_rel_attack_rate[end]
	scatter(UK_rel_attack_rate)
end

# ╔═╡ 3a629928-5bd6-11eb-1c64-4387b5e719da
md"
### Julia functions

We create Julia functions to recreate the analysis outlined above.
"

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

# ╔═╡ 8c5009b0-5b3a-11eb-056d-cdc747f38701
function find_group_rel_attack_rate(T,params)
	K = create_K_matrix(T,params,0.)
	F = eigen(K)
	return abs.(F.vectors[:,end])./abs.(F.vectors[end,end])
end

# ╔═╡ c87577f2-5b3c-11eb-0643-01a1cdbcd023
function find_leading_eval_K(T,params,r)
	K = create_K_matrix(T,params,r)
	F = eigen(K)
	leading_eval = Real(F.values[end])
end

# ╔═╡ 5f401a52-5b3b-11eb-225a-59634b976264
function find_exp_growth_rate(T,params)
	try
		f = r -> find_leading_eval_K(T,params,r) - 1.
		find_zero(f, (0.00001, 1.))
	catch 
		-1.
	end
end

# ╔═╡ d8e23d6e-5bd7-11eb-2a7c-930d706d315d
md"
### Univariate independent priors for each parameter
And gather evidence from the literature to inform a **univariate** independent prior for each unknown parameter.
"

# ╔═╡ 2fa1d2a2-5bdf-11eb-2d80-c119edf7835c
md"
*$\beta_0$ baseline transmissibility.* Given that the age-dependent effects probably depress the true $R_0$ compared to a baseline of no age-dependent effects we have a high prior for $\beta_0$.  
"

# ╔═╡ 40b9e162-5bdd-11eb-3b63-ff5069580271
begin
	d_β₀ = Gamma(5,2.5/5)
	plt_beta = plot(d_β₀,lab="β₀",title = "β₀ prior");
	# nothing
end

# ╔═╡ 55b21b60-5bde-11eb-3a59-f71bbb49f27c
md"
*Latency, incubation, infectious period and presymptomatic transmission.* There are reports of high rates of pre-symptomatic transmission which informs our **assumption** that individuals in the pre-symptomatic phase are as infectious as the symptomatic phase. The % of transmission before and after onset is likely to heavily on aggressiveness of isolation measures.
"

# ╔═╡ c153bbba-5be5-11eb-2d49-35259eb41749
#Parameter priors
begin
	d_meanlatent =Gamma(20,3.1/20)
	d_Pduration =Gamma(15,2/15)
	d_meaninf = Gamma(15,5/15)
	d_meanAinf= Gamma(15,5/15)
	plot(d_meanlatent,lab="mean latent",xlabel = "Mean value of period",ylabel = "Density",title = "Prior distributions for means",legend = :right)
	plot!(d_Pduration,lab="mean pre-symptomatic inf.")
	plot!(d_meaninf,lab="mean post-symptomatic inf.")
end

# ╔═╡ 655ecd36-5bf9-11eb-3595-c51b2d215cf9
md"
*Relative infectiousness of asymptomatic infections, relative symptomatic rate by age, relative susceptibility by age.* We prior each age group by the Polletti et al follow up on contacts study.
"

# ╔═╡ 5244de94-5bfd-11eb-2019-f5c7a45ababe
begin
	dispersion = 0.2
	d_symp_0_19 = Beta(55*dispersion,(304-55)*dispersion)
	d_symp_20_39 = Beta(119*dispersion,(531-119)*dispersion)
	d_symp_40_59 = Beta(306*dispersion,(1002-306)*dispersion)
	d_symp_60_79 = Beta(294*dispersion,(829-294)*dispersion)
	d_symp_80_ = Beta(102*dispersion,(158-102)*dispersion)


	plot(d_symp_0_19,lab= "0-19",xlabel = "prob. symptomatic",ylabel = "Density",
			title = "Priors for age-dependent symptomatic prob.")
	plot!(d_symp_20_39,lab= "20-39")
	plot!(d_symp_40_59,lab= "40-59")
	plot!(d_symp_60_79,lab= "60-79")
	plot!(d_symp_80_,lab= "80+")	
end

# ╔═╡ 2cb525ce-5bff-11eb-2b5e-cb03a5d0b271
begin
	d_ϵ = Beta(10,25-10)
	plot(d_ϵ,lab= "",title = "Prior for relative infectiousness of asymptomatics")
end

# ╔═╡ 96f1e65c-5bff-11eb-3576-45b4c8a8c5da
begin
	d_children_sus = Beta((34/1.47)*dispersion*2,(100-(34/1.47))*dispersion*2)
	d_adults_sus = Beta(1+ 0.5*dispersion*100/1.47,1 + 0.5*dispersion*(100 - 100/1.47) )
	plot(d_children_sus,lab = "0-14",title = "Priors for relative susceptibility compared to 65+",
	xlabel = "rel. sus.",ylabel = "Density",legend = :right)
	plot!(d_adults_sus,lab = "15-64")
end

# ╔═╡ 9cc36794-5c00-11eb-3639-adfb2ece3cb8
md"
### Drawing from the univariate priors and eliminating prior combinations that fail to match observables
Below we list the criteria that observables need to match to be accepted into our joint prior.
"

# ╔═╡ 2c30f1d4-5bdf-11eb-3473-a7582367c79e
incubation_mean_range = [2.5,7.5]

# ╔═╡ 692218d8-5be5-11eb-0e04-d560bba66e95
presymptomatic_transmission_prop_range = [0.4,0.7]

# ╔═╡ 880c2fe6-5c02-11eb-1be4-9f8418aa859e
r_range = [log(2)/5,log(2)/3]

# ╔═╡ 9a02fa9a-5c02-11eb-16e1-e1d83982e27f
mean_abs_err_attack_range = 2.

# ╔═╡ 3466631e-5c04-11eb-3d18-bfe5e088b3df
md"
We create a function for drawing from the univariate priors.
"

# ╔═╡ b5ccc8ca-5c02-11eb-1412-33123caa8c4a
function draw_parameter_set()
	δ = sort!([rand(d_symp_0_19) for a = 1:4])
	append!(δ,sort!([rand(d_symp_20_39) for a = 5:8]))
	append!(δ,sort!([rand(d_symp_40_59) for a = 9:12]))
	append!(δ,sort!([rand(d_symp_60_79) for a = 13:15]))
	append!(δ,rand(d_symp_80_))
	
	σ =  sort!([rand(d_children_sus) for a = 1:3])
	append!(σ, sort!([rand(d_adults_sus) for a = 4:15]))
	append!(σ,1.)	
	# sort!(σ)

	θ = (β₀ = rand(d_β₀),
		α = 1/rand(d_meanlatent),
		ϵ= rand(d_ϵ),
		δ = δ,
		υ = zeros(16),
		αP = 1/rand(d_Pduration),
		γA = 1/rand(d_meanAinf),
		γM = 1/rand(d_meaninf),
		γV = 1/5,
		σ = σ)
	return θ
end

# ╔═╡ 5ae8f506-5c04-11eb-0927-3752b906b091
md"
And a function for accepting or rejecting the prior set drawn.
"

# ╔═╡ 7110b0da-5c04-11eb-2126-69f99575a96b
θ₀ = draw_parameter_set()

# ╔═╡ e88d5f64-5c04-11eb-04ee-1d965ffb2837
function accept_reject_params(θ)
	r = find_exp_growth_rate(rescaled_T,θ)
	v = find_group_rel_attack_rate(rescaled_T,θ)
	mean_abs_err = mean(abs.(v[5:end] .- UK_rel_attack_rate[5:end]))
	isinincubrange = (1/θ.α)  + (1/θ.αP) >= incubation_mean_range[1] && (1/θ.α)  + (1/θ.αP) <= incubation_mean_range[2] 
	ispresymp_range = (1/θ.αP)/((1/θ.αP) + (1/θ.γM))  >= presymptomatic_transmission_prop_range[1] && (1/θ.αP)/((1/θ.αP) + (1/θ.γM))  <= presymptomatic_transmission_prop_range[2] 
	isr_range = r >= r_range[1] && r <= r_range[2]
	is_abs_err_range = mean_abs_err <= mean_abs_err_attack_range
	return isinincubrange && ispresymp_range && isr_range && is_abs_err_range
end
	

# ╔═╡ 652a4a14-5c05-11eb-30b7-7fb982654d0c
sum([accept_reject_params(draw_parameter_set()) for k = 1:1000])

# ╔═╡ c1e48092-5c07-11eb-293d-2da52426657d
begin
	accepted_parameter_set = []
	for k = 1:10
		θ = draw_parameter_set()
		if accept_reject_params(θ)
			append!(accepted_parameter_set,[θ])
		end
	end
end

# ╔═╡ 712f46d8-5c42-11eb-3cd7-8547406ef6c9
begin
	accepted_parameter_set_matrix = zeros(length(accepted_parameter_set),6+16+16)
	for i = 1:size(accepted_parameter_set_matrix,1)
		accepted_parameter_set_matrix[i,:] .= vcat(accepted_parameter_set[i].β₀,
													accepted_parameter_set[i].α,
													accepted_parameter_set[i].ϵ,
													accepted_parameter_set[i].αP,
													accepted_parameter_set[i].γA,
													accepted_parameter_set[i].γM,
													accepted_parameter_set[i].δ,
													accepted_parameter_set[i].σ)
	end
	accepted_parameter_set_matrix = Matrix(accepted_parameter_set_matrix')
end	

# ╔═╡ 63e45682-5c45-11eb-20d2-372938f4eefc
D = accepted_parameter_set_matrix[1:37,:]#Remove the row of ones for rel. susceptibility of 75+ yos

# ╔═╡ c5608da2-5c44-11eb-03be-a3b3b35672da
det(log.(D)*log.(D)')

# ╔═╡ f0472210-5c42-11eb-241a-191e0011b084
open("drawn_priors.csv", "w") do io
           writedlm(io, accepted_parameter_set_matrix,',')
       end

# ╔═╡ 4c695f7a-5c45-11eb-3805-a5878cd30ab2
d_log_priors = fit_mle(MvNormal,log.(D))

# ╔═╡ 92151488-5c43-11eb-1259-5fb43d168eda
heatmap(d_log_priors.Σ)

# ╔═╡ d8b3d354-5c45-11eb-1007-61a768b1240e
JLD2.@save("d_log_priors.jld2",d_log_priors)

# ╔═╡ 14d88b7a-5c46-11eb-1327-319ee8250de7
JLD2.@save("accepted_parameters.jld2",accepted_parameter_set)

# ╔═╡ 37f2604a-5c46-11eb-0cc8-2d96ffb62b71
r_prior = [find_exp_growth_rate(rescaled_T,θ) for θ in accepted_parameter_set]

# ╔═╡ 629cae5e-5c46-11eb-2875-d786b01c1c64
histogram(log(2)./r_prior,norm = :pdf,bins = 50)

# ╔═╡ 9b270fda-5c46-11eb-1408-9fbefe624912
v_prior = [find_group_rel_attack_rate(rescaled_T,θ) for θ in accepted_parameter_set]

# ╔═╡ c02f7e8c-5c46-11eb-0343-cf09ac7c3584
begin
bar(mean(v_prior),yerr = 3*std(v_prior))
scatter!(UK_rel_attack_rate)
end


# ╔═╡ Cell order:
# ╠═d72ad6f2-5b0a-11eb-2ebd-fb584a8b4696
# ╠═1d163474-5c3a-11eb-23f0-53daebc2852c
# ╟─efbe0910-5b09-11eb-38b8-854a30129d2e
# ╟─76b31d52-5b14-11eb-228b-bbcb9cf17fee
# ╟─2507bcba-5b3e-11eb-02f8-df5974baad20
# ╟─83b39838-5bd2-11eb-0ed2-4518b6194b28
# ╟─2a6a1870-5b40-11eb-2a50-83af2c5e752b
# ╟─9f84377a-5bd2-11eb-0657-2314b0925317
# ╠═d4440036-5b40-11eb-17d6-c52e21f9150c
# ╟─16773d86-5b3e-11eb-324b-514d22072f47
# ╟─47ee4d2e-5b38-11eb-221e-97edf8e0af66
# ╟─970828f2-5bd4-11eb-36f5-2b64c6d793df
# ╟─8bf233da-5bd5-11eb-3e64-31299446486d
# ╟─29f311be-5c01-11eb-1036-6953abca540a
# ╠═43d37af8-5c01-11eb-1b2a-bd4b071568b0
# ╟─3a629928-5bd6-11eb-1c64-4387b5e719da
# ╟─89a91d60-5b36-11eb-0efc-6324bdf17549
# ╟─c2c359a8-5b39-11eb-04b2-a9d9ad3b4526
# ╟─8c5009b0-5b3a-11eb-056d-cdc747f38701
# ╟─c87577f2-5b3c-11eb-0643-01a1cdbcd023
# ╟─5f401a52-5b3b-11eb-225a-59634b976264
# ╟─d8e23d6e-5bd7-11eb-2a7c-930d706d315d
# ╟─2fa1d2a2-5bdf-11eb-2d80-c119edf7835c
# ╠═40b9e162-5bdd-11eb-3b63-ff5069580271
# ╟─55b21b60-5bde-11eb-3a59-f71bbb49f27c
# ╠═c153bbba-5be5-11eb-2d49-35259eb41749
# ╟─655ecd36-5bf9-11eb-3595-c51b2d215cf9
# ╠═5244de94-5bfd-11eb-2019-f5c7a45ababe
# ╠═2cb525ce-5bff-11eb-2b5e-cb03a5d0b271
# ╠═96f1e65c-5bff-11eb-3576-45b4c8a8c5da
# ╠═9cc36794-5c00-11eb-3639-adfb2ece3cb8
# ╠═2c30f1d4-5bdf-11eb-3473-a7582367c79e
# ╠═692218d8-5be5-11eb-0e04-d560bba66e95
# ╠═880c2fe6-5c02-11eb-1be4-9f8418aa859e
# ╠═9a02fa9a-5c02-11eb-16e1-e1d83982e27f
# ╠═3466631e-5c04-11eb-3d18-bfe5e088b3df
# ╠═b5ccc8ca-5c02-11eb-1412-33123caa8c4a
# ╟─5ae8f506-5c04-11eb-0927-3752b906b091
# ╠═7110b0da-5c04-11eb-2126-69f99575a96b
# ╠═e88d5f64-5c04-11eb-04ee-1d965ffb2837
# ╠═652a4a14-5c05-11eb-30b7-7fb982654d0c
# ╠═c1e48092-5c07-11eb-293d-2da52426657d
# ╠═712f46d8-5c42-11eb-3cd7-8547406ef6c9
# ╠═63e45682-5c45-11eb-20d2-372938f4eefc
# ╠═c5608da2-5c44-11eb-03be-a3b3b35672da
# ╠═f0472210-5c42-11eb-241a-191e0011b084
# ╠═4c695f7a-5c45-11eb-3805-a5878cd30ab2
# ╠═92151488-5c43-11eb-1259-5fb43d168eda
# ╠═d8b3d354-5c45-11eb-1007-61a768b1240e
# ╠═14d88b7a-5c46-11eb-1327-319ee8250de7
# ╠═37f2604a-5c46-11eb-0cc8-2d96ffb62b71
# ╠═629cae5e-5c46-11eb-2875-d786b01c1c64
# ╠═9b270fda-5c46-11eb-1408-9fbefe624912
# ╠═c02f7e8c-5c46-11eb-0343-cf09ac7c3584
