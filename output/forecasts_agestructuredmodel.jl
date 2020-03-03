push!(LOAD_PATH, joinpath(homedir(),"GitHub/KenyaCoV/src"))
using Plots,Parameters,Distributions,DifferentialEquations,JLD2,DataFrames,StatsPlots,FileIO,MAT,RecursiveArrayTools
import KenyaCoV
using LinearAlgebra:eigen
using Statistics: median, quantile

include("forecast_functions.jl")
gr()#Plotting frontend
"""
Define uncertainty of parameter estimates
"""
d_incubation = LogNormal(log(4.8),0.25) #Liu et al
mean(d_incubation)
(quantile(d_incubation,0.025),median(d_incubation),quantile(d_incubation,0.975))
d_R₀ = Gamma(100,2.92/100) ##Liu et al
mean(d_R₀)
(quantile(d_R₀,0.025),median(d_R₀),quantile(d_R₀,0.975))

"""
Load age structured data
"""

u0,P,P_dest = KenyaCoV.model_ingredients_from_data("data/data_for_age_structuredmodel.jld2",
                                            "data/flight_numbers.csv",
                                            "data/projected_global_prevelance.csv")
u0[KenyaCoV.ind_nairobi_as,5,4] = 5#Five initial infecteds in Nairobi in the 20-24 age group
prob = KenyaCoV.create_KenyaCoV_non_neg_prob(u0,(0.,365.),P)
"""
Declare the treatment/isolation rates considered
"""
treatment_rates = [0.,1/21.,1/14,1/7]


"""
SCENARIO A:
* Age-specific susceptibilties calculated from Chinese case data
* Asymptomatics are not infectious
* 20% of infections are symptomatics
"""
P.dt = 0.25;
P.ext_inf_rate = 0.;
P.ϵ = 0.
P.δ = 0.2
P.γ = 1/2.5
P.σ = 1/rand(d_incubation)
P.β = rand(d_R₀)*P.γ/(P.δ + P.ϵ*(1-P.δ))

results_A = run_scenario(P,prob,1000,treatment_rates)
@save joinpath(homedir(),"Github/KenyaCoVOutputs/results_A.jld2") results_A

"""
SCENARIO B:
* No age-dependent susceptibility
* Asymptomatics are not infectious
* 20% of infections are symptomatics
"""
P.χ = ones(KenyaCoV.n_a)
results_B = run_scenario(P,prob,1000,treatment_rates)
@save joinpath(homedir(),"Github/KenyaCoVOutputs/results_B.jld2") results_B
@save "output/results_B.jld2" results_B

"""
SCENARIO A2:
* Age-specific susceptibilties calculated from Chinese case data
* Asymptomatics are randomly infectious between 0-50% of diseased individuals
* 20% of infections are symptomatics
"""
u0,P,P_dest = KenyaCoV.model_ingredients_from_data("data/data_for_age_structuredmodel.jld2",
                                            "data/flight_numbers.csv",
                                            "data/projected_global_prevelance.csv")
u0[KenyaCoV.ind_nairobi_as,5,4] = 5#Five initial infecteds in Nairobi in the 20-24 age group

P.dt = 0.25;
P.ext_inf_rate = 0.;
P.ϵ = rand(Uniform(0.,0.5))
P.δ = 0.2
P.γ = 1/2.5
P.σ = 1/rand(d_incubation)
P.β = rand(d_R₀)*P.γ/(P.δ + P.ϵ*(1-P.δ))

prob = KenyaCoV.create_KenyaCoV_non_neg_prob(u0,(0.,365.),P)

results_A_rand_inf = run_scenario(P,prob,2000,treatment_rates,randomise_params_and_infectiousness)
@save "output/results_A_rand_inf.jld2" results_A_rand_inf



scatter(C_age[:,1],yerror = (C_age[:,2],C_age[:,3]))
boxplot(peaks[:,:]')
plot(1:365,B[12,:,1],ribbon = (B[12,:,2],B[12,:,3]))
