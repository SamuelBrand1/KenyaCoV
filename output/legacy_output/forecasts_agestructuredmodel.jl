push!(LOAD_PATH, joinpath(homedir(),"GitHub/KenyaCoV/src"))
using Plots,Parameters,Distributions,DifferentialEquations,JLD2,DataFrames,StatsPlots,FileIO,MAT,RecursiveArrayTools
import KenyaCoV
using LinearAlgebra:eigen
using Statistics: median, quantile

gr()#Plotting frontend


"""
Load age structured data, define initial state and declare the KenyaCoV problem
"""

u0,P,P_dest = KenyaCoV.model_ingredients_from_data("data/data_for_age_structuredmodel.jld2",
                                            "data/flight_numbers.csv",
                                            "data/projected_global_prevelance.csv")
@load "data/agemixingmatrix_china.jld2" M_China
@load "data/reporting_rate_for_MERS_like_scenario.jld2" reporting_rate

u0[KenyaCoV.ind_nairobi_as,5,4] = 5#Five initial infecteds in Nairobi in the 20-24 age group
prob = KenyaCoV.create_KenyaCoV_non_neg_prob(u0,(0.,365.),P)
"""
Declare the treatment/isolation rates considered, and a callback that limits intervention
"""
treatment_rates = [(0.,1),(0.,0.75),(0.,0.5),(1/7.,1.),(1/7,0.75),(1/7,0.5),(1/3.5,1.),(1/3.5,0.75),(1/3.5,0.5)]

function isolation_limit(u,t,integrator) # Isolation can only continue until 10000 detected cases
  integrator.p.isolating_detecteds && sum(u[:,:,8]) > 1e4
end
function affect_isolation_limit!(integrator)
  integrator.p.τ = 0
  integrator.p.isolating_detecteds = false
end
cb_iso_limit = DiscreteCallback(isolation_limit,affect_isolation_limit!)




"""
SCENARIO A1:
* Age-specific susceptibilties calculated from Chinese case data
* Asymptomatics are not infectious
* 20% of infections are symptomatics
"""
P.dt = 0.25;
P.ext_inf_rate = 0.;
P.ϵ = 0.
P.δ = 0.2
P.γ = 1/2.5
P.σ = 1/rand(KenyaCoV.d_incubation)
P.β = rand(KenyaCoV.d_R₀)*P.γ/(P.δ + P.ϵ*(1-P.δ))

results_A1 = KenyaCoV.run_scenario(P,prob,10,treatment_rates,cb_iso_limit)

@save joinpath(homedir(),"Github/KenyaCoVOutputs/results_A1.jld2") results_A1

"""
SCENARIO B1:
* No age-specific susceptibilties --- its disease difference
* Asymptomatics are 50% as infectious as symptomatics
* 20% of infections are symptomatics
"""
u0,P,P_dest = KenyaCoV.model_ingredients_from_data("data/data_for_age_structuredmodel.jld2",
                                            "data/flight_numbers.csv",
                                            "data/projected_global_prevelance.csv")
#Redistribute susceptibility JUST to rescale infectiousness so we get the correct R₀/r
P.χ = ones(KenyaCoV.n_a)
sus_matrix = repeat(P.χ,1,KenyaCoV.n_a)
eigs, = eigen(sus_matrix.*M_China)
R₀_scale = Real(eigs[end])
P.χ = ones(KenyaCoV.n_a)/R₀_scale

P.dt = 0.25;
P.ext_inf_rate = 0.;
P.ϵ = 0.5
P.δ = 0.2
P.γ = 1/2.5
P.σ = 1/rand(KenyaCoV.d_incubation)
P.β = rand(KenyaCoV.d_R₀)*P.γ/(P.δ + P.ϵ*(1-P.δ))
results_B1 = KenyaCoV.run_scenario(P,prob,1000,treatment_rates,cb_iso_limit)
@save joinpath(homedir(),"Github/KenyaCoVOutputs/results_B1.jld2") results_B1
results_B1 = 0



"""
SCENARIO A2:
* No age-dependent susceptibility
* Asymptomatics are not infectious
* 80% of infections are symptomatics
"""
u0,P,P_dest = KenyaCoV.model_ingredients_from_data("data/data_for_age_structuredmodel.jld2",
                                            "data/flight_numbers.csv",
                                            "data/projected_global_prevelance.csv")

P.dt = 0.25;
P.ext_inf_rate = 0.;
P.ϵ = 0.
P.δ = 0.8
P.γ = 1/2.5
P.σ = 1/rand(KenyaCoV.d_incubation)
P.β = rand(KenyaCoV.d_R₀)*P.γ/(P.δ + P.ϵ*(1-P.δ))

results_A2 = KenyaCoV.run_scenario(P,prob,1000,treatment_rates,cb_iso_limit)
@save joinpath(homedir(),"Github/KenyaCoVOutputs/results_A2.jld2") results_A2


"""
SCENARIO B2:
* No age-specific susceptibilties --- its disease difference
* Asymptomatics are 50% as infectious as symptomatics
* 80% of infections are symptomatics
"""
u0,P,P_dest = KenyaCoV.model_ingredients_from_data("data/data_for_age_structuredmodel.jld2",
                                            "data/flight_numbers.csv",
                                            "data/projected_global_prevelance.csv")
#Redistribute susceptibility JUST to rescale infectiousness so we get the correct R₀/r
P.χ = ones(KenyaCoV.n_a)
sus_matrix = repeat(P.χ,1,KenyaCoV.n_a)
eigs, = eigen((P.β/P.γ)*sus_matrix.*P.M)
R₀_scale = Real(eigs[end])
P.χ = ones(KenyaCoV.n_a)/R₀_scale

P.dt = 0.25;
P.ext_inf_rate = 0.;
P.ϵ = 0.5
P.δ = 0.8
P.γ = 1/2.5
P.σ = 1/rand(KenyaCoV.d_incubation)
P.β = rand(KenyaCoV.d_R₀)*P.γ/(P.δ + P.ϵ*(1-P.δ))
results_B2 = KenyaCoV.run_scenario(P,prob,1000,treatment_rates,cb_iso_limit)
@save joinpath(homedir(),"Github/KenyaCoVOutputs/results_B2.jld2") results_B2
results_B2 = 0
