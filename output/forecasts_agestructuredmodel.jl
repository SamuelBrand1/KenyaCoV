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
u0[KenyaCoV.ind_nairobi_as,5,4] = 5#Five initial infecteds in Nairobi in the 20-24 age group
prob = KenyaCoV.create_KenyaCoV_non_neg_prob(u0,(0.,365.),P)
"""
Declare the treatment/isolation rates considered, and a callback that limits intervention
"""
treatment_rates = [0.,1/7.,1/3.5,1]

function isolation_limit(u,t,integrator) # Isolation can only continue until 10000 detected cases
  integrator.p.isolating_detecteds && sum(u[:,:,8]) > 1e4
end
function affect_isolation_limit!(integrator)
  integrator.p.τ = 0
  integrator.p.isolating_detecteds = false
end
cb_iso_limit = DiscreteCallback(isolation_limit,affect_isolation_limit!)




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
P.σ = 1/rand(KenyaCoV.d_incubation)
P.β = rand(KenyaCoV.d_R₀)*P.γ/(P.δ + P.ϵ*(1-P.δ))

results_A = KenyaCoV.run_scenario(P,prob,1000,treatment_rates,cb_iso_limit)

@save joinpath(homedir(),"Github/KenyaCoVOutputs/results_A.jld2") results_A

"""
SCENARIO B:
* No age-dependent susceptibility
* Asymptomatics are not infectious
* 20% of infections are symptomatics
"""
P.χ = ones(KenyaCoV.n_a)
results_B = KenyaCoV.run_scenario(P,prob,1000,treatment_rates,cb_iso_limit)
@save joinpath(homedir(),"Github/KenyaCoVOutputs/results_B.jld2") results_B
"""
SCENARIO C:
* Age-specific susceptibilties calculated from Chinese case data
* Asymptomatics are 25% as infectious as symptomatics
* 20% of infections are symptomatics
"""

u0,P,P_dest = KenyaCoV.model_ingredients_from_data("data/data_for_age_structuredmodel.jld2",
                                            "data/flight_numbers.csv",
                                            "data/projected_global_prevelance.csv")
P.dt = 0.25;
P.ext_inf_rate = 0.;
P.ϵ = 0.25
P.δ = 0.2
P.γ = 1/2.5
P.σ = 1/rand(KenyaCoV.d_incubation)
P.β = rand(KenyaCoV.d_R₀)*P.γ/(P.δ + P.ϵ*(1-P.δ))
results_C = KenyaCoV.run_scenario(P,prob,1000,treatment_rates,cb_iso_limit)
@save joinpath(homedir(),"Github/KenyaCoVOutputs/results_C.jld2") results_C
results_C = 0

"""
SCENARIO D:
* Age-specific susceptibilties calculated from Chinese case data
* Asymptomatics are 100% as infectious as symptomatics
* 20% of infections are symptomatics
"""

u0,P,P_dest = KenyaCoV.model_ingredients_from_data("data/data_for_age_structuredmodel.jld2",
                                            "data/flight_numbers.csv",
                                            "data/projected_global_prevelance.csv")
P.dt = 0.25;
P.ext_inf_rate = 0.;
P.ϵ = 1.
P.δ = 0.2
P.γ = 1/2.5
P.σ = 1/rand(KenyaCoV.d_incubation)
P.β = rand(KenyaCoV.d_R₀)*P.γ/(P.δ + P.ϵ*(1-P.δ))
results_D = KenyaCoV.run_scenario(P,prob,1000,treatment_rates,cb_iso_limit)

@save joinpath(homedir(),"Github/KenyaCoVOutputs/results_D.jld2") results_D
results_D = 0
"""
SCENARIO E:
* Age-specific susceptibilties calculated from Chinese case data
* Asymptomatics are 0% as infectious as symptomatics
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
results_E = KenyaCoV.run_scenario(P,prob,1000,treatment_rates,cb_iso_limit)

@save joinpath(homedir(),"Github/KenyaCoVOutputs/results_E.jld2") results_E
results_E = 0
"""
SCENARIO F:
* Age-specific susceptibilties calculated from Chinese case data
* Asymptomatics are 100% as infectious as symptomatics
* 80% of infections are symptomatics
"""

u0,P,P_dest = KenyaCoV.model_ingredients_from_data("data/data_for_age_structuredmodel.jld2",
                                            "data/flight_numbers.csv",
                                            "data/projected_global_prevelance.csv")
P.dt = 0.25;
P.ext_inf_rate = 0.;
P.ϵ = 1.
P.δ = 0.8
P.γ = 1/2.5
P.σ = 1/rand(KenyaCoV.d_incubation)
P.β = rand(KenyaCoV.d_R₀)*P.γ/(P.δ + P.ϵ*(1-P.δ))
results_F = KenyaCoV.run_scenario(P,prob,1000,treatment_rates,cb_iso_limit)

@save joinpath(homedir(),"Github/KenyaCoVOutputs/results_F.jld2") results_F
results_F = 0
