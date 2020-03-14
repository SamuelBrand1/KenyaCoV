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


@load "data/detection_rates_for_different_taus.jld2" d_0 d_01 d_025 d_05 d_1

#Initial infecteds
u0[KenyaCoV.ind_nairobi_as,5,3] = 20
u0[KenyaCoV.ind_nairobi_as,5,4] = 5
prob = KenyaCoV.create_KenyaCoV_non_neg_prob(u0,(0.,365.),P)
"""
Declare the treatment/isolation rates considered, and a callback that limits intervention
"""
treatment_rates = [(0.,1),(0.,0.5),(1/7.,1.),(1/7,0.5),(1/3.5,1.),(1/3.5,0.5)]

function isolation_limit(u,t,integrator) # Isolation can only continue until 5000 detected cases
  integrator.p.isolating_detecteds && sum(u[:,:,8]) > 5e3
end
function affect_isolation_limit!(integrator)
  integrator.p.τ = 0
  integrator.p.isolating_detecteds = false
end
cb_iso_limit = DiscreteCallback(isolation_limit,affect_isolation_limit!)



"""
SCENARIO 1:
* No age-specific susceptibilties --- its disease difference
* Asymptomatics are 0% as infectious as symptomatics
"""
u0,P,P_dest = KenyaCoV.model_ingredients_from_data("data/data_for_age_structuredmodel.jld2",
                                            "data/flight_numbers.csv",
                                            "data/projected_global_prevelance.csv")
#Redistribute susceptibility JUST to rescale infectiousness so we get the correct R₀/r
P.χ = ones(KenyaCoV.n_a)
P.rel_detection_rate = d_0
P.dt = 0.25;
P.ext_inf_rate = 0.;
P.ϵ = 0.
P.γ = 1/2.5
R₀_scale = KenyaCoV.calculate_R₀_scale(P)
P.χ = ones(KenyaCoV.n_a)/R₀_scale


P.σ = 1/rand(KenyaCoV.d_incubation)
P.β = rand(KenyaCoV.d_R₀)*P.γ
results_1 = KenyaCoV.run_scenario(P,prob,1000,treatment_rates,cb_iso_limit)
@save joinpath(homedir(),"Github/KenyaCoVOutputs/results_1.jld2") results_1


println("Finished 1")



"""
SCENARIO 2:
* No age-specific susceptibilties --- its disease difference
* Asymptomatics are 10% as infectious as symptomatics
"""
u0,P,P_dest = KenyaCoV.model_ingredients_from_data("data/data_for_age_structuredmodel.jld2",
                                            "data/flight_numbers.csv",
                                            "data/projected_global_prevelance.csv")
#Redistribute susceptibility JUST to rescale infectiousness so we get the correct R₀/r
P.χ = ones(KenyaCoV.n_a)
P.rel_detection_rate = d_01
P.dt = 0.25;
P.ext_inf_rate = 0.;
P.ϵ = 0.1
P.γ = 1/2.5
R₀_scale = KenyaCoV.calculate_R₀_scale(P)
P.χ = ones(KenyaCoV.n_a)/R₀_scale


P.σ = 1/rand(KenyaCoV.d_incubation)
P.β = rand(KenyaCoV.d_R₀)*P.γ
results_2 = KenyaCoV.run_scenario(P,prob,1000,treatment_rates,cb_iso_limit)
@save joinpath(homedir(),"Github/KenyaCoVOutputs/results_2.jld2") results_2


println("Finished 2")
"""
SCENARIO 3:
* No age-specific susceptibilties --- its disease difference
* Asymptomatics are 25% as infectious as symptomatics
"""
u0,P,P_dest = KenyaCoV.model_ingredients_from_data("data/data_for_age_structuredmodel.jld2",
                                            "data/flight_numbers.csv",
                                            "data/projected_global_prevelance.csv")
#Redistribute susceptibility JUST to rescale infectiousness so we get the correct R₀/r
P.χ = ones(KenyaCoV.n_a)
P.rel_detection_rate = d_025
P.dt = 0.25;
P.ext_inf_rate = 0.;
P.ϵ = 0.25
P.γ = 1/2.5
R₀_scale = KenyaCoV.calculate_R₀_scale(P)
P.χ = ones(KenyaCoV.n_a)/R₀_scale


P.σ = 1/rand(KenyaCoV.d_incubation)
P.β = rand(KenyaCoV.d_R₀)*P.γ
results_3 = KenyaCoV.run_scenario(P,prob,1000,treatment_rates,cb_iso_limit)
@save joinpath(homedir(),"Github/KenyaCoVOutputs/results_3.jld2") results_3


println("Finished 3")
