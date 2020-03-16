push!(LOAD_PATH, joinpath(homedir(),"GitHub/KenyaCoV/src"))
using Plots,Parameters,Distributions,DifferentialEquations,JLD2,DataFrames,StatsPlots,FileIO,MAT,RecursiveArrayTools
using Revise
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

"""
Declare the treatment/isolation rates considered, and a callback that limits intervention
"""
treatment_rates = [(0.,1),(0.,0.5),(1/7.,1.),(1/7,0.5),(1/3.5,1.),(1/3.5,0.5)]

function isolation_limit(u,t,integrator) # Isolation can only continue until 1000 detected cases
  integrator.p.isolating_detecteds && sum(u[:,:,8]) > 1e3
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
P.β = 2.5*P.γ
#Choosing initial conditions - scenario 1
rand_chosen_R,eq_case_distribution,fifth_gen_case_distribution = KenyaCoV.calculate_R₀(P)
num_initial_detected = ceil(Int64,sum(fifth_gen_case_distribution.*d_0))
detected_age_group = argmax(fifth_gen_case_distribution.*d_0)
undetected_age_profile = round.(Int64,fifth_gen_case_distribution.*(1 .- d_0))

u0[KenyaCoV.ind_nairobi_as,detected_age_group,4] = num_initial_detected
u0[KenyaCoV.ind_nairobi_as,:,3] = undetected_age_profile
prob = KenyaCoV.create_KenyaCoV_non_neg_prob(u0,(0.,365.),P)

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
P.β = 2.5*P.γ

#Choosing initial conditions - scenario 2
rand_chosen_R,eq_case_distribution,fifth_gen_case_distribution = KenyaCoV.calculate_R₀(P)
num_initial_detected = ceil(Int64,sum(fifth_gen_case_distribution.*d_0))
detected_age_group = argmax(fifth_gen_case_distribution.*d_0)
undetected_age_profile = round.(Int64,fifth_gen_case_distribution.*(1 .- d_0))

u0[KenyaCoV.ind_nairobi_as,detected_age_group,4] = num_initial_detected
u0[KenyaCoV.ind_nairobi_as,:,3] = undetected_age_profile
prob = KenyaCoV.create_KenyaCoV_non_neg_prob(u0,(0.,365.),P)

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
P.β = 2.5*P.γ

#Choosing initial conditions - scenario 3
rand_chosen_R,eq_case_distribution,fifth_gen_case_distribution = KenyaCoV.calculate_R₀(P)
num_initial_detected = ceil(Int64,sum(fifth_gen_case_distribution.*d_0))
detected_age_group = argmax(fifth_gen_case_distribution.*d_0)
undetected_age_profile = round.(Int64,fifth_gen_case_distribution.*(1 .- d_0))

u0[KenyaCoV.ind_nairobi_as,detected_age_group,4] = num_initial_detected
u0[KenyaCoV.ind_nairobi_as,:,3] = undetected_age_profile
prob = KenyaCoV.create_KenyaCoV_non_neg_prob(u0,(0.,365.),P)

results_3 = KenyaCoV.run_scenario(P,prob,1000,treatment_rates,cb_iso_limit)
@save joinpath(homedir(),"Github/KenyaCoVOutputs/results_3.jld2") results_3


println("Finished 3")


println("Finished 2")
"""
SCENARIO 4:
* No age-specific susceptibilties --- its disease difference
* Asymptomatics are 50% as infectious as symptomatics
"""
u0,P,P_dest = KenyaCoV.model_ingredients_from_data("data/data_for_age_structuredmodel.jld2",
                                            "data/flight_numbers.csv",
                                            "data/projected_global_prevelance.csv")
#Redistribute susceptibility JUST to rescale infectiousness so we get the correct R₀/r
P.χ = ones(KenyaCoV.n_a)
P.rel_detection_rate = d_05
P.dt = 0.25;
P.ext_inf_rate = 0.;
P.ϵ = 0.5
P.γ = 1/2.5
R₀_scale = KenyaCoV.calculate_R₀_scale(P)
P.χ = ones(KenyaCoV.n_a)/R₀_scale


P.σ = 1/rand(KenyaCoV.d_incubation)
P.β = 2.5*P.γ

#Choosing initial conditions - scenario 3
rand_chosen_R,eq_case_distribution,fifth_gen_case_distribution = KenyaCoV.calculate_R₀(P)
num_initial_detected = ceil(Int64,sum(fifth_gen_case_distribution.*d_05))
detected_age_group = argmax(fifth_gen_case_distribution.*d_05)
undetected_age_profile = round.(Int64,fifth_gen_case_distribution.*(1 .- d_05))

u0[KenyaCoV.ind_nairobi_as,detected_age_group,4] = num_initial_detected
u0[KenyaCoV.ind_nairobi_as,:,3] = undetected_age_profile
prob = KenyaCoV.create_KenyaCoV_non_neg_prob(u0,(0.,365.),P)

results_4 = KenyaCoV.run_scenario(P,prob,1000,treatment_rates,cb_iso_limit)
@save joinpath(homedir(),"Github/KenyaCoVOutputs/results_4.jld2") results_4


println("Finished 4")

"""
SCENARIO 5:
* No age-specific susceptibilties --- its disease difference
* Asymptomatics are 100% as infectious as symptomatics
"""
u0,P,P_dest = KenyaCoV.model_ingredients_from_data("data/data_for_age_structuredmodel.jld2",
                                            "data/flight_numbers.csv",
                                            "data/projected_global_prevelance.csv")
#Redistribute susceptibility JUST to rescale infectiousness so we get the correct R₀/r
P.χ = ones(KenyaCoV.n_a)
P.rel_detection_rate = d_1
P.dt = 0.25;
P.ext_inf_rate = 0.;
P.ϵ = 1.
P.γ = 1/2.5
R₀_scale = KenyaCoV.calculate_R₀_scale(P)
P.χ = ones(KenyaCoV.n_a)/R₀_scale


P.σ = 1/rand(KenyaCoV.d_incubation)
P.β = 2.5*P.γ

#Choosing initial conditions - scenario 3
rand_chosen_R,eq_case_distribution,fifth_gen_case_distribution = KenyaCoV.calculate_R₀(P)
num_initial_detected = ceil(Int64,sum(fifth_gen_case_distribution.*d_1))
detected_age_group = argmax(fifth_gen_case_distribution.*d_1)
undetected_age_profile = round.(Int64,fifth_gen_case_distribution.*(1 .- d_1))

u0[KenyaCoV.ind_nairobi_as,detected_age_group,4] = num_initial_detected
u0[KenyaCoV.ind_nairobi_as,:,3] = undetected_age_profile
prob = KenyaCoV.create_KenyaCoV_non_neg_prob(u0,(0.,365.),P)

results_5 = KenyaCoV.run_scenario(P,prob,1000,treatment_rates,cb_iso_limit)
@save joinpath(homedir(),"Github/KenyaCoVOutputs/results_5.jld2") results_5


println("Finished 5")


"""
SOCIAL DISTANCING SCENARIOS
Load all mixing matrices
"""

@load "data/agemixingmatrix_Kenya_all_types.jld2" M_Kenya M_Kenya_ho M_Kenya_other M_Kenya_school M_Kenya_work
reduced_treatment_rates = [(0.,1),(1/3.5,0.5)]

"""
SCENARIO 2 --- SOCIAL DISTANCING:
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
P.β = 2.5*P.γ


#Choosing initial conditions - scenario 2
rand_chosen_R,eq_case_distribution,fifth_gen_case_distribution = KenyaCoV.calculate_R₀(P)
num_initial_detected = ceil(Int64,sum(fifth_gen_case_distribution.*d_01))
detected_age_group = argmax(fifth_gen_case_distribution.*d_01)
undetected_age_profile = round.(Int64,fifth_gen_case_distribution.*(1 .- d_01))


#Implement social distancing
KenyaCoV.calculate_R₀(P)
P.M = M_Kenya_ho .+ M_Kenya_work
KenyaCoV.calculate_R₀(P)


u0[KenyaCoV.ind_nairobi_as,detected_age_group,4] = num_initial_detected
u0[KenyaCoV.ind_nairobi_as,:,3] = undetected_age_profile
prob = KenyaCoV.create_KenyaCoV_non_neg_prob(u0,(0.,365.),P)

results_2S = KenyaCoV.run_scenario(P,prob,1000,reduced_treatment_rates,cb_iso_limit)
@save joinpath(homedir(),"Github/KenyaCoVOutputs/results_2S.jld2") results_2S


println("Finished 2 social distancing")
"""
SCENARIO 3 --- Social distancing:
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
P.β = 2.5*P.γ




#Choosing initial conditions - scenario 3
rand_chosen_R,eq_case_distribution,fifth_gen_case_distribution = KenyaCoV.calculate_R₀(P)
num_initial_detected = ceil(Int64,sum(fifth_gen_case_distribution.*d_025))
detected_age_group = argmax(fifth_gen_case_distribution.*d_025)
undetected_age_profile = round.(Int64,fifth_gen_case_distribution.*(1 .- d_025))

#Implement social distancing
KenyaCoV.calculate_R₀(P)
P.M = M_Kenya_ho .+ M_Kenya_work
KenyaCoV.calculate_R₀(P)


u0[KenyaCoV.ind_nairobi_as,detected_age_group,4] = num_initial_detected
u0[KenyaCoV.ind_nairobi_as,:,3] = undetected_age_profile
prob = KenyaCoV.create_KenyaCoV_non_neg_prob(u0,(0.,365.),P)

results_3S = KenyaCoV.run_scenario(P,prob,1000,reduced_treatment_rates,cb_iso_limit)
@save joinpath(homedir(),"Github/KenyaCoVOutputs/results_3S.jld2") results_3S


println("Finished 3 Social distancing")
