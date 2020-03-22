push!(LOAD_PATH, joinpath(homedir(),"GitHub/KenyaCoV/src"))
using Plots,Parameters,Distributions,DifferentialEquations,JLD2,DataFrames,StatsPlots,FileIO,MAT,RecursiveArrayTools
using Revise
import KenyaCoV
using LinearAlgebra:eigen
using Statistics: median, quantile

gr()#Plotting frontend
"""
IMPERFECT social distancing + spatial lockdown + decrease in effective contacts by 50% over 30 days
All households reduce contact outside household or workplace by 75%.
School contact rates set to zero, workplace contact rates reduced by 25%.
Household contact rates assumed to increase by 25%.
All spatial movements reduced to 0.1% of time spent in other areas
"""

"""
Load age structured data, define initial state and declare the KenyaCoV problem
"""

u0,P,P_dest = KenyaCoV.model_ingredients_from_data("data/data_for_age_structuredmodel.jld2",
                                            "data/flight_numbers.csv",
                                            "data/projected_global_prevelance.csv")
T_normal = deepcopy(P.T)

@load "data/detection_rates_for_different_taus.jld2" d_0 d_01 d_025 d_05 d_1
@load "data/susceptibility_rates.jld2" σ

#Initial infecteds
"""
Load age mixing matrices (these are all in to (row) from (col) format)
"""

@load "data/agemixingmatrix_Kenya_all_types.jld2" M_Kenya M_Kenya_ho M_Kenya_other M_Kenya_school M_Kenya_work



"""
Declare the treatment/isolation rates considered, and a callback that limits intervention
"""
treatment_rates = [(0.,1),(0.,0.5),(1/7.,1.),(1/7,0.5),(1/3.5,1.),(1/3.5,0.5)]
reduced_treatment_rates = [(1/3.5,0.5)]
function ramp_down(t)
    if t < 60.
        return (1-t/60) + 0.5*t/60
    else
        return 0.5
    end
end

function isolation_limit(u,t,integrator) # Isolation can only continue until 1000 detected cases
  integrator.p.isolating_detecteds && sum(u[:,:,8]) > 1e3
end
function affect_isolation_limit!(integrator)
  integrator.p.τ = 0
  integrator.p.isolating_detecteds = false
end
cb_iso_limit = DiscreteCallback(isolation_limit,affect_isolation_limit!)

function social_distancing_limit(u,t,integrator) # Isolation can only continue until 1000 detected cases
  integrator.p.lockdown && t > 90.
end
function affect_social_distancing_limit!(integrator)
  integrator.p.M = M_Kenya
  integrator.p.T = T_normal
  integrator.p.lockdown = false
  integrator.p.c_t = t -> 1.
end

cb_SD_limit = DiscreteCallback(social_distancing_limit,affect_social_distancing_limit!)

both_cbs = CallbackSet(cb_iso_limit,cb_SD_limit)


"""
SCENARIO 1 --- Imperfect social distancing + spatial lockdown for 3 months:
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


#Implement social distancing + spatial lockdown
KenyaCoV.calculate_R₀(P)
P.M = M_Kenya_ho*(1.25) .+ M_Kenya_work*(0.75) .+ M_Kenya_other*(0.25)
KenyaCoV.calculate_R₀(P)
KenyaCoV.transportstructure_params!(P,[0.001 for i = 1:KenyaCoV.n_wa],P_dest)
P.c_t = ramp_down

u0[KenyaCoV.ind_nairobi_as,detected_age_group,4] = num_initial_detected
u0[KenyaCoV.ind_nairobi_as,:,3] = undetected_age_profile
prob = KenyaCoV.create_KenyaCoV_non_neg_prob(u0,(0.,2*365.),P)

results_1SI_SL_DC_three_months = KenyaCoV.run_scenario(P,prob,1000,reduced_treatment_rates,both_cbs)
@save joinpath(homedir(),"Github/KenyaCoVOutputs/results_1SI_SL_DC_three_months.jld2") results_1SI_SL_DC_three_months


println("Finished 1 imperfect social distancing + spatial lockdown for 3 months")



"""
SCENARIO 2 --- Imperfect social distancing + spatial lockdown:
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


#Implement social distancing + spatial lockdown
KenyaCoV.calculate_R₀(P)
P.M = M_Kenya_ho*(1.25) .+ M_Kenya_work*(0.75) .+ M_Kenya_other*(0.25)
KenyaCoV.calculate_R₀(P)
KenyaCoV.transportstructure_params!(P,[0.001 for i = 1:KenyaCoV.n_wa],P_dest)
P.c_t = ramp_down


u0[KenyaCoV.ind_nairobi_as,detected_age_group,4] = num_initial_detected
u0[KenyaCoV.ind_nairobi_as,:,3] = undetected_age_profile
prob = KenyaCoV.create_KenyaCoV_non_neg_prob(u0,(0.,2*365.),P)

results_2SI_SL_DC_three_months = KenyaCoV.run_scenario(P,prob,1000,reduced_treatment_rates,both_cbs)
@save joinpath(homedir(),"Github/KenyaCoVOutputs/results_2SI_SL_DC_three_months.jld2") results_2SI_SL_DC_three_months


println("Finished 2 imperfect social distancing + spatial lockdown for 3 months")
"""
SCENARIO 3 --- social distancing + spatial lockdown:
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

#Implement social distancing + spatial lockdown
KenyaCoV.calculate_R₀(P)
P.M = M_Kenya_ho*(1.25) .+ M_Kenya_work*(0.75) .+ M_Kenya_other*(0.25)
KenyaCoV.calculate_R₀(P)
KenyaCoV.transportstructure_params!(P,[0.001 for i = 1:KenyaCoV.n_wa],P_dest)
P.c_t = ramp_down

u0[KenyaCoV.ind_nairobi_as,detected_age_group,4] = num_initial_detected
u0[KenyaCoV.ind_nairobi_as,:,3] = undetected_age_profile
prob = KenyaCoV.create_KenyaCoV_non_neg_prob(u0,(0.,2*365.),P)

results_3SI_SL_DC_three_months = KenyaCoV.run_scenario(P,prob,1000,reduced_treatment_rates,both_cbs)
@save joinpath(homedir(),"Github/KenyaCoVOutputs/results_3SI_SL_DC_three_months.jld2") results_3SI_SL_DC_three_months


println("Finished 3 imperfect social distancing + spatial lockdown for 3 months")


"""
SCENARIO 4 --- social distancing + spatial lockdown for 3 months:
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




#Choosing initial conditions - scenario 4
rand_chosen_R,eq_case_distribution,fifth_gen_case_distribution = KenyaCoV.calculate_R₀(P)
num_initial_detected = ceil(Int64,sum(fifth_gen_case_distribution.*d_05))
detected_age_group = argmax(fifth_gen_case_distribution.*d_05)
undetected_age_profile = round.(Int64,fifth_gen_case_distribution.*(1 .- d_05))

#Implement social distancing + spatial lockdown
KenyaCoV.calculate_R₀(P)
P.M = M_Kenya_ho*(1.25) .+ M_Kenya_work*(0.75) .+ M_Kenya_other*(0.25)
KenyaCoV.calculate_R₀(P)
KenyaCoV.transportstructure_params!(P,[0.001 for i = 1:KenyaCoV.n_wa],P_dest)
P.c_t = ramp_down


u0[KenyaCoV.ind_nairobi_as,detected_age_group,4] = num_initial_detected
u0[KenyaCoV.ind_nairobi_as,:,3] = undetected_age_profile
prob = KenyaCoV.create_KenyaCoV_non_neg_prob(u0,(0.,2*365.),P)

results_4SI_SL_DC_three_months = KenyaCoV.run_scenario(P,prob,1000,reduced_treatment_rates,both_cbs)
@save joinpath(homedir(),"Github/KenyaCoVOutputs/results_4SI_SL_DC_three_months.jld2") results_4SI_SL_DC_three_months


println("Finished 4 imperfect social distancing + spatial lockdown for 3 months")


"""
SCENARIO 5 --- social distancing + spatial lockdown for 3 months:
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




#Choosing initial conditions - scenario 5
rand_chosen_R,eq_case_distribution,fifth_gen_case_distribution = KenyaCoV.calculate_R₀(P)
num_initial_detected = ceil(Int64,sum(fifth_gen_case_distribution.*d_1))
detected_age_group = argmax(fifth_gen_case_distribution.*d_1)
undetected_age_profile = round.(Int64,fifth_gen_case_distribution.*(1 .- d_1))

#Implement social distancing + spatial lockdown
KenyaCoV.calculate_R₀(P)
P.M = M_Kenya_ho*(1.25) .+ M_Kenya_work*(0.75) .+ M_Kenya_other*(0.25)
KenyaCoV.calculate_R₀(P)
KenyaCoV.transportstructure_params!(P,[0.001 for i = 1:KenyaCoV.n_wa],P_dest)
P.c_t = ramp_down


u0[KenyaCoV.ind_nairobi_as,detected_age_group,4] = num_initial_detected
u0[KenyaCoV.ind_nairobi_as,:,3] = undetected_age_profile
prob = KenyaCoV.create_KenyaCoV_non_neg_prob(u0,(0.,2*365.),P)

results_5SI_SL_DC_three_months = KenyaCoV.run_scenario(P,prob,1000,reduced_treatment_rates,both_cbs)
@save joinpath(homedir(),"Github/KenyaCoVOutputs/results_5SI_SL_DC_three_months.jld2") results_5SI_SL_DC_three_months


println("Finished 5 imperfect social distancing + spatial lockdown for 3 months")
