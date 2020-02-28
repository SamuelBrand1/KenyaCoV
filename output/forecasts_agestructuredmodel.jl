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



"""
SCENARIO A:
* Liu et al estimates
* 14 days intervention on diseased cases
* Age-specific susceptibilties calculated from Chinese case data
* Asymptomatics are not infectious
* 20% of infections are symptomatics
-----------------------
"""
P.dt = 0.25;
P.ext_inf_rate = 0.;
# P.τ = 1/21;
P.ϵ = 0.
P.δ = 0.2
P.γ = 1/2.5
P.σ = 1/rand(d_incubation)
P.β = rand(d_R₀)*P.γ/(P.δ + P.ϵ*(1-P.δ))
u0[KenyaCoV.ind_nairobi_as,5,4] = 5#Five initial infecteds in Nairobi in the 20-24 age group
prob = KenyaCoV.create_KenyaCoV_non_neg_prob(u0,(0.,365.),P)

treatment_rates = [0.]
results_A = run_scenario(P,prob,10,treatment_rates)

sims = run_simulations(P,prob,10,0.)
z = incidence_from_sims(sims)
sims[1][2]





scatter(C_age[:,1],yerror = (C_age[:,2],C_age[:,3]))
boxplot(peaks[:,:]')
plot(1:365,B[12,:,1],ribbon = (B[12,:,2],B[12,:,3]))

"""
SCENARIO A:
* Liu et al estimates
* Age-specific susceptibilties calculated from Chinese case data
* Asymptomatics are not infectious
* 20% of infections are symptomatics
-----------------------

7 days intervention on diseased cases
"""
P.dt = 0.25;
P.ext_inf_rate = 0.;
P.τ = 1/7;
P.ϵ = 0.
P.δ = 0.2
P.γ = 1/2.5
P.σ = 1/5.2
P.β = 2.2*P.γ/(P.δ + P.ϵ*(1-P.δ))


u0[KenyaCoV.ind_nairobi_as,5,4] = 5#Five initial infecteds in Nairobi in the 20-24 age group
prob = KenyaCoV.create_KenyaCoV_non_neg_prob(u0,(0.,365.),P)
ensemble_prob = EnsembleProblem(prob,output_func = output_daily_prev_and_cum_infections)
sims_intervention_7 = solve(ensemble_prob,FunctionMap(),dt = P.dt,trajectories = 1000)

"""
Analysis of scenario A results
"""
area_asymptomatics_scenA_base,area_diseased_scenA_base = cases_from_sims(sims)
area_asymptomatics_scenA_14days,area_diseased_scenA_14days = cases_from_sims(sims_intervention_14)
area_asymptomatics_scenA_7days,area_diseased_scenA_7days = cases_from_sims(sims_intervention_7)
#Saving the output elsewhere because it is a bit large for Github
@save joinpath(homedir(),"Github/KenyaCoVOutputs/cases_scenario_A.jld2") area_asymptomatics_scenA_base area_asymptomatics_scenA_14days area_asymptomatics_scenA_7days area_diseased_scenA_base area_diseased_scenA_14days area_diseased_scenA_7days

total_cases_D = hcat(area_diseased_scenA_base[:,end],area_diseased_scenA_14days[:,end],area_diseased_scenA_7days[:,end])
total_cases_A = hcat(area_asymptomatics_scenA_base[:,end],area_asymptomatics_scenA_14days[:,end],area_asymptomatics_scenA_7days[:,end])

total_cases_collated = hcat(area_diseased_scenA_base[:,end],area_asymptomatics_scenA_base[:,end],
                            area_diseased_scenA_14days[:,end],area_asymptomatics_scenA_14days[:,end],
                            area_diseased_scenA_7days[:,end],area_asymptomatics_scenA_7days[:,end])

plt = boxplot(total_cases_collated,lab = ["Symptomatic" "Asymptomatic" "" "" "" ""],
                    ylabel = "Total cases",
                    xticks = ([1.5,3.5,5.5],["No intervention" "Av. 14 days to intervention" "Av. 7 days to intervention"]),
                    color = [:red :blue],
                    title = "Scenario A: Asymps uninfectious, 20% symptomatic rate ")
savefig(plt,"output/ScenarioA_cases.pdf")


"""
SCENARIO B:
* Liu et al estimates
* 14 days intervention on diseased cases
* Age-specific susceptibilties calculated from Chinese case data
* Asymptomatics are 10% as infectious as diseased
* 20% of infections are symptomatics
"""
P.dt = 0.25;
P.ext_inf_rate = 0.;
P.τ = 1/14;
P.ϵ = 0.1
P.δ = 0.2
P.γ = 1/2.5
P.σ = 1/5.2
P.β = 2.2*P.γ/(P.δ + P.ϵ*(1-P.δ))


u0[KenyaCoV.ind_nairobi_as,5,4] = 5#Five initial infecteds in Nairobi in the 20-24 age group
prob = KenyaCoV.create_KenyaCoV_non_neg_prob(u0,(0.,365.),P)
ensemble_prob = EnsembleProblem(prob,output_func = output_daily_prev_and_cum_infections)
sims_intervention_scenario4 = solve(ensemble_prob,FunctionMap(),dt = P.dt,trajectories = 1000)
sims_intervention_scenarioB_14days = sims_intervention_scenario4
@save joinpath(homedir(),"Github/KenyaCoVOutputs/sims_sofar.jld2") sims sims_intervention_14 sims_intervention_7 sims_intervention_scenarioB_14days
@load joinpath(homedir(),"Github/KenyaCoVOutputs/sims_sofar.jld2") sims sims_intervention_14 sims_intervention_7 sims_intervention_scenarioB_14days
