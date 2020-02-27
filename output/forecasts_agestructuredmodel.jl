push!(LOAD_PATH, joinpath(homedir(),"GitHub/KenyaCoV/src"))
using Plots,Parameters,Distributions,DifferentialEquations,JLD2,DataFrames,StatsPlots,FileIO,MAT,RecursiveArrayTools
# using Revise
import KenyaCoV
using LinearAlgebra:eigen
using Statistics: quantile

gr()
"""
Load age structured data
"""

u0,P,P_dest = KenyaCoV.model_ingredients_from_data("data/data_for_age_structuredmodel.jld2",
                                            "data/flight_numbers.csv",
                                            "data/projected_global_prevelance.csv")

"""
Output functions: Target the peak timing for each county, cases by each county,
timing of total peak and total cases
"""

function output_daily_prev_and_cum_infections(sol,i)
    times = sol.prob.tspan[1]:1:sol.prob.tspan[end]
    z = [sol(t)[:,:,[3,4,7,8]] for t in times]
    return z,false
end

"""
Analysis functions:
"""
function cases_from_sims(sims)
    total_asymptomatics = [sum(sim[end][:,:,3]) for sim in sims]
    total_diseased = [sum(sim[end][:,:,4]) for sim in sims]
    area_asymptomatics = zeros(1000,KenyaCoV.n_wa)
    area_diseased = zeros(1000,KenyaCoV.n_wa)
    for (i,sim) in enumerate(sims)
        for j in 1:KenyaCoV.n_wa
            area_asymptomatics[i,j] = sum(sim[end][j,:,3])
            area_diseased[i,j] = sum(sim[end][j,:,4])
        end
    end
    area_asymptomatics = hcat(area_asymptomatics,total_asymptomatics)
    area_diseased = hcat(area_diseased,total_diseased)
    return area_asymptomatics,area_diseased
end


# plot(t,y,ribbon = (σ₋,σ₊),fillalpha=0.3,color=:red,lw = 3)




"""
SCENARIO A:
* Liu et al estimates
* Age-specific susceptibilties calculated from Chinese case data
* Asymptomatics are not infectious
* 20% of infections are symptomatics
------------------------------------
NO INTERVENTION
"""
P.dt = 0.25;
P.ext_inf_rate = 0.;
P.τ = 0.;
P.ϵ = 0.
P.δ = 0.2
P.γ = 1/2.5
P.σ = 1/5.2
P.β = 2.2*P.γ/(P.δ + P.ϵ*(1-P.δ))
sus_matrix = repeat(P.χ,1,KenyaCoV.n_a)
eigs, = eigen((P.δ + P.ϵ*(1-P.δ))*(P.β/P.γ)*sus_matrix.*P.M)
Kenya_R₀ = Real(eigs[end])#R₀ for Kenya

u0[KenyaCoV.ind_nairobi_as,5,4] = 5#Five initial infecteds in Nairobi in the 20-24 age group
prob = KenyaCoV.create_KenyaCoV_non_neg_prob(u0,(0.,365.),P)
ensemble_prob = EnsembleProblem(prob,output_func = output_daily_prev_and_cum_infections)
sims = solve(ensemble_prob,FunctionMap(),dt = P.dt,trajectories = 1000)



"""
SCENARIO A:
* Liu et al estimates
* 14 days intervention on diseased cases
* Age-specific susceptibilties calculated from Chinese case data
* Asymptomatics are not infectious
* 20% of infections are symptomatics
-----------------------
14 days intervention on diseased cases
"""
P.dt = 0.25;
P.ext_inf_rate = 0.;
P.τ = 1/14;
P.ϵ = 0.
P.δ = 0.2
P.γ = 1/2.5
P.σ = 1/5.2
P.β = 2.2*P.γ/(P.δ + P.ϵ*(1-P.δ))


u0[KenyaCoV.ind_nairobi_as,5,4] = 5#Five initial infecteds in Nairobi in the 20-24 age group
prob = KenyaCoV.create_KenyaCoV_non_neg_prob(u0,(0.,365.),P)
ensemble_prob = EnsembleProblem(prob,output_func = output_daily_prev_and_cum_infections)
sims_intervention_14 = solve(ensemble_prob,FunctionMap(),dt = P.dt,trajectories = 1000)


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
