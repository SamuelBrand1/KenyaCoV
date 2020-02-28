push!(LOAD_PATH, joinpath(homedir(),"GitHub/KenyaCoV/src"))
using Plots,Parameters,Distributions,DifferentialEquations,JLD2,DataFrames,StatsPlots,FileIO,MAT,RecursiveArrayTools
# using Revise
import KenyaCoV
using LinearAlgebra:eigen
using Statistics: median, quantile
gr()
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
Output functions: Target the peak timing for each county, cases by each county,
timing of total peak and total cases
"""

function output_daily_and_final_incidence(sol,i)
    times = sol.prob.tspan[1]:1:sol.prob.tspan[end]
    z = [sum(sol(t)[:,:,7:8],dims = 2)[:,1,:]  for t in times]
    return (z,sol[end][:,:,7:8]),false
end

function randomise_params(prob,i,repeat)
    _P = deepcopy(prob.p)
    _P.σ = 1/rand(d_incubation)
    _P.β = rand(d_R₀)*_P.γ/(_P.δ + _P.ϵ*(1-_P.δ))
    return remake(prob,p=_P)
end

"""
Simulation functions
"""
function run_simulations(P::KenyaCoV.CoVParameters_AS,prob,n_traj,τ)
    P.τ = τ
    ensemble_prob = EnsembleProblem(prob,
                                    prob_func = randomise_params,
                                    output_func = output_daily_and_final_incidence)
    return solve(ensemble_prob,FunctionMap(),dt = P.dt,trajectories = n_traj)
end

function run_scenario(P::KenyaCoV.CoVParameters_AS,prob,n_traj,treatment_rates)
    results = []
    for τ in treatment_rates
        sims = run_simulations(P,prob,n_traj,τ)
        analysisdata = incidence_from_sims(sims)
        push!(results,analysisdata)
    end
    return results
end

"""
Analysis functions:
"""
function incidence_from_sims(sims)
    n = length(sims)
    m = length(sims[1][1])
    inc_A_data = zeros(21,m-1,n)
    inc_D_data = zeros(21,m-1,n)
    final_case_data = zeros(21,16,n)
    for (k,sim) in enumerate(sims)
        for i = 1:20,t = 1:(m-1)
            inc_A_data[i,t,k] = sim[1][t+1][i,1] - sim[1][t][i,1]#Daily incidence rather than cumulative
            inc_D_data[i,t,k] = sim[1][t+1][i,2] - sim[1][t][i,2]
        end
        for t = 1:(m-1)
            inc_A_data[21,t,k] = sum(sim[1][t+1][:,1] .- sim[1][t][:,1])
            inc_D_data[21,t,k] = sum(sim[1][t+1][:,2] .- sim[1][t][:,2])
        end
        for i =1:20,a=1:16
            final_case_data[i,a,k] = sum(sim[2][i,a,:])
        end
        for a = 1:16
            final_case_data[21,a,k] = sum(sim[2][:,a,:])
        end
    end
    peak_times = zeros(n,21)
    inc_A_conf_intvs = zeros(21,m-1,3)
    inc_D_conf_intvs = zeros(21,m-1,3)
    final_case_intvs = zeros(21,16,3)
    for i = 1:21,t = 1:(m-1)
        inc_A_conf_intvs[i,t,1] = quantile(inc_A_data[i,t,:],0.5)
        inc_A_conf_intvs[i,t,2] = inc_A_conf_intvs[i,t,1] - quantile(inc_A_data[i,t,:],0.025) #Lower conf. int.
        inc_A_conf_intvs[i,t,3] =  quantile(inc_A_data[i,t,:],0.975) - inc_A_conf_intvs[i,t,1] #Upper conf. int.
        inc_D_conf_intvs[i,t,1] = quantile(inc_D_data[i,t,:],0.5)
        inc_D_conf_intvs[i,t,2] = inc_D_conf_intvs[i,t,1] - quantile(inc_D_data[i,t,:],0.025) #Lower conf. int.
        inc_D_conf_intvs[i,t,3] =  quantile(inc_D_data[i,t,:],0.975) - inc_D_conf_intvs[i,t,1] #Upper conf. int.
    end

    for i = 1:21, k =1:n
        peak_times[k,i] = argmax(inc_A_data[i,:,k])
    end

    for i = 1:21,a=1:16
        final_case_intvs[i,a,1] = quantile(final_case_data[i,a,:],0.5)
        final_case_intvs[i,a,2] = quantile(final_case_data[i,a,:],0.5) - quantile(final_case_data[i,a,:],0.025)
        final_case_intvs[i,a,3] = quantile(final_case_data[i,a,:],0.975) - quantile(final_case_data[i,a,:],0.5)
    end



    return inc_A_conf_intvs,inc_D_conf_intvs,peak_times,final_case_intvs
end




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

treatment_rates = [0.,1/21.,1/14,1/7]
results_A = run_scenario(P,prob,1000,treatment_rates)

# @save joinpath(homedir(),"Github/KenyaCoVOutputs/results_A.jld2") results_A


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
