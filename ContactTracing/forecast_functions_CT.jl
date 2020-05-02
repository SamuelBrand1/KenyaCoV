
using Plots,Parameters,Distributions,DifferentialEquations,JLD2,DataFrames,StatsPlots,FileIO,MAT,RecursiveArrayTools
import KenyaCoV_CT
using LinearAlgebra:eigen
using Statistics: median, quantile

###########
function randomise_params(prob,i,repeat)
    _P = deepcopy(prob.p)
    _P.σ = 1/rand(d_incubation)
    _P.β = rand(d_R₀)*_P.γ/(_P.δ + _P.ϵ*(1-_P.δ))
    return remake(prob,p=_P)
end
function output_daily_and_final_incidence(sol,i)
    times = sol.prob.tspan[1]:1:sol.prob.tspan[end] #time step = 1 day
    cumIs = [sum(sol(t)[:,:,9:11],dims = 2)[:,1,:]  for t in times] # only cumulative A, M and V
    I=[[sum(sol(t)[wa,:,4:6],dims = 1)[:,1,:] #=+ sum(sol(t)[wa,:,5],dims = 1)[:,1,:] + sum(sol(t)[wa,:,6],dims = 1)[:,1,:]=#  for wa=1:20]  for t in times]      # I = A+M+V
    return [cumIs,I,sol[end][:,:,9:11]],false # save z (time series with only incidence with no age structure), and save the final distribution (end) age/space but no time
end

function run_set_scenarios(folder,session,scenarios,β,ϵ,n_traj)
    if !isdir(folder)   mkdir(folder)   end
    global sims_vector=[]
    for i=1:size(scenarios,1)
        sc_nb=session*100+scenarios[i]
        u0,P,P_dest = KenyaCoV_CT.model_ingredients_from_data("data/data_for_age_structuredmodel.jld2","data/flight_numbers.csv","data/projected_global_prevelance.csv")
        #Redistribute susceptibility JUST to rescale infectiousness so we get the correct R₀/r
        P.χ = ones(KenyaCoV_CT.n_a)
        @load "data/detection_rates_for_different_taus.jld2" d_1
        P.rel_detection_rate = d_1
        P.dt = 0.25;
        P.ext_inf_rate = 0.;
        P.ϵ = ϵ#1.#Asymp are as infectious as symptomatics
        #Set the susceptibility vector --- just to specify the correct R₀
        sus_matrix = repeat(P.χ,1,17)
        R_A = P.ϵ*((1/P.σ₂) + (1/P.γ) ) #effective duration of asymptomatic
        R_M = (P.ϵ/P.σ₂) + (P.ϵ_D/P.γ) #effective duration of mild
        R_V = (P.ϵ/P.σ₂) + (P.ϵ_V/P.τ) #effective duration of severe
        R_vector = [(1-P.rel_detection_rate[a])*R_A + P.rel_detection_rate[a]*(1-P.hₐ[a])*R_M + P.rel_detection_rate[a]*P.hₐ[a]*R_V for a = 1:17]
        inf_matrix = repeat(R_vector',17,1)

        eigs, = eigen(sus_matrix.*P.M.*inf_matrix)
        max_eigval = Real(eigs[end])
        P.χ = ones(KenyaCoV_CT.n_a)/max_eigval
        P.β = β #rand(KenyaCoV.d_R₀) #Choose R₀ randomly from 2-3 range
        u0[KenyaCoV_CT.ind_nairobi_as,8,3] = 5      #Initial infecteds (P) in Nairobi in the 30-34 age group
        #P.MS_strategy=MS_strategies[i]
        prob = KenyaCoV_CT.create_KenyaCoV_non_neg_prob(u0,(0.,365.),P)
        print("Simulating session ",session," sc ",sc_nb,"   -   ")
        #@time sims = solve(EnsembleProblem(prob#=,prob_func=randomise_params=#,output_func = output_daily_and_final_incidence),
        #                    FunctionMap(),dt=P.dt,trajectories=n_traj)
        #@time sims = [init(prob,output_func = output_daily_and_final_incidence,Tsit5(),dt=P.dt) for traj=1:n_traj]
        #@time sims = init(prob,FunctionMap(),dt=P.dt)#,output_func = output_daily_and_final_incidence)#)
                            #integrator = init(prob,Tsit5();dt=1//2^(4),tstops=[0.5])
        #@time sims_vector=[init(prob,FunctionMap(),dt=P.dt,tstops=[365.])      for traj=1:n_traj]
        integrator=init(prob,FunctionMap(),dt=P.dt,tstops=[365.],#=save_idxs = [9,10,11],=#saveat=1:1:365)
        @time sims_vector=[ solve!(integrator).u     for traj=1:n_traj]
        #println(scenarios[i],"  :  ",size(sims_vector,1))
        #display(plot([sum(sims_vector[1][1][t]) for t=1:366]))
        @save folder*"sims_sc"*string(sc_nb)*".jld2" sims_vector
    end
end


function run_consensus_simulations_CT(folder,session,scenarios,β,ϵ,n_traj)
    if !isdir(folder)   mkdir(folder)   end
    global sims_vector=[]
    for i=1:size(scenarios,1)
        sc_nb=session*100+scenarios[i]
        u0,P,P_dest = KenyaCoV_CT.model_ingredients_from_data("data/data_for_age_structuredmodel.jld2","data/flight_numbers.csv","data/projected_global_prevelance.csv")
        P.χ = ones(KenyaCoV_CT.n_a)
        @load "data/detection_rates_for_different_taus.jld2" d_1
        P.rel_detection_rate = d_1
        P.dt = 0.25;        P.ext_inf_rate = 0.;        P.ϵ = ϵ#1.#Asymp are as infectious as symptomatics
        #Set the susceptibility vector --- just to specify the correct R₀
        sus_matrix = repeat(P.χ,1,17)
        R_A = P.ϵ*((1/P.σ₂) + (1/P.γ) ) #effective duration of asymptomatic
        R_M = (P.ϵ/P.σ₂) + (P.ϵ_D/P.γ) #effective duration of mild
        R_V = (P.ϵ/P.σ₂) + (P.ϵ_V/P.τ) #effective duration of severe
        R_vector = [(1-P.rel_detection_rate[a])*R_A + P.rel_detection_rate[a]*(1-P.hₐ[a])*R_M + P.rel_detection_rate[a]*P.hₐ[a]*R_V for a = 1:17]
        inf_matrix = repeat(R_vector',17,1)

        eigs, = eigen(sus_matrix.*P.M.*inf_matrix)
        max_eigval = Real(eigs[end])
        P.χ = ones(KenyaCoV_CT.n_a)/max_eigval
        P.β = β #rand(KenyaCoV.d_R₀) #Choose R₀ randomly from 2-3 range
        u0[KenyaCoV_CT.ind_nairobi_as,8,3] = 5      #Initial infecteds (P) in Nairobi in the 30-34 age group
        prob = KenyaCoV_CT.create_KenyaCoV_non_neg_prob(u0,(0.,365.),P)
        print("Simulating session ",session," sc ",sc_nb,"   -   ")
        ensemble_prob = EnsembleProblem(prob#=, prob_func = consensus_randomise_params=#, output_func = output_daily_and_final_incidence)
        return solve(ensemble_prob,FunctionMap(),dt = P.dt,callback = KenyaCoV_CT.callback_CT,trajectories = n_traj)
        sims_vector=[sims.u[i]      for i=1:size(sims.u,1)]
        @save folder*"sims_sc"*string(sc_nb)*".jld2" sims_vector
    end
end
