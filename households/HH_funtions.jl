
using Plots,Parameters,Distributions,DifferentialEquations,JLD2,DataFrames,StatsPlots,FileIO,MAT,RecursiveArrayTools
import HH_KenyaCoV
using LinearAlgebra:eigen
using Statistics: median, quantile

####### Matt's death probabilities
#=Sympt_2_hosp=[3.8 3.8 2.6 2.6 2.8 2.8 2.7 2.7 5.4 5.4 12.6 12.6 19.7 19.7 28.7 28.7 27.3 27.3 27.3 27.3 27.3] ./ 100
Sympt_2_critcal=[0.01 0.01 0.02 0.02 0.08 0.08 0.18 0.18 0.34 0.34 1.5 1.5 5.4 5.4 12.4 12.4 19.3 19.3 19.3 19.3 19.3] ./ 100
Hosp_2_Death=[3.8 3.8 3.8 3.8 3.8 3.8 3.8 4.0 4.5 5.6 7.8 11.3 16.9 23.2 29.1 34.8 53.5 53.5 53.5 53.5 53.5] ./ 100
μ₁=Sympt_2_hosp .+ Hosp_2_Death
#μ₁[1:16]=#
####### China's death rates  [http://weekly.chinacdc.cn/en/article/id/e53946e2-c6c4-41e9-9a9b-fea8db1a8f51]
μ₁=[.0,.0,.2,.2,.2,.2,.2,.2,.4,.4,1.3,1.3,3.6,3.6,8,#=8,14.8=#11.4] ./100 #11.4=(8+14.8)/2
μ₁=.2
###########
function randomise_params(prob,i,repeat)
    _P = deepcopy(prob.p)
    _P.σ = 1/rand(d_incubation)
    _P.β = rand(d_R₀)*_P.γ/(_P.δ + _P.ϵ*(1-_P.δ))
    return remake(prob,p=_P)
end
function output_daily_and_final_incidence(sol,i)
    times = sol.prob.tspan[1]:1:sol.prob.tspan[end] #time step = 1 day
    cumIs = [sum(sol(t)[:,:,7:8],dims = 2)[:,1,:]  for t in times] # only cumulative IA and ID
    Q=[[sum(sol(t)[wa,:,5],dims = 1)[:,1,:] #=+ sum(sol(t)[wa,:,11],dims = 1)[:,1,:]=#  for wa=1:20]  for t in times]      # Q = Q ##+ QS
    I=[[sum(sol(t)[wa,:,3],dims = 1)[:,1,:] + sum(sol(t)[wa,:,4],dims = 1)[:,1,:]  for wa=1:20]  for t in times]      # I = IA + ID
    ###to delete when running large sims
    S=[[sum(sol(t)[wa,:,1],dims = 1)[:,1,:]  for wa=1:20]  for t in times]      # S
    R=[[sum(sol(t)[wa,:,6],dims = 1)[:,1,:]  for wa=1:20]  for t in times]      # R
    #IQ=[[sum(sol(t)[wa,:,10],dims = 1)[:,1,:]  for wa=1:20]  for t in times]
    #cumIC=[[sum(sol(t)[wa,:,12],dims = 1)[:,1,:]  for wa=1:20]  for t in times] # cumulative infected isolated after contact
    #cumC=[[sum(sol(t)[wa,:,13],dims = 1)[:,1,:]  for wa=1:20]  for t in times] # cumulative contacts (all) isolated after contact
    cumDeaths=[[sum(sol(t)[wa,:,9],dims = 1)[:,1,:]  for wa=1:20]  for t in times] # cumulative contacts (all) isolated after contact
    return [cumIs,cumDeaths,I,Q,S,R,sol[end][:,:,7:8],sol[end][:,:,9]],false # save z (time series with only incidence with no age structure), and save the final distribution (end) age/space but no time
end

########### WITH CT_DELAY
#function run_set_scenarios(folder,session,scenarios,ϵ,σ,γ,δ,β,r_R₀,α,n_traj,CT_Imin_list,CT_dur_list,CT_delay_list,dt)
function run_set_scenarios(folder,session,scenarios,ϵ,σ,γ,δ,βˢ,βᶠ,r_R₀,n_traj,dt)
    if !isdir(folder)   mkdir(folder)   end
    for i=1:size(scenarios,1)#size(CT_dur_list,1)
        sc_nb=session*100+scenarios[i]
        u0,P,P_dest = HH_KenyaCoV.model_ingredients_from_data("data/data_for_age_structuredmodel.jld2","data/flight_numbers.csv","data/projected_global_prevelance.csv")
        u0[HH_KenyaCoV.ind_nairobi_as,5,4] = 5#Five initial infecteds in Nairobi in the 20-24 age group
        P.dt = dt;     P.ext_inf_rate = 0.;    P.ϵ = ϵ;        P.δ = δ;      P.γ = γ;      P.σ = σ;
        P.βˢ = βˢ;   P.βᶠ = βᶠ;
        P.τ=1/2.;       # P.κ=12;        P.κₘ=10;     P.Δₜ=10;   P.κ_per_event4=100; P.IDs_cfirst=true
        #P.CT_Imin=CT_Imin_list[i];  P.α =α;       P.μ₁=μ₁[1:KenyaCoV_CT.n_a];     P.CT_delay=CT_delay_list[i];  P.CT_dur=CT_dur_list[i]
        P.μ₁=μ₁
        #for wa=1:KenyaCoV_CT.n_a, a=1:KenyaCoV_CT.n_a       P.Mₚ[wa,a]=P.M[wa,a]/sum(P.M[wa,:]);    end

        prob = HH_KenyaCoV.create_KenyaCoV_non_neg_prob(u0,(0.,365.),P)


        #println("\nScenario "*string(sc_nb)*" Running ",n_traj," sims for α,CT_dur,CT_delay=",P.α,",",P.CT_dur,P.CT_delay)
        @time sims = solve(EnsembleProblem(prob#=,prob_func=randomise_params=#,output_func = output_daily_and_final_incidence),
                            FunctionMap(),dt=P.dt,trajectories=n_traj)
        sims_vector=[]
        for i=1:size(sims.u,1)
            push!(sims_vector, sims.u[i])
        end
        #sessionParams=KenyaCoV_CT.SessionParams(sc_nb=sc_nb,n_traj=n_traj,R₀=r_R₀,α=α,κ_per_event4=100,IDs_cfirst=P.IDs_cfirst,
        #                            dt=dt,ext_inf_rate=P.ext_inf_rate,ϵ=ϵ,δ=δ,γ=γ,σ=σ,β=β,τ=1/2.,κ=12,κₘ=10,Δₜ=10,CT_Imin=P.CT_Imin,CT_dur=P.CT_dur,CT_delay=P.CT_delay)
        @save folder*"sims_sc"*string(sc_nb)*".jld2" #=sessionParams=# sims_vector
    end
end
