
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
    cumIs = [sum(sol(t)[:,:,7:9],dims = 2)[:,1,:]  for t in times] # only cumulative IA ID and IQ
    Q=[[sum(sol(t)[wa,:,5],dims = 1)[:,1,:] + sum(sol(t)[wa,:,11],dims = 1)[:,1,:]  for wa=1:20]  for t in times]      # Q = Q + QS
    I=[[sum(sol(t)[wa,:,3],dims = 1)[:,1,:] + sum(sol(t)[wa,:,4],dims = 1)[:,1,:] + sum(sol(t)[wa,:,10],dims = 1)[:,1,:]  for wa=1:20]  for t in times]      # I = IA + ID + IQ
    ###to delete when running large sims
    S=[[sum(sol(t)[wa,:,1],dims = 1)[:,1,:]  for wa=1:20]  for t in times]      # S
    R=[[sum(sol(t)[wa,:,6],dims = 1)[:,1,:]  for wa=1:20]  for t in times]      # R
    IQ=[[sum(sol(t)[wa,:,10],dims = 1)[:,1,:]  for wa=1:20]  for t in times]
    cumIC=[[sum(sol(t)[wa,:,12],dims = 1)[:,1,:]  for wa=1:20]  for t in times] # cumulative infected isolated after contact
    cumC=[[sum(sol(t)[wa,:,13],dims = 1)[:,1,:]  for wa=1:20]  for t in times] # cumulative contacts (all) isolated after contact
    #cumDeaths=[[sum(sol(t)[wa,:,14],dims = 1)[:,1,:]  for wa=1:20]  for t in times] # cumulative contacts (all) isolated after contact
    return [cumIs,cumC,cumIC,I,IQ,Q,S,R,sol[end][:,:,7:9]],false # save z (time series with only incidence with no age structure), and save the final distribution (end) age/space but no time
end

###########
function run_one_scenario(folder,ϵ,σ,γ,δ,β,r_R₀,first_sc_nb,τₚ,n_traj,CT_Imin_list,CT_dur_list,dt)
    if !isdir(folder)   mkdir(folder)   end
    u0,P,P_dest = KenyaCoV_CT.model_ingredients_from_data("data/data_for_age_structuredmodel.jld2","data/flight_numbers.csv","data/projected_global_prevelance.csv")
    u0[KenyaCoV_CT.ind_nairobi_as,5,4] = 5#Five initial infecteds in Nairobi in the 20-24 age group
    P.dt = dt;     P.ext_inf_rate = 0.;    P.ϵ = ϵ;        P.δ = δ;      P.γ = γ;      P.σ = σ;    P.β = β
    P.τ=1/2.;        P.κ=12;        P.κₘ=10;     P.Δₜ=10;   P.κ_per_event4=100; P.IDs_cfirst=true
    for wa=1:KenyaCoV_CT.n_a, a=1:KenyaCoV_CT.n_a       P.Mₚ[wa,a]=P.M[wa,a]/sum(P.M[wa,:]);    end

    prob = KenyaCoV_CT.create_KenyaCoV_non_neg_prob(u0,(0.,365.),P)
    sc_nb=first_sc_nb
      P.τₚ =τₚ;  P.CT_Imin=CT_Imin_list[1];P.CT_dur=CT_dur_list[1]
        println("\nSession "*string(sc_nb)*" Running ",n_traj," sims for τₚ,IDs_cfirst=",P.τₚ,",",P.IDs_cfirst)
        @time sims = solve(EnsembleProblem(prob#=,prob_func=randomise_params=#,output_func = output_daily_and_final_incidence),
                            FunctionMap(),dt=P.dt,trajectories=n_traj)
        println()
        global sims_vector=[]
        for i=1:size(sims.u,1)
            push!(sims_vector, sims.u[i])
        end
        sessionParams=KenyaCoV_CT.SessionParams(sc_nb=sc_nb,n_traj=n_traj,R₀=r_R₀,τₚ=τₚ,κ_per_event4=100,IDs_cfirst=P.IDs_cfirst,
                                    dt=dt,ext_inf_rate=P.ext_inf_rate,ϵ=ϵ,δ=δ,γ=γ,σ=σ,β=β,τ=1/2.,κ=12,κₘ=10,Δₜ=10)
        @save folder*"sims_sc"*string(sc_nb)*".jld2" sessionParams sims_vector
        sc_nb+=1
end
function run_set_scenarios(folder,ϵ,σ,γ,δ,β,r_R₀,first_sc_nb,τₚ,n_traj,CT_Imin,CT_dur_list,dt)
    if !isdir(folder)   mkdir(folder)   end
    sc_nb=first_sc_nb
    for CT_dur in CT_dur_list
        u0,P,P_dest = KenyaCoV_CT.model_ingredients_from_data("data/data_for_age_structuredmodel.jld2","data/flight_numbers.csv","data/projected_global_prevelance.csv")
        u0[KenyaCoV_CT.ind_nairobi_as,5,4] = 5#Five initial infecteds in Nairobi in the 20-24 age group
        P.dt = dt;     P.ext_inf_rate = 0.;    P.ϵ = ϵ;        P.δ = δ;      P.γ = γ;      P.σ = σ;    P.β = β
        P.τ=1/2.;        P.κ=12;        P.κₘ=10;     P.Δₜ=10;   P.κ_per_event4=100; P.IDs_cfirst=true
        P.CT_Imin=CT_Imin;  P.τₚ =τₚ;
        for wa=1:KenyaCoV_CT.n_a, a=1:KenyaCoV_CT.n_a       P.Mₚ[wa,a]=P.M[wa,a]/sum(P.M[wa,:]);    end

        prob = KenyaCoV_CT.create_KenyaCoV_non_neg_prob(u0,(0.,365.),P)

        P.CT_dur=CT_dur
        println("\nSession "*string(sc_nb)*" Running ",n_traj," sims for τₚ,CT_dur=",P.τₚ,",",P.CT_dur)
        @time sims = solve(EnsembleProblem(prob#=,prob_func=randomise_params=#,output_func = output_daily_and_final_incidence),
                            FunctionMap(),dt=P.dt,trajectories=n_traj)
        global sims_vector=[]
        for i=1:size(sims.u,1)
            push!(sims_vector, sims.u[i])
        end
        sessionParams=KenyaCoV_CT.SessionParams(sc_nb=sc_nb,n_traj=n_traj,R₀=r_R₀,τₚ=τₚ,κ_per_event4=100,IDs_cfirst=P.IDs_cfirst,
                                    dt=dt,ext_inf_rate=P.ext_inf_rate,ϵ=ϵ,δ=δ,γ=γ,σ=σ,β=β,τ=1/2.,κ=12,κₘ=10,Δₜ=10,CT_Imin=P.CT_Imin,CT_dur=P.CT_dur)
        @save folder*"sims_sc"*string(sc_nb)*".jld2" sessionParams sims_vector
        sc_nb+=1
    end
end
