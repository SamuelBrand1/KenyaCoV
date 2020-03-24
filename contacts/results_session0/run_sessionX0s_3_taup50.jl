push!(LOAD_PATH, "./src")
push!(LOAD_PATH, "./contacts")
using Plots,Parameters,Distributions,DifferentialEquations,JLD2,DataFrames,StatsPlots,FileIO,MAT,RecursiveArrayTools
import KenyaCoV_contacts
using LinearAlgebra:eigen
using Statistics: median, quantile

###########
function randomise_params(prob,i,repeat)
    _P = deepcopy(prob.p)
    #_P.σ = 1/rand(d_incubation)
    #_P.β = rand(d_R₀)*_P.γ/(_P.δ + _P.ϵ*(1-_P.δ))
    return remake(prob,p=_P)
end
function run_saveJLD2daily_and_final_incidence(P::KenyaCoV_contacts.CoVParameters_AS,prob,n_traj,τₚ_list,session)
    folder="./contacts/results_session"*string(Int(floor(session/10)))*"0s/results_session"*string(session)*"/"
    #folder="W:/BACKUP MakingContacts/2020-03-16_v5/results_session"*string(session)*"/"
    if !isdir(folder)   mkdir(folder)   end
    i=0
    #τₚ=τₚ_list ###############################
    for τₚ in τₚ_list
        #i+=1
        P.τₚ =τₚ
        println("\nSession"*string(session)*" Running ",n_traj," sims for τₚ=",τₚ)
        @time sims = solve(EnsembleProblem(prob#=,prob_func=randomise_params=#,output_func = output_daily_and_final_incidence),
                            FunctionMap(),dt=P.dt,trajectories=n_traj)
        sims_vector=[]
        for i=1:size(sims.u,1)
            push!(sims_vector, sims.u[i])
        end
        @save folder*"sims"*string(n_traj)*"_taup"*string(τₚ)*".jld2" sims_vector
    end
end
function output_daily_and_final_incidence(sol,i)
    times = sol.prob.tspan[1]:1:sol.prob.tspan[end] #time step = 1 day
    cumulatives = [sum(sol(t)[:,:,7:9],dims = 2)[:,1,:]  for t in times] # only cumulative IA ID and IQ
    Q=[[sum(sol(t)[wa,:,5],dims = 1)[:,1,:] + sum(sol(t)[wa,:,11],dims = 1)[:,1,:]  for wa=1:20]  for t in times]      # Q = Q + QS
    I=[[sum(sol(t)[wa,:,3],dims = 1)[:,1,:] + sum(sol(t)[wa,:,4],dims = 1)[:,1,:] + sum(sol(t)[wa,:,10],dims = 1)[:,1,:]  for wa=1:20]  for t in times]      # I = IA + ID + IQ

    ###to delete when running large sims
    S=[[sum(sol(t)[wa,:,1],dims = 1)[:,1,:]  for wa=1:20]  for t in times]      # S
    R=[[sum(sol(t)[wa,:,6],dims = 1)[:,1,:]  for wa=1:20]  for t in times]      # R
    cumIC=[[sum(sol(t)[wa,:,12],dims = 1)[:,1,:]  for wa=1:20]  for t in times] # cumulative infected isolated after contact

    return [cumulatives, Q, I, sol[end][:,:,7:9],S,R,cumIC],false # save z (time series with only incidence with no age structure), and save the final distribution (end) age/space but no time
end

###########
function run_session(session,Κ_max_capacity_KENYA,Κ_max_capacity_Nairobi,Κ_max_capacity_Kilifi,κ_per_event4,τₚ_list,n_traj,ϵ,σ,γ,δ,β,Δₜ,stop_Q,IDs_contacted_first)
    u0,P,P_dest = KenyaCoV_contacts.model_ingredients_from_data("data/data_for_age_structuredmodel.jld2","data/flight_numbers.csv","data/projected_global_prevelance.csv")
    u0[KenyaCoV_contacts.ind_nairobi_as,5,4] = 5#Five initial infecteds in Nairobi in the 20-24 age group
    P.dt = 0.5;     P.ext_inf_rate = 0.;    P.ϵ = ϵ;        P.δ = δ;      P.γ = γ;      P.σ = σ;    P.β = β
    P.τ=1/2.;        P.κ=20;        P.κₘ=Δₜ;     P.Δₜ=Δₜ;     P.stop_Q=stop_Q;      P.IDs_contacted_first=IDs_contacted_first;
    for wa=1:KenyaCoV_contacts.n_a, a=1:KenyaCoV_contacts.n_a       P.Mₚ[wa,a]=P.M[wa,a]/sum(P.M[wa,:]);    end

    P.Κ_max_capacity=[Κ_max_capacity_KENYA for e in P.Κ_max_capacity]
    P.Κ_max_capacity[4]=Κ_max_capacity_Nairobi;        P.Κ_max_capacity[12]=Κ_max_capacity_Kilifi
    P.κ_per_event4=κ_per_event4

    prob = KenyaCoV_contacts.create_KenyaCoV_non_neg_prob(u0,(0.,365.),P)
    #run_saveJLD2daily_and_final_incidence(P,prob,n_traj,τₚ_list,session)
    run_combined_scenarios(P,prob,n_traj,τₚ_list,session)
end


function run_combined_scenarios(folder,ϵ,σ,γ,δ,β,R₀,Κ_max_capacity_KENYA,Κ_max_capacity4,first_sc_nb,τₚ_list,stop_Q_list,Κ_max_capacity12_list,κ_per_event4_list,IDs_cfirst_list,n_traj,dt)
    if !isdir(folder)   mkdir(folder)   end
    u0,P,P_dest = KenyaCoV_contacts.model_ingredients_from_data("data/data_for_age_structuredmodel.jld2","data/flight_numbers.csv","data/projected_global_prevelance.csv")
    u0[KenyaCoV_contacts.ind_nairobi_as,5,4] = 5#Five initial infecteds in Nairobi in the 20-24 age group
    P.dt = dt;     P.ext_inf_rate = 0.;    P.ϵ = ϵ;        P.δ = δ;      P.γ = γ;      P.σ = σ;    P.β = β #β = r_R₀*γ/(δ + ϵ*(1-δ))
    P.τ=1/2.;        P.κ=12;        P.κₘ=10;     P.Δₜ=10;
    for wa=1:KenyaCoV_contacts.n_a, a=1:KenyaCoV_contacts.n_a       P.Mₚ[wa,a]=P.M[wa,a]/sum(P.M[wa,:]);    end
    P.Κ_max_capacity=[Κ_max_capacity_KENYA for e in P.Κ_max_capacity];      P.Κ_max_capacity[4]=Κ_max_capacity4;


    prob = KenyaCoV_contacts.create_KenyaCoV_non_neg_prob(u0,(0.,365.),P)
    sc_nb=first_sc_nb
    for stop_Q in stop_Q_list, τₚ in τₚ_list, Κ_max_capacity12 in Κ_max_capacity12_list, κ_per_event4 in κ_per_event4_list, IDs_cfirst in IDs_cfirst_list
        P.stop_Q=stop_Q;  P.τₚ =τₚ;  P.Κ_max_capacity[12]=Κ_max_capacity12;  P.κ_per_event4=κ_per_event4;  P.IDs_cfirst=IDs_cfirst
        println("\nSession "*string(sc_nb)*" Running ",n_traj," sims for τₚ,stop_Q,Κ_max_capacity12,κ_per_event4,IDs_cfirst=",P.τₚ,",",P.stop_Q,",",P.Κ_max_capacity[12],",",P.κ_per_event4,",",P.IDs_cfirst)
        @time sims = solve(EnsembleProblem(prob#=,prob_func=randomise_params=#,output_func = output_daily_and_final_incidence),
                            FunctionMap(),dt=P.dt,trajectories=n_traj)
        sims_vector=[]
        for i=1:size(sims.u,1)
            push!(sims_vector, sims.u[i])
        end
        sessionParams=KenyaCoV_contacts.SessionParams(sc_nb=sc_nb,n_traj=n_traj,τₚ=τₚ,stop_Q=stop_Q,Κ_max_capacity12=Κ_max_capacity12,κ_per_event4=κ_per_event4,IDs_cfirst=IDs_cfirst,
                                    dt=dt,ext_inf_rate=P.ext_inf_rate,ϵ=ϵ,δ=δ,γ=γ,σ=σ,β=β,τ=1/2.,κ=12,κₘ=10,Δₜ=10,Κ_max_capacity=Κ_max_capacity_KENYA,Κ_max_capacity4=Κ_max_capacity4)
        @save folder*"sims_sc"*string(sc_nb)*".jld2" sessionParams sims_vector
        #@save folder*"sims"*string(n_traj)*"_sc"*string(sc_nb)*".jld2" sims_vector
        #=file = matopen(folder*"sims"*string(n_traj)*"_sc"*string(sc_nb)*".mat","w");
        write(file, "sessionParams", [sc_nb,n_traj,τₚ,stop_Q,Κ_max_capacity12,κ_per_event4,IDs_cfirst,dt,P.ext_inf_rate,ϵ,δ,γ,σ,β,1/2.,12,10,10,Κ_max_capacity_KENYA,Κ_max_capacity4]);
        write(file, "sims_vector", sims_vector);close(file)=#
        sc_nb+=1
    end
end
###########
d_incubation = LogNormal(log(4.8),0.25) #Liu et al
d_R₀ = Gamma(100,2.92/100) ##Liu et al
r_R₀=3
ϵ = .3#rand(Uniform(0.,0.5))
σ = .2#1/rand(d_incubation)
γ = 1/2.5
δ = 0.2
β = r_R₀*γ/(δ + ϵ*(1-δ))#β = rand(d_R₀)*γ/(δ + ϵ*(1-δ))
println("\nϵ,σ,γ,δ,β=",ϵ,",",σ,",",γ,",",δ,",",β,"\n")

###########
Κ_max_capacity_KENYA=0; Κ_max_capacity4=0#Κ_max_capacity_Nairobi
τₚ_list=[.5]
stop_Q_list=[true,false]    #****5 do we stop detecting after we stopped tracing
Κ_max_capacity12_list=[0,1e3,1e4,1e5,1e6]#Κ_max_capacity_Kilifi_list
κ_per_event4_list=[50,100,200]
IDs_cfirst_list=[true,false]
n_traj=50

folder="./contacts/results_session0/";    first_sc_nb=500
run_combined_scenarios(folder,ϵ,σ,γ,δ,β,r_R₀,Κ_max_capacity_KENYA,Κ_max_capacity4,first_sc_nb,τₚ_list,stop_Q_list,Κ_max_capacity12_list,κ_per_event4_list,IDs_cfirst_list,n_traj,.5)
