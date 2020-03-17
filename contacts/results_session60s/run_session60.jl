push!(LOAD_PATH, "./src")
push!(LOAD_PATH, "./contacts")
using Plots,Parameters,Distributions,DifferentialEquations,JLD2,DataFrames,StatsPlots,FileIO,MAT,RecursiveArrayTools
import KenyaCoV_contacts
using LinearAlgebra:eigen
using Statistics: median, quantile

function randomise_params(prob,i,repeat)
    _P = deepcopy(prob.p)
    #_P.σ = 1/rand(d_incubation)
    #_P.β = rand(d_R₀)*_P.γ/(_P.δ + _P.ϵ*(1-_P.δ))
    return remake(prob,p=_P)
end
function run_saveMAT(P::KenyaCoV_contacts.CoVParameters_AS,prob,n_traj,τₚ_list,session)
    folder="./contacts/results_session60s/results_session"*string(session)*"/"
    #folder="W:/BACKUP MakingContacts/2020-03-16_v5/results_session"*string(session)*"/"
    mkdir(folder)
    for τₚ in τₚ_list
        P.τₚ =τₚ
        println("Session"*string(session)*" Running ",n_traj," sims for τₚ=",τₚ)
        @time sims = solve(EnsembleProblem(prob#=,prob_func=randomise_params=##=,output_func = output_daily_and_final_incidence=#),
                            FunctionMap(),dt=P.dt,trajectories=n_traj)
        sims_vector=sims.u[1].t
        for i=1:size(sims.u,1)
            sims_vector=[sims_vector sims.u[i].u]
        end
        file = matopen(folder*"sims"*string(n_traj)*"_taup"*string(τₚ)*".mat", "w")
        write(file, "sims", sims_vector)
        close(file)
    end
end

function run_saveJLD2daily_and_final_incidence(P::KenyaCoV_contacts.CoVParameters_AS,prob,n_traj,τₚ_list,session)
    folder="./contacts/results_session60s/results_session"*string(session)*"/"
    #folder="W:/BACKUP MakingContacts/2020-03-16_v5/results_session"*string(session)*"/"
    mkdir(folder)
    for τₚ in τₚ_list
        P.τₚ =τₚ
        println("Session"*string(session)*" Running ",n_traj," sims for τₚ=",τₚ)
        global sims = solve(EnsembleProblem(prob#=,prob_func=randomise_params=#,output_func = output_daily_and_final_incidence),
                            FunctionMap(),dt=P.dt,trajectories=n_traj)
        sims_vector=[]#sims.u[1].t
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
    #S=[[sum(sol(t)[wa,:,1],dims = 1)[:,1,:]  for wa=1:20]  for t in times]      # S
    #R=[[sum(sol(t)[wa,:,6],dims = 1)[:,1,:]  for wa=1:20]  for t in times]      # R

    return [cumulatives, Q, I, sol[end][:,:,7:9]#=,S,R=#],false # save z (time series with only incidence with no age structure), and save the final distribution (end) age/space but no time
end

"""
Define uncertainty of parameter estimates
"""
d_incubation = LogNormal(log(4.8),0.25) #Liu et al
mean(d_incubation)
(quantile(d_incubation,0.025),median(d_incubation),quantile(d_incubation,0.975))
d_R₀ = Gamma(100,2.92/100) ##Liu et al
mean(d_R₀)
(quantile(d_R₀,0.025),median(d_R₀),quantile(d_R₀,0.975))



n_traj=1000
"""
SCENARIOS
"""
u0,P,P_dest = KenyaCoV_contacts.model_ingredients_from_data("data/data_for_age_structuredmodel.jld2","data/flight_numbers.csv","data/projected_global_prevelance.csv")
u0[KenyaCoV_contacts.ind_nairobi_as,5,4] = 5#Five initial infecteds in Nairobi in the 20-24 age group
P.dt = 0.5;     P.ext_inf_rate = 0.;    P.ϵ = rand(Uniform(0.,0.5))
P.δ = 0.2;      P.γ = 1/2.5;            P.σ = 1/rand(d_incubation)
P.β = rand(d_R₀)*P.γ/(P.δ + P.ϵ*(1-P.δ))
P.τ=1/3.
P.κ=10
P.κₘ=7;     P.Δₜ=7
P.κ_per_event4=50
τₚ_list=[.0,.25,.5,.75,.9]
#P.Κ_max_capacity=[1000 for e in P.Κ_max_capacity]
#P.Κ_max_capacity[KenyaCoV_contacts.ind_nairobi_as]=1e3
#P.Κ_max_capacity[KenyaCoV_contacts.ind_mombasa_as]=1e2
P.Κ_max_capacity[12]=1e3
session_nb=60
P.κ_per_event4=50

for wa=1:KenyaCoV_contacts.n_a, a=1:KenyaCoV_contacts.n_a       P.Mₚ[wa,a]=P.M[wa,a]/sum(P.M[wa,:]);    end
prob = KenyaCoV_contacts.create_KenyaCoV_non_neg_prob(u0,(0.,365.),P)
results_sessions = run_saveJLD2daily_and_final_incidence(P,prob,n_traj,τₚ_list,session_nb)

#########
u0,P,P_dest = KenyaCoV_contacts.model_ingredients_from_data("data/data_for_age_structuredmodel.jld2","data/flight_numbers.csv","data/projected_global_prevelance.csv")
u0[KenyaCoV_contacts.ind_nairobi_as,5,4] = 5#Five initial infecteds in Nairobi in the 20-24 age group
P.dt = 0.5;     P.ext_inf_rate = 0.;    P.ϵ = rand(Uniform(0.,0.5))
P.δ = 0.2;      P.γ = 1/2.5;            P.σ = 1/rand(d_incubation)
P.β = rand(d_R₀)*P.γ/(P.δ + P.ϵ*(1-P.δ))
P.τ=1/3.
P.κ=10
P.κₘ=7;     P.Δₜ=7
P.κ_per_event4=50
τₚ_list=[.25,.5,.75,.9]

P.Κ_max_capacity[12]=5e3
session_nb=61
P.κ_per_event4=50

for wa=1:KenyaCoV_contacts.n_a, a=1:KenyaCoV_contacts.n_a       P.Mₚ[wa,a]=P.M[wa,a]/sum(P.M[wa,:]);    end
prob = KenyaCoV_contacts.create_KenyaCoV_non_neg_prob(u0,(0.,365.),P)
results_sessions = run_saveJLD2daily_and_final_incidence(P,prob,n_traj,τₚ_list,session_nb)

############
u0,P,P_dest = KenyaCoV_contacts.model_ingredients_from_data("data/data_for_age_structuredmodel.jld2","data/flight_numbers.csv","data/projected_global_prevelance.csv")
u0[KenyaCoV_contacts.ind_nairobi_as,5,4] = 5#Five initial infecteds in Nairobi in the 20-24 age group
P.dt = 0.5;     P.ext_inf_rate = 0.;    P.ϵ = rand(Uniform(0.,0.5))
P.δ = 0.2;      P.γ = 1/2.5;            P.σ = 1/rand(d_incubation)
P.β = rand(d_R₀)*P.γ/(P.δ + P.ϵ*(1-P.δ))
P.τ=1/3.
P.κ=10
P.κₘ=7;     P.Δₜ=7
P.κ_per_event4=50
τₚ_list=[.25,.5,.75,.9]

P.Κ_max_capacity[12]=1e4
session_nb=62
P.κ_per_event4=50

for wa=1:KenyaCoV_contacts.n_a, a=1:KenyaCoV_contacts.n_a       P.Mₚ[wa,a]=P.M[wa,a]/sum(P.M[wa,:]);    end
prob = KenyaCoV_contacts.create_KenyaCoV_non_neg_prob(u0,(0.,365.),P)
results_sessions = run_saveJLD2daily_and_final_incidence(P,prob,n_traj,τₚ_list,session_nb)

##########

###########
u0,P,P_dest = KenyaCoV_contacts.model_ingredients_from_data("data/data_for_age_structuredmodel.jld2","data/flight_numbers.csv","data/projected_global_prevelance.csv")
u0[KenyaCoV_contacts.ind_nairobi_as,5,4] = 5#Five initial infecteds in Nairobi in the 20-24 age group
P.dt = 0.5;     P.ext_inf_rate = 0.;    P.ϵ = rand(Uniform(0.,0.5))
P.δ = 0.2;      P.γ = 1/2.5;            P.σ = 1/rand(d_incubation)
P.β = rand(d_R₀)*P.γ/(P.δ + P.ϵ*(1-P.δ))
P.τ=1/3.
P.κ=10
P.κₘ=7;     P.Δₜ=7
P.κ_per_event4=50
τₚ_list=[.25,.5,.75,.9]
#P.Κ_max_capacity=[1000 for e in P.Κ_max_capacity]
#P.Κ_max_capacity[KenyaCoV_contacts.ind_nairobi_as]=1e3
#P.Κ_max_capacity[KenyaCoV_contacts.ind_mombasa_as]=1e2
P.Κ_max_capacity[12]=1e3
session_nb=65
P.κ_per_event4=100

for wa=1:KenyaCoV_contacts.n_a, a=1:KenyaCoV_contacts.n_a       P.Mₚ[wa,a]=P.M[wa,a]/sum(P.M[wa,:]);    end
prob = KenyaCoV_contacts.create_KenyaCoV_non_neg_prob(u0,(0.,365.),P)
results_sessions = run_saveJLD2daily_and_final_incidence(P,prob,n_traj,τₚ_list,session_nb)

#########
u0,P,P_dest = KenyaCoV_contacts.model_ingredients_from_data("data/data_for_age_structuredmodel.jld2","data/flight_numbers.csv","data/projected_global_prevelance.csv")
u0[KenyaCoV_contacts.ind_nairobi_as,5,4] = 5#Five initial infecteds in Nairobi in the 20-24 age group
P.dt = 0.5;     P.ext_inf_rate = 0.;    P.ϵ = rand(Uniform(0.,0.5))
P.δ = 0.2;      P.γ = 1/2.5;            P.σ = 1/rand(d_incubation)
P.β = rand(d_R₀)*P.γ/(P.δ + P.ϵ*(1-P.δ))
P.τ=1/3.
P.κ=10
P.κₘ=7;     P.Δₜ=7
P.κ_per_event4=50
τₚ_list=[.25,.5,.75,.9]

P.Κ_max_capacity[12]=5e3
session_nb=66
P.κ_per_event4=100

for wa=1:KenyaCoV_contacts.n_a, a=1:KenyaCoV_contacts.n_a       P.Mₚ[wa,a]=P.M[wa,a]/sum(P.M[wa,:]);    end
prob = KenyaCoV_contacts.create_KenyaCoV_non_neg_prob(u0,(0.,365.),P)
results_sessions = run_saveJLD2daily_and_final_incidence(P,prob,n_traj,τₚ_list,session_nb)

############
u0,P,P_dest = KenyaCoV_contacts.model_ingredients_from_data("data/data_for_age_structuredmodel.jld2","data/flight_numbers.csv","data/projected_global_prevelance.csv")
u0[KenyaCoV_contacts.ind_nairobi_as,5,4] = 5#Five initial infecteds in Nairobi in the 20-24 age group
P.dt = 0.5;     P.ext_inf_rate = 0.;    P.ϵ = rand(Uniform(0.,0.5))
P.δ = 0.2;      P.γ = 1/2.5;            P.σ = 1/rand(d_incubation)
P.β = rand(d_R₀)*P.γ/(P.δ + P.ϵ*(1-P.δ))
P.τ=1/3.
P.κ=10
P.κₘ=7;     P.Δₜ=7
P.κ_per_event4=50
τₚ_list=[.25,.5,.75,.9]

P.Κ_max_capacity[12]=1e4
session_nb=67
P.κ_per_event4=100

for wa=1:KenyaCoV_contacts.n_a, a=1:KenyaCoV_contacts.n_a       P.Mₚ[wa,a]=P.M[wa,a]/sum(P.M[wa,:]);    end
prob = KenyaCoV_contacts.create_KenyaCoV_non_neg_prob(u0,(0.,365.),P)
results_sessions = run_saveJLD2daily_and_final_incidence(P,prob,n_traj,τₚ_list,session_nb)
