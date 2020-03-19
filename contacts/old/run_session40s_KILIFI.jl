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
    #folder="./contacts/results_session40s/results_session"*string(session)*"/"
    #folder="W:/BACKUP MakingContacts/2020-03-16_v5/results_session"*string(session)*"/"
    folder="./contacts/results_session40s/results_session"*string(session)*"/"
    mkdir(folder)
    for τₚ in τₚ_list
        P.τₚ =τₚ
        println("Session"*string(session)*" Running ",n_traj," sims for τₚ=",τₚ)
        @time sims = solve(EnsembleProblem(prob#=,prob_func=randomise_params=#),FunctionMap(),dt=P.dt,trajectories=n_traj)
        sims_vector=sims.u[1].t
        for i=1:size(sims.u,1)
            sims_vector=[sims_vector sims.u[i].u]
        end
        file = matopen(folder*"sims"*string(n_traj)*"_taup"*string(τₚ)*".mat", "w")
        write(file, "sims", sims_vector)
        close(file)
    end
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

"""
SCENARIO I
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
τₚ_list=[.25,.5,.75,.9]#[.0,.25,.5,.75,.9]
#P.Κ_max_capacity=[1000 for e in P.Κ_max_capacity]
#P.Κ_max_capacity[KenyaCoV_contacts.ind_nairobi_as]=1e3
#P.Κ_max_capacity[KenyaCoV_contacts.ind_mombasa_as]=1e2
P.Κ_max_capacity[12]=1e3
session_nb=46
n_traj=50

for wa=1:KenyaCoV_contacts.n_a, a=1:KenyaCoV_contacts.n_a       P.Mₚ[wa,a]=P.M[wa,a]/sum(P.M[wa,:]);    end
prob = KenyaCoV_contacts.create_KenyaCoV_non_neg_prob(u0,(0.,365.),P)
results_sessions = run_saveMAT(P,prob,n_traj,τₚ_list,session_nb)

"""
SCENARIO II
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
τₚ_list=[.25,.5,.75,.9]#[.0,.25,.5,.75,.9]
#P.Κ_max_capacity=[1000 for e in P.Κ_max_capacity]
#P.Κ_max_capacity[KenyaCoV_contacts.ind_nairobi_as]=1e3
#P.Κ_max_capacity[KenyaCoV_contacts.ind_mombasa_as]=1e2
P.Κ_max_capacity[12]=5e3
session_nb=47
#n_traj=5

for wa=1:KenyaCoV_contacts.n_a, a=1:KenyaCoV_contacts.n_a       P.Mₚ[wa,a]=P.M[wa,a]/sum(P.M[wa,:]);    end
prob = KenyaCoV_contacts.create_KenyaCoV_non_neg_prob(u0,(0.,365.),P)
results_sessions = run_saveMAT(P,prob,n_traj,τₚ_list,session_nb)

"""
SCENARIO III
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
τₚ_list=[.25,.5,.75,.9]#[.0,.25,.5,.75,.9]
#P.Κ_max_capacity=[1000 for e in P.Κ_max_capacity]
#P.Κ_max_capacity[KenyaCoV_contacts.ind_nairobi_as]=1e3
#P.Κ_max_capacity[KenyaCoV_contacts.ind_mombasa_as]=1e2
P.Κ_max_capacity[12]=1e4
session_nb=48
#n_traj=5

for wa=1:KenyaCoV_contacts.n_a, a=1:KenyaCoV_contacts.n_a       P.Mₚ[wa,a]=P.M[wa,a]/sum(P.M[wa,:]);    end
prob = KenyaCoV_contacts.create_KenyaCoV_non_neg_prob(u0,(0.,365.),P)
results_sessions = run_saveMAT(P,prob,n_traj,τₚ_list,session_nb)
