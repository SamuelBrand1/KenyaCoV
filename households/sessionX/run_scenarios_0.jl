push!(LOAD_PATH, "./src")
push!(LOAD_PATH, "./households")
include("../HH_funtions.jl")

###########
r_R₀=2.5; δ=.05
ϵ = .3; σ = .2; γ = 1/2.5; β = r_R₀*γ/(δ + ϵ*(1-δ))
println("\nR₀,ϵ,δ=",r_R₀,",",ϵ,",",δ,"\n")

###########
session=43
scenarios=[0]#[50,53,56,59]
n_traj=5

βˢ=r_R₀*γ/(δ + ϵ*(1-δ))+.2#1.5
βᶠ=r_R₀*γ/(δ + ϵ*(1-δ))-.2


folder="./households/session"*string(session)*"_"*string(n_traj)*"sims/";
#run_set_scenarios(folder,session,scenarios,ϵ,σ,γ,δ,β,r_R₀,α,n_traj,CT_Imin_list,CT_dur_list,CT_delay_list,.5)
run_set_scenarios(folder,session,scenarios,ϵ,σ,γ,δ,βˢ,βᶠ,r_R₀,n_traj,.5)
