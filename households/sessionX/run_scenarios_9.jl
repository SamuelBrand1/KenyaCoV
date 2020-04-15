push!(LOAD_PATH, "./src")
push!(LOAD_PATH, "./contacts")
include("../contacttracing_funtionswithdeaths.jl")

###########
r_R₀=2.5; δ=.05
ϵ = .3; σ = .2; γ = 1/2.5; β = r_R₀*γ/(δ + ϵ*(1-δ))
println("\nR₀,ϵ,δ=",r_R₀,",",ϵ,",",δ,"\n")

###########
session=30
scenarios=[90,93,96,99]
n_traj=50

α=[.9 for i=1:20]
CT_Imin_list=[[5 for j=1:20] for i=1:size(scenarios,1)];
CT_dur_list=[[30*(i-1) for j=1:20] for i=1:size(scenarios,1)];
CT_delay_list=[zeros(20) for j=1:size(scenarios,1)];

folder="./contacts/session"*string(session)*"_"*string(n_traj)*"sims/";
run_set_scenarios(folder,session,scenarios,ϵ,σ,γ,δ,β,r_R₀,α,n_traj,CT_Imin_list,CT_dur_list,CT_delay_list,.5)
