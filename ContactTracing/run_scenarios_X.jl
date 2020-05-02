push!(LOAD_PATH, "./src")
    push!(LOAD_PATH, "./ContactTracing")
    include("forecast_functions_CT.jl")

###########
β = 2.5
ϵ = 1.
session=01
scenarios=[1,2]#[i  for i=1:90]#[8]
n_traj=5
nb_months=1

folder="./ContactTracing/session"*string(session)*"_"*string(n_traj)*"sims/";
run_set_scenarios(folder,session,scenarios,β,ϵ,n_traj)

p=plot();for s=9:11   plot!(p,[sum(sims_vector[1][t][:,:,s]) for t=1:365]#=,ylims=(1,1e5)=#)  end;display(p)
println("cumI  =",sum(sims_vector[1][end][:,:,9:11]))
    println("Asympt=",sum(sims_vector[1][end][:,:,9]))
    println("Sympt =",sum(sims_vector[1][end][:,:,10]))
    println("Severe=",sum(sims_vector[1][end][:,:,11]))
    println("H=",sum(sims_vector[1][end][:,:,12]))

plot([sims_vector[1][t][:,:,2] for t=1:365])
for i in sims_vector[1]
  display(plot!(sims_vector[1].u,vars=(0,1),legend=false))
end
