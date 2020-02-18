
push!(LOAD_PATH, "/Users/Sam/GitHub/KenyaCoV/src")
using Plots,Parameters,Distributions,DifferentialEquations,JLD2,DataFrames
using Revise
import KenyaCoV



"""
Load a completely susceptible population differentiated by county and urban/rural.
Also, load the optimised mixing matrix for spatial transmission
"""
# #Load data and completely susceptible Population
# u0,P,transport_matrix = model_ingredients_from_data("./src/2009_National_Estimates_of_Rural_and_Urban_Populations_by_County.csv",
#Load data completely susceptible Population
u0,P,transport_matrix = KenyaCoV.model_ingredients_from_data("data/combined_population_estimates.csv",
                                                             "data/optimal_transition_matrix.jld2",
                                                            "data/optimal_movement_matrix.jld2",
                                                            "data/flight_numbers.csv",
                                                            "data//projected_global_prevelance.csv")
"""
Example of methods that modify underlying parameters
"""
#This method modifies the parameter set for changing the mixing structure
# KenyaCoV.transportstructure_params!(P,[0.001 for i = 1:KenyaCoV.n],transport_matrix)
#You can modify other parameters directly
P.τ = 0. #e.g. no treatment rate
for (i,p) in enumerate(P.global_prev)
    P.global_prev[i] = 0.
end
P.global_prev
P.μ₁ = 0.
#Define initial conditions by modifying the completely susceptible population
P.dt = 1.
u0[30,4,1] += 1#One diseased in Nairobi

prob = KenyaCoV.create_KenyaCoV_non_neg_prob(u0,(0.,365.),P)
prob_tl = KenyaCoV.create_KenyaCoV_prob(u0,(0.,365.),P)
@time sol_tl = solve(prob_tl,SimpleTauLeaping(),dt = P.dt)

@time sol = solve(prob,FunctionMap(),dt = P.dt)

f = findall(sol.u[end] .< 0 )
f_tl = findall(sol_tl.u[end] .<0)
