
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
#Define initial conditions by modifying the completely susceptible population
P.dt = 0.25
u0[30,4,1] += 1#One diseased in Nairobi

prob = KenyaCoV.create_KenyaCoV_non_neg_prob(u0,(0.,365.),P)

@time sol = solve(prob,FunctionMap(),dt = P.dt)

f = findall(sol.u[end] .< 0 )
sol.u[end][f]
#Create a JumpProblem which you can solve --- needs DifferentialEquations module for the solvers
jump_prob_tl = create_KenyaCoV_prob(u0,(0.,365.),P)
#Go straight to solution using solver compiled in the KenyaCoV module
@time sol_tl = KenyaCoV.solve_KenyaCoV_prob(u0,(0.,365.),P,0.5)
[any(x .< 0) for x in sol_tl.u]
f = findall(sol_tl.u[end] .< 0 )
