push!(LOAD_PATH, "/Users/Sam/GitHub/KenyaCoV/src")
using Plots,Parameters,Distributions,DifferentialEquations,JLD2,DataFrames
using Revise
using KenyaCoV

@load "data/optimal_transition_matrix.jld2" T_opt
@load "data/optimal_movement_matrix.jld2" P_opt
sum(P_opt)
#Load data and completely susceptible Population
u0,P,transport_matrix = KenyaCoV.model_ingredients_from_data("data/combined_population_estimates.csv","data/optimal_transition_matrix.jld2","data/optimal_movement_matrix.jld2" )
#This method modifies the parameter set for changing the movement structure
KenyaCoV.transportstructure_params!(P,[0.001 for i = 1:KenyaCoV.n],transport_matrix)
#Then you can modify other parameters
P.τ = 1/1. #e.g. increased treatment rate
#Define initial conditions by modifying the completely susceptible population
u0[30,3,1] += 1#One asymptomatic in Nairobi

#Create a JumpProblem which you can solve --- needs DifferentialEquations module for the solvers
jump_prob_tl = create_KenyaCoV_prob(u0,(0.,365.),P)
#Go straight to solution using solver compiled in the KenyaCoV module
@time sol_tl = solve_KenyaCoV_prob(u0,(0.,60.),P,0.5)
