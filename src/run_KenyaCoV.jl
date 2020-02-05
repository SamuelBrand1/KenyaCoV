push!(LOAD_PATH, "/Users/Sam/GitHub/KenyaCoV/src")
using Plots,Parameters,Distributions,DifferentialEquations,JLD2,DataFrames
using Revise
using KenyaCoV

"""
Load a completely susceptible population differentiated by county and urban/rural.
Also, load a default parameter set including the optimised mixing matrix (T), and time spent away (ρ)
The mixing matrix is related to the transport matrix by the choice of ρ, so I also load transport_matrix
so different choices of ρ can be tested
"""

#Load data completely susceptible Population
u0,P,transport_matrix = KenyaCoV.model_ingredients_from_data("data/combined_population_estimates.csv","data/optimal_transition_matrix.jld2","data/optimal_movement_matrix.jld2" )
"""
Example of methods that modify underlying parameters
"""
#This method modifies the parameter set for changing the mixing structure
# KenyaCoV.transportstructure_params!(P,[0.001 for i = 1:KenyaCoV.n],transport_matrix)
#You can modify other parameters directly
P.τ = 0. #e.g. no treatment rate

#Define initial conditions by modifying the completely susceptible population
u0[30,3,1] += 1#One asymptomatic in Nairobi

#Create a JumpProblem which you can solve --- needs DifferentialEquations module for the solvers
jump_prob_tl = create_KenyaCoV_prob(u0,(0.,365.),P)
#Go straight to solution using solver compiled in the KenyaCoV module
@time sol_tl = solve_KenyaCoV_prob(u0,(0.,60.),P,0.5)
