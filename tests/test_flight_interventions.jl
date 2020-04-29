# In this script we run the model for a month
# for a few different flight interventions.
# In each case we assume the prevalence in the
# external population is 1e-6. Before the intervention,
# we assume 500 daily visitors to Mombassa and
# 1000 daily visitors to Nairobi.

push!(LOAD_PATH, "C:/Users/Joseph Hilton/Documents/GitHub/KenyaCoV/src")
using Plots,Parameters,Distributions
using KenyaCoV

# First do case with no controls:
u0,P,transport_matrix = model_ingredients_from_data("./src/2009_National_Estimates_of_Rural_and_Urban_Populations_by_County.csv","./src/flight_numbers_no_controls.csv","./src/projected_global_prevelance.csv",0.01)
#This method modifies the parameter set for changing the movement structure
KenyaCoV.transportstructure_params!(P,0.001,transport_matrix)
#Then you can modify other parameters
P.τ = 1/1. #e.g. increased treatment rate
#Define initial conditions by modifying the completely susceptible population
u0[30,3,1] += 1#One asymptomatic in Nairobi

#Create a JumpProblem which you can solve --- needs DifferentialEquations module for the solvers
jump_prob_tl = create_KenyaCoV_prob(u0,(0.,360.),P)
#Go straight to solution using solver compiled in the KenyaCoV module
@time sol_tl_no_controls = solve_KenyaCoV_prob(u0,(0.,360.),P,0.25)

# Now do no flights whatsoever:
u0,P,transport_matrix = model_ingredients_from_data("./src/2009_National_Estimates_of_Rural_and_Urban_Populations_by_County.csv","./src/flight_numbers_complete_shutdown.csv","./src/projected_global_prevelance.csv",0.01)
#This method modifies the parameter set for changing the movement structure
KenyaCoV.transportstructure_params!(P,0.001,transport_matrix)
#Then you can modify other parameters
P.τ = 1/1. #e.g. increased treatment rate
#Define initial conditions by modifying the completely susceptible population
u0[30,3,1] += 1#One asymptomatic in Nairobi

#Create a JumpProblem which you can solve --- needs DifferentialEquations module for the solvers
jump_prob_tl = create_KenyaCoV_prob(u0,(0.,360.),P)
#Go straight to solution using solver compiled in the KenyaCoV module
@time sol_tl_complete_shutdown = solve_KenyaCoV_prob(u0,(0.,360.),P,0.25)

# Flights halved after one week:
u0,P,transport_matrix = model_ingredients_from_data("./src/2009_National_Estimates_of_Rural_and_Urban_Populations_by_County.csv","./src/flight_numbers_halved_after_week.csv","./src/projected_global_prevelance.csv",0.01)
#This method modifies the parameter set for changing the movement structure
KenyaCoV.transportstructure_params!(P,0.001,transport_matrix)
#Then you can modify other parameters
P.τ = 1/1. #e.g. increased treatment rate
#Define initial conditions by modifying the completely susceptible population
u0[30,3,1] += 1#One asymptomatic in Nairobi

#Create a JumpProblem which you can solve --- needs DifferentialEquations module for the solvers
jump_prob_tl = create_KenyaCoV_prob(u0,(0.,360.),P)
#Go straight to solution using solver compiled in the KenyaCoV module
@time sol_tl_halved_after_week = solve_KenyaCoV_prob(u0,(0.,360.),P,0.25)

# Flights quartered after one week:
u0,P,transport_matrix = model_ingredients_from_data("./src/2009_National_Estimates_of_Rural_and_Urban_Populations_by_County.csv","./src/flight_numbers_quartered_after_week.csv","./src/projected_global_prevelance.csv",0.01)
#This method modifies the parameter set for changing the movement structure
KenyaCoV.transportstructure_params!(P,0.001,transport_matrix)
#Then you can modify other parameters
P.τ = 1/1. #e.g. increased treatment rate
#Define initial conditions by modifying the completely susceptible population
u0[30,3,1] += 1#One asymptomatic in Nairobi

#Create a JumpProblem which you can solve --- needs DifferentialEquations module for the solvers
jump_prob_tl = create_KenyaCoV_prob(u0,(0.,360.),P)
#Go straight to solution using solver compiled in the KenyaCoV module
@time sol_tl_quartered_after_week = solve_KenyaCoV_prob(u0,(0.,360.),P,0.25)

# Flights stopped completely after one week:
u0,P,transport_matrix = model_ingredients_from_data("./src/2009_National_Estimates_of_Rural_and_Urban_Populations_by_County.csv","./src/flight_numbers_shutdown_after_week.csv","./src/projected_global_prevelance.csv",0.01)
#This method modifies the parameter set for changing the movement structure
KenyaCoV.transportstructure_params!(P,0.001,transport_matrix)
#Then you can modify other parameters
P.τ = 1/1. #e.g. increased treatment rate
#Define initial conditions by modifying the completely susceptible population
u0[30,3,1] += 1#One asymptomatic in Nairobi

#Create a JumpProblem which you can solve --- needs DifferentialEquations module for the solvers
jump_prob_tl = create_KenyaCoV_prob(u0,(0.,360.),P)
#Go straight to solution using solver compiled in the KenyaCoV module
@time sol_tl_shutdown_after_week = solve_KenyaCoV_prob(u0,(0.,360.),P,0.25)
