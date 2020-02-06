
#Test --- if ϵ (infectiousness of asymptomatics is zero) and external rates are zero then nothing should happen despite 1000 initial asymptomatics

function test_no_infecteds()
    #Load data completely susceptible Population
    u0,P,transport_matrix = KenyaCoV.model_ingredients_from_data("data/combined_population_estimates.csv","data/optimal_transition_matrix.jld2","data/optimal_movement_matrix.jld2" )
    """
    Example of methods that modify underlying parameters
    """
    #This method modifies the parameter set for changing the mixing structure
    # KenyaCoV.transportstructure_params!(P,[0.001 for i = 1:KenyaCoV.n],transport_matrix)
    #You can modify other parameters directly
    P.τ = 0. #e.g. no treatment rate
    P.ϵ = 0.
    P.ext_mom = 0.
    P.ext_nai = 0.
    #Define initial conditions by modifying the completely susceptible population
    u0[30,3,1] += 1000#One asymptomatic in Nairobi

    #Create a JumpProblem which you can solve --- needs DifferentialEquations module for the solvers
    jump_prob_tl = create_KenyaCoV_prob(u0,(0.,365.),P)
    #Go straight to solution using solver compiled in the KenyaCoV module
    @time sol_tl = solve_KenyaCoV_prob(u0,(0.,365.),P,1.)

    return ~any([sum(reshape(sol_tl[end],n,n_s,2)[i,7:8,1:2]) for i = 1:47] .!= 0)
end
