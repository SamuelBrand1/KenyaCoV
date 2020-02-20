
#Test --- if ϵ (infectiousness of asymptomatics is zero) and external rates are zero then nothing should happen despite 1000 initial asymptomatics

function test_no_infecteds()
    #Load data completely susceptible Population
    u0,P,transport_matrix = KenyaCoV.model_ingredients_from_data("data/combined_population_estimates.csv",
                                                                 "data/optimal_transition_matrix.jld2",
                                                                "data/optimal_movement_matrix.jld2",
                                                                "data/flight_numbers.csv",
                                                                "data/projected_global_prevelance.csv")
    #You can modify other parameters directly
    P.τ = 0. #e.g. no treatment rate
    P.ϵ = 0.
    P.into_mom .= 0
    P.into_nai .= 0
    #Define initial conditions by modifying the completely susceptible population
    u0[30,3,1] += 1000#One asymptomatic in Nairobi

    #Go straight to solution using solver compiled in the KenyaCoV module
    sol_tl = KenyaCoV.solve_KenyaCoV_prob(u0,(0.,365.),P,1.)

    return ~any([sum(reshape(sol_tl[end],n,n_s,2)[i,7:8,1:2]) for i = 1:47] .!= 0)
end

function stays_nonnegative()
    u0,P,transport_matrix = KenyaCoV.model_ingredients_from_data("data/combined_population_estimates.csv",
                                                                 "data/optimal_transition_matrix.jld2",
                                                                 "data/optimal_movement_matrix.jld2",
                                                                 "data/flight_numbers.csv",
                                                                 "data/projected_global_prevelance.csv")
    P.τ = 0. #e.g. no treatment rate
    for (i,p) in enumerate(P.global_prev)
        P.global_prev[i] = 0.
    end
    # P.global_prev
    P.μ₁ = 0.
    P.ϵ = 1.
    P.β = 2.2*P.γ
    #Define initial conditions by modifying the completely susceptible population
    P.dt = 1.
    u0[30,4,1] += 10#One diseased in Nairobi
    prob = KenyaCoV.create_KenyaCoV_non_neg_prob(u0,(0.,365.),P)
    sol = solve(prob,FunctionMap(),dt = P.dt)
    ~any([any(u .< 0) for u in sol.u])
end
# function run_epidemic_infs_everywhere()
#     u0,P,transport_matrix = KenyaCoV.model_ingredients_from_data("data/combined_population_estimates.csv",
#                                                                  "data/optimal_transition_matrix.jld2",
#                                                                 "data/optimal_movement_matrix.jld2",
#                                                                 "data/flight_numbers.csv",
#                                                                 "data/projected_global_prevelance.csv")
#
#     u0[30,4,1] += 100#One diseased in Nairobi
#
#     #Create a JumpProblem which you can solve --- needs DifferentialEquations module for the solvers
#     jump_prob_tl = create_KenyaCoV_prob(u0,(0.,365.),P)
#     #Go straight to solution using solver compiled in the KenyaCoV module
#     sol_tl = solve_KenyaCoV_prob(u0,(0.,365.),P,1.)
# end
