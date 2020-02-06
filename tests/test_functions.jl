
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
    sol_tl = solve_KenyaCoV_prob(u0,(0.,365.),P,1.)

    return ~any([sum(reshape(sol_tl[end],n,n_s,2)[i,7:8,1:2]) for i = 1:47] .!= 0)
end

function stays_nonnegative()
    u0,P,transport_matrix = KenyaCoV.model_ingredients_from_data("data/combined_population_estimates.csv",
                                                                 "data/optimal_transition_matrix.jld2",
                                                                "data/optimal_movement_matrix.jld2",
                                                                "data/flight_numbers.csv",
                                                                "data/projected_global_prevelance.csv")

    u0[30,4,1] += 100#One diseased in Nairobi

    #Create a JumpProblem which you can solve --- needs DifferentialEquations module for the solvers
    jump_prob_tl = create_KenyaCoV_prob(u0,(0.,365.),P)
    #Go straight to solution using solver compiled in the KenyaCoV module
    sol_tl = solve_KenyaCoV_prob(u0,(0.,365.),P,1.)
    return ~any([any(x .< 0) for x in sol_tl.u])
end
u0,P,transport_matrix = KenyaCoV.model_ingredients_from_data("data/combined_population_estimates.csv",
                                                             "data/optimal_transition_matrix.jld2",
                                                            "data/optimal_movement_matrix.jld2",
                                                            "data/flight_numbers.csv",
                                                            "data/projected_global_prevelance.csv")

u0[30,4,1] = 10#One diseased in Nairobi
P.β = 2.5*P.γ
P.μ₁ = 0.
jump_tl = create_KenyaCoV_prob(u0,(0.,365.),P)
sol = solve(jump_tl,SimpleTauLeaping(),dt = 1.)
#Go straight to solution using solver compiled in the KenyaCoV module


[any(x .< 0) for x in sol.u]
x = reshape(sol.u[end],n,n_s,2)
f = findall(x .<0)
sum(x) - sum(u0)
x[f[2]]
