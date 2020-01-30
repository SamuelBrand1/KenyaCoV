push!(LOAD_PATH, "/Users/Sam/GitHub/KenyaCoV/src")
using Plots,Parameters,Distributions
using KenyaCoV
u0,P,transport_matrix = model_ingredients_from_data("./src/2009_National_Estimates_of_Rural_and_Urban_Populations_by_County.csv",0.01)
KenyaCoV.transportstructure_params!(P,0.001,transport_matrix)
P.τ = 1/1.
u0[30,3,1] += 9#One asymptomatic in Nairobi
jump_prob_tl = create_KenyaCoV_prob(u0,(0.,365.),P)
@time sol_tl = solve_KenyaCoV_prob(u0,(0.,365.),P)
