push!(LOAD_PATH, "/Users/Sam/GitHub/KenyaCoV/src")
using Plots,Parameters,Distributions,DifferentialEquations,KenyaCoV
u0,P,transport_matrix = model_ingredients_from_data("./src/2009_National_Estimates_of_Rural_and_Urban_Populations_by_County.csv",0.01);
KenyaCoV.transportstructure_params!(P,0.001,transport_matrix)
P.τ = 0.
P.β = 2.5*(P.γ)
