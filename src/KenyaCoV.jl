__precompile__(true)

module KenyaCoV

using DifferentialEquations,
    DataFrames,
    Parameters,
    LinearAlgebra,
    JLD2,
    Distributions,
    SparseArrays,
    CSV,
    RecursiveArrayTools,
    Dates,
    Plots,
    DelimitedFiles

export model_ingredients_from_data,
        transportstructure_params!,
        create_KenyaCoV_prob,
        solve_KenyaCoV_prob
n = 47 #global defining number of counties
n_s = 12 #global defining number of states
n_t = 16 #global defining number of events per location
n_a = 17 #global defining number of age categories
n_wa = 20 #global defining number of wider area groupings
n_ta = 8 #global defining number of events per location and age group
mobile_age_indices = 5:11; #This assumes that 16-49 year ols move around and others don't
immobile_age_indices = [1,2,3,4,12,13,14,15,16,17]

ind_mombasa = 28
ind_nairobi = 30
ind_mombasa_as = 12
ind_nairobi_as = 4

# index_as = CartesianIndices((1:n_wa, 1:n_a,1:n_s))
# index_as_events = CartesianIndices((1:n_wa, 1:n_a,1:n_ta))
# linear_as = LinearIndices((1:n_wa, 1:n_a,1:n_s))
# linear_as_events = LinearIndices((1:n_wa, 1:n_a,1:n_ta))


index_as = CartesianIndices((1:n, 1:n_a,1:n_s))
index_as_events = CartesianIndices((1:n, 1:n_a,1:n_ta))
linear_as = LinearIndices((1:n, 1:n_a,1:n_s))
linear_as_events = LinearIndices((1:n, 1:n_a,1:n_ta))

@load "data/posterior_distribution_R0.jld2" posterior_R₀


include("types.jl");
include("dynamics.jl")
include("transmissionmodel.jl");
include("hospitalisation_data.jl")
include("hospital_model.jl")
include("forecast_functions.jl");



end # module
