__precompile__(true)

module KenyaCoV

using DifferentialEquations,
    DataFrames,
    Parameters,
    LinearAlgebra,
    JLD2,
    Distributions,
    SparseArrays

export n,n_t,n_s,n_a,
        model_ingredients_from_data,
        transportstructure_params!,
        create_KenyaCoV_prob,
        solve_KenyaCoV_prob
n = 47 #global defining number of counties
n_s = 9 #global defining number of state
n_t = 16 #global defining number of events per location
n_a = 16 #global defining number of age categories
n_wa = 20 #global defining number of wider area groupings

ind_mombasa = 28
ind_nairobi = 30

include("kenya_data.jl");
include("gravity_model.jl");
include("types.jl");
include("regularjumps.jl");
include("agestructurejumps.jl")
include("transmissionmodel.jl");



end # module
