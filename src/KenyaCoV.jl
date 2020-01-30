__precompile__(true)

module KenyaCoV

using DifferentialEquations,
    DataFrames,
    Parameters,
    LinearAlgebra

export n,n_t,n_s,
        model_ingredients_from_data,
        transportstructure_params!,
        create_KenyaCoV_prob,
        solve_KenyaCoV_prob
n = 47 #global defining number of areas
n_s = 9 #global defining number of state
n_t = 16 #global defining number of events per location

include("kenya_data.jl");
include("gravity_model.jl");
include("types.jl");
include("regularjumps.jl");
include("transmissionmodel.jl");



end # module
