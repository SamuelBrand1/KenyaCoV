__precompile__(true)

module KenyaCoV_contacts

using DifferentialEquations,
    DataFrames,
    Parameters,
    LinearAlgebra,
    JLD2,
    Distributions,
    SparseArrays,
    Random                                                                      #****

export n,n_t,n_s,n_a,
        model_ingredients_from_data,
        transportstructure_params!,
        create_KenyaCoV_prob,
        solve_KenyaCoV_prob, L_IQ
n = 47 #global defining number of counties
#n_s = 9 #global defining number of state
n_s = 14 #global defining number of state                                       #****2      #**** Added one state IQ
n_t = 16 #global defining number of events per location
n_a = 16 #global defining number of age categories
n_wa = 20 #global defining number of wider area groupings
#n_ta = 8 #global defining number of events per location and age group
n_ta = 16 #global defining number of events per location and age group          #****2      #**** Added two events E->IQ and IQ->R  + modeified Iᴰ->H  to IQ->H
mobile_age_indices = 5:11; #This assumes that 16-49 year olds move around and others don't
immobile_age_indices = [1,2,3,4,12,13,14,15,16]

ind_mombasa = 28
ind_nairobi = 30
ind_mombasa_as = 12
ind_nairobi_as = 4
#L_IQ=[]
index_as = CartesianIndices((1:n_wa, 1:n_a,1:n_s))
index_as_events = CartesianIndices((1:n_wa, 1:n_a,1:n_ta))
linear_as = LinearIndices((1:n_wa, 1:n_a,1:n_s))
linear_as_events = LinearIndices((1:n_wa, 1:n_a,1:n_ta))


include("../src/kenya_data.jl");                                                #****
include("../src/gravity_model.jl");
include("contacttracing_types.jl");
include("../src/regularjumps.jl");
include("contacttracing_jumps_withTmax.jl")
include("../src/transmissionmodel.jl");


end # module