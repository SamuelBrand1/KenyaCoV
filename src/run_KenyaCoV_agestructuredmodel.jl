push!(LOAD_PATH, joinpath(homedir(),"GitHub/KenyaCoV/src"))
using Plots,Parameters,Distributions,DifferentialEquations,JLD2,DataFrames
using Revise
import KenyaCoV
using LinearAlgebra:eigen


u0,P,P_dest = KenyaCoV.model_ingredients_from_data("data/data_for_age_structuredmodel.jld2",
                                            "data/flight_numbers.csv",
                                            "data/projected_global_prevelance.csv")
"""
Can adjust β to match a desired R₀ value by evaluating the leading eigenvalue of the age mixing matrix
"""
eigs, = eigen(P.M)
max_eigval = Real(eigs[end])
P.β = 2.2*P.γ/max_eigval
