push!(LOAD_PATH, joinpath(homedir(),"GitHub/KenyaCoV/src"))
using Plots,Parameters,Distributions,DifferentialEquations,JLD2,DataFrames
using Revise
import KenyaCoV
using LinearAlgebra:eigen


"""
Load age structured data
"""

u0,P,P_dest = KenyaCoV.model_ingredients_from_data("data/data_for_age_structuredmodel.jld2",
                                            "data/flight_numbers.csv",
                                            "data/projected_global_prevelance.csv")
"""
Can adjust β to match a desired R₀ value by evaluating the leading eigenvalue
"""
sus_matrix = repeat(P.χ,1,16)
eigs, = eigen(sus_matrix.*P.M)
max_eigval = Real(eigs[end])
P.β = 2.2*P.γ/max_eigval
# @time KenyaCoV.rates(P.poi_rates,u0,P,1.)
#dc = P.dc
# @time KenyaCoV.nonneg_tauleap(deepcopy(u0),u0,P,0.)

u0[KenyaCoV.ind_nairobi_as,5,4] += 10#10 diseased

prob = KenyaCoV.create_KenyaCoV_non_neg_prob(u0,(0.,365.),P)
@time sol = solve(prob,FunctionMap(),dt = 0.25)

prob_ode = KenyaCoV.create_KenyaCoV_ode_prob(u0,(0.,365.),P)

@time sol_ode = solve(prob_ode,Tsit5())
