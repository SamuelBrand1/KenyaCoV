push!(LOAD_PATH, joinpath(homedir(),"GitHub/KenyaCoV/src"))
using Plots,Parameters,Distributions,DifferentialEquations,JLD2,DataFrames,CSV,RecursiveArrayTools
using Revise
import KenyaCoV
using LinearAlgebra:eigen,normalize


"""
Load age structured data
"""

u0,P,P_dest = KenyaCoV.model_ingredients_from_data("data/data_for_age_structuredmodel_with_counties.jld2",
                                            "data/flight_numbers.csv",
                                            "data/projected_global_prevelance.csv")
N = sum(u0[:,:,1])


@load "data/agemixingmatrix_china.jld2" M_China
@load "data/agemixingmatrix_Kenya_norestrictions.jld2" M_Kenya
@load "data/agemixingmatrix_Kenya_homeonly.jld2" M_Kenya_ho
@load "data/detection_rates_for_different_epsilons_model2.jld2" d_0 d_01 d_025 d_05 d_1

counties = CSV.read("data/2019_census_age_pyramids_counties.csv")
names = counties.county
Nairobi_index = findfirst(counties.county .== "Nairobi")
Mombassa_index = findfirst(counties.county .== "Mombasa")
Kwale_index = findfirst(counties.county .== "Kwale")
Kilifi_index = findfirst(counties.county .== "Kilifi")
Mandera_index  = findfirst(counties.county .== "Mandera")

χ_zhang = vcat(0.34*ones(3),ones(10),1.47*ones(4))

P.χ = copy(χ_zhang)
P.rel_detection_rate = d_1
P.dt = 0.25;
P.ext_inf_rate = 0.;
P.ϵ = 1.
#Set the susceptibility vector --- just to specify the correct R₀
sus_matrix = repeat(χ_zhang,1,17)
R_A = P.ϵ*((1/P.σ₂) + (1/P.γ) ) #effective duration of asymptomatic
R_M = (P.ϵ/P.σ₂) + (P.ϵ_D/P.γ) #effective duration of mild
R_V = (P.ϵ/P.σ₂) + (P.ϵ_V/P.τ) #effective duration of severe
R_vector = [(1-P.rel_detection_rate[a])*R_A + P.rel_detection_rate[a]*(1-P.hₐ[a])*R_M + P.rel_detection_rate[a]*P.hₐ[a]*R_V for a = 1:17]
inf_matrix = repeat(R_vector',17,1)

eigs_china, = eigen(sus_matrix.*M_China.*inf_matrix)
max_eigval_china = Real(eigs_china[end])
eigs_kenya, = eigen(sus_matrix.*M_Kenya.*inf_matrix)
max_eigval_Kenya = Real(eigs_kenya[end])
multiplier_for_kenya = max_eigval_Kenya/max_eigval_china
P.χ .= χ_zhang ./max_eigval_china #This rescales everything so β is the same as R₀ for China

# u0[Mombassa_index,8,3] = 10 #10 initial pre-symptomatics in Mombasa
# u0[Mandera_index,8,3] = 5 #5 initial pre-symptomatics in Mandera



P.dt = 0.25



"""
Run model

"""
P.β = 3.
u0[Nairobi_index,8,3] = 30 #10 initial pre-symptomatics in Nairobi

prob = KenyaCoV.create_KenyaCoV_non_neg_prob(u0,(0.,60.),P)
@time sol = solve(prob,FunctionMap(),dt = P.dt)
sims_test = KenyaCoV.run_consensus_simulations(P,prob,10,CallbackSet())
output = extract_information_from_simulations(sims_test);

plt_HU,plt_ICU = plot_ranked_bars_health_usage(output," (test)",names)

model_str =
"""
This is a test String
for model description.
"""
