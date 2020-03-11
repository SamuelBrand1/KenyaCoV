push!(LOAD_PATH, joinpath(homedir(),"GitHub/KenyaCoV/src"))
using Plots,Parameters,Distributions,DifferentialEquations,JLD2,DataFrames,StatsPlots,FileIO,MAT,RecursiveArrayTools
import KenyaCoV
using LinearAlgebra:eigen
using Statistics: median, quantile

"""
Load age structured data, define initial state and declare the KenyaCoV problem
"""

u0,P,P_dest = KenyaCoV.model_ingredients_from_data("data/data_for_age_structuredmodel.jld2",
                                            "data/flight_numbers.csv",
                                            "data/projected_global_prevelance.csv")
@load "data/agemixingmatrix_china.jld2" M_China
@load "data/reporting_rate_for_MERS_like_scenario.jld2" reporting_rate

function R_kenya_prediction(R₀,r::Vector{Float64},χ::Vector{Float64},δ,ϵ_A,ϵ_D,τ,γ)
    sus_matrix = repeat(χ,1,KenyaCoV.n_a)
    rep_matrix = repeat(δ*r' + ϵ_A*(1 .- (δ*r')),KenyaCoV.n_a,1)
    rep_matrix_control = repeat(ϵ_D*(γ/(τ+γ))*δ*r' + ϵ_A*(1 .- (δ*r')),KenyaCoV.n_a,1)
    eigs, = eigen(rep_matrix.*sus_matrix.*M_China)
    R₀_scale_china = Real(eigs[end])/γ
    β = R₀/R₀_scale_china
    eigs, = eigen(rep_matrix_control.*sus_matrix.*P.M)
    R₀_scale_kenya = Real(eigs[end])/γ
    return β*R₀_scale_kenya
end
R_kenya_prediction(1,ones(16),ones(16),0.2,1,1.,0.,1/2.5)

x_range = 1:(-0.01):0.5
R_kenya_tau_0 = [R_kenya_prediction(1,ones(16),P.χ,0.2,0.,x,0,1/2.5) for x in x_range]
R_kenya_tau_7 = [R_kenya_prediction(1,ones(16),P.χ,0.2,0.,x,1/7,1/2.5) for x in x_range]
R_kenya_tau_3_5 = [R_kenya_prediction(1,ones(16),P.χ,0.2,0.,x,1/3.5,1/2.5) for x in x_range]


plot((1 .- collect(x_range))*100,R_kenya_tau_0,lw=3,lab = "No isolation")
plot!((1 .- collect(x_range))*100,R_kenya_tau_7,lw=3,lab = "Av. 1 week to isolation")
plot!((1 .- collect(x_range))*100,R_kenya_tau_3_5,lw=3,lab = "Av. 1/2 week to isolation")
plot!([0,50],[1/2.5,1/2.5],lw=3,lab = "Threshold for Chinese reproductive ratio = 2.5",ls= :dash,color=:black)
xlabel!("Percentage reduction in transmissibility of clinical cases")
ylabel!("Relative reproductive ratio for Kenya")
title!("Effectiveness of intervention: SARS-like scenario")
savefig("plotting/eff_SARS_like.png")

e_range = 0:0.01:0.5
MERS_R_tau_0 = [R_kenya_prediction(1,reporting_rate,ones(16),0.2,e,1.,0,1/2.5) for e in e_range]
MERS_R_tau_max = [R_kenya_prediction(1,reporting_rate,ones(16),0.2,e,0.5,1/3.5,1/2.5) for e in e_range]

MERS_R_tau_0_alt = [R_kenya_prediction(1,reporting_rate,ones(16),0.8,e,1.,0,1/2.5) for e in e_range]
MERS_R_tau_max_alt = [R_kenya_prediction(1,reporting_rate,ones(16),0.8,e,0.5,1/3.5,1/2.5) for e in e_range]

plot(e_range,MERS_R_tau_0,lab = "No intervention",lw=3,color = :red,legend = :bottomright)
plot!(e_range,MERS_R_tau_0_alt,lab = "Baseline 80% case detection",lw=3,color=:blue)
plot!(e_range,MERS_R_tau_max,lab = "Max. intervention",lw=3,color = :blue)
plot!(e_range,MERS_R_tau_max_alt,lab = "Baseline 80% case detection - max. intervention",lw=3,color = :blue,ls = :dot)
plot!([0,0.5],[1/2.5,1/2.5],lw=3,lab = "Threshold for Chinese reproductive ratio = 2.5",ls= :dash,color=:black)
xlabel!("Relative infectiousness of undetected cases")
ylabel!("Relative reproductive ratio for Kenya")
title!("Effectiveness of intervention: MERS-like scenario")

savefig("plotting/eff_MERS_like.png")
