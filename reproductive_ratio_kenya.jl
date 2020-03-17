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
@load "data/detection_rates_for_different_taus.jld2" d_0 d_01 d_025 d_05 d_1


P.χ = ones(KenyaCoV.n_a)
P.rel_detection_rate = d_0
P.dt = 0.25;
P.ext_inf_rate = 0.;
P.ϵ = 0.
P.γ = 1/2.5
P.σ = 1/rand(KenyaCoV.d_incubation)
P.β = 2.5*P.γ
R₀_scale = KenyaCoV.calculate_R₀_scale(P)
P.χ = ones(KenyaCoV.n_a)/R₀_scale
KenyaCoV.calculate_R₀(P) #Before

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
