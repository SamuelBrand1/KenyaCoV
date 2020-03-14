push!(LOAD_PATH, joinpath(homedir(),"GitHub/KenyaCoV/src"))
using Plots,Parameters,Distributions,JLD2,DataFrames,StatsPlots,FileIO,MAT
using Statistics: median, quantile

treatment_rates = [(0.,1),(0.,0.5),(1/7.,1.),(1/7,0.5),(1/3.5,1.),(1/3.5,0.5)]

"""
Plots ---
1) Early growth for both scenarios
2) Peaks
3) Size and age distribution
"""
results_1

@load joinpath(homedir(),"Github/KenyaCoVOutputs/results_1.jld2") results_1
@load joinpath(homedir(),"Github/KenyaCoVOutputs/results_2.jld2") results_2
@load joinpath(homedir(),"Github/KenyaCoVOutputs/results_3.jld2") results_3

include("plotting_functions.jl");
res_group = [results_1]

plt = plot_total_incidence(res_group,treatment_rates[2],5)
title!(plt,"Av. 1 week to isolation, 25% reduction in infectiousness")
savefig(plt,"plotting/incidence_with_intervention.png")

plot_total_incidence(results_1,treatment_rates)

"""
Spatial plot
"""
inc_D_at_60 = results_A[1][1][1:20,60,1]
ordering = sortperm(inc_D_at_60,rev = true)

plt_no_intervention = plot_incidence_spatial(results_A,treatment_rates,1,ordering)
savefig(plt_no_intervention,"plotting/rate_of_early_growth_A.png")

plt_incidence_by_age_noi = plot_total_incidence_by_age(results_A,treatment_rates,1)
plt_incidence_by_age_1wk = plot_total_incidence_by_age(results_A,treatment_rates,2)
plt_incidence_by_age_halfwk = plot_total_incidence_by_age(results_A,treatment_rates,3)
plt_incidence_by_age_oneday = plot_total_incidence_by_age(results_A,treatment_rates,4)

savefig(plt_incidence_by_age_noi,"plotting/total_cases_by_age_noi.png")
savefig(plt_incidence_by_age_1wk,"plotting/total_cases_by_age_1wk.png")
savefig(plt_incidence_by_age_halfwk,"plotting/total_cases_by_age_halfwk.png")
savefig(plt_incidence_by_age_oneday,"plotting/total_cases_by_age_oneday.png")
