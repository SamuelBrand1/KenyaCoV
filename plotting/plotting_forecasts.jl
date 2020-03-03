push!(LOAD_PATH, joinpath(homedir(),"GitHub/KenyaCoV/src"))
using Plots,Parameters,Distributions,JLD2,DataFrames,StatsPlots,FileIO,MAT
using Statistics: median, quantile

treatment_rates = [0.,1/7,1/3.5,1]
include("plotting_functions.jl");

@load joinpath(homedir(),"Github/KenyaCoVOutputs/results_A.jld2") results_A
plt = plot_total_incidence(results_A,treatment_rates)
savefig(plt,"plotting/early_growth_A.png")


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

savefig(plt_incidence_by_age,"plotting/total_cases_by_age_and_treatment.png")
plt = plot_total_incidence_by_treatment(results_A,treatment_rates)
savefig(plt,"plotting/total_cases_by_treatment.png")
