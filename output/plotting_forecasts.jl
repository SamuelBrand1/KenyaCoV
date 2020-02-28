push!(LOAD_PATH, joinpath(homedir(),"GitHub/KenyaCoV/src"))
using Plots,Parameters,Distributions,JLD2,DataFrames,StatsPlots,FileIO,MAT
using Statistics: median, quantile
include("plotting_functions.jl");
@load joinpath(homedir(),"Github/KenyaCoVOutputs/results_A.jld2") results_A

plt_incidence = plot_incidence_timeseries(results_A,treatment_rates)
savefig(plt_incidence,"plotting/daily_incidence_A.png")
plt_incidence_spatial_two_week = plot_incidence_spatial(results_A,treatment_rates,3)
savefig(plt_incidence_spatial_two_week,"plotting/daily_incidence_spatial_A_two_week.png")

plt_incidence_by_age = plot_total_incidence_by_age(results_A,treatment_rates)
savefig(plt_incidence_by_age,"plotting/total_cases_by_age_and_treatment.png")
plt = plot_total_incidence_by_treatment(results_A,treatment_rates)
savefig(plt,"plotting/total_cases_by_treatment.png")
