# push!(LOAD_PATH, joinpath(homedir(),"/Documents/Covid-19/jl_models/src"))
push!(LOAD_PATH, joinpath(homedir(),"GitHub/KenyaCoV/src"))

using Parameters,Distributions,DifferentialEquations,JLD2,DataFrames,FileIO,DelimitedFiles,RecursiveArrayTools,Plots
using CSV,ExcelFiles,DataFrames
using Statistics: median, quantile
using LinearAlgebra: eigen
using Dates
using StatsPlots

#Load data
# @load joinpath(homedir(),"Documents/Covid-19/jl_models/KenyaCoVOutputs/sims_consensus_new_baseline_vs2.jld2") sims_baseline
#@load joinpath(homedir(),"Documents/Covid-19/jl_models/KenyaCoVOutputs/sims_consensus_end_lockdown.jld2") sims_end_regional_lockdown
#@load joinpath(homedir(),"Documents/Covid-19/jl_models/KenyaCoVOutputs/sims_consensus_open_schools_june.jld2") sims_open_schools_june
# @load joinpath(homedir(),"Documents/Covid-19/jl_models/KenyaCoVOutputs/sims_consensus_open_schools_august.jld2") sims_open_schools_august


@load joinpath(homedir(),"Github/KenyaCoVOutputs/sims_consensus_new_baseline_30_perc_reduction.jld2") sims_baseline
@load joinpath(homedir(),"Github/KenyaCoVOutputs/sims_consensus_end_lockdown_30_perc_reduction.jld2") sims_end_regional_lockdown
@load joinpath(homedir(),"Github/KenyaCoVOutputs/sims_consensus_open_schools_june_30_perc_reduction.jld2") sims_open_schools_june
@load joinpath(homedir(),"Github/KenyaCoVOutputs/sims_consensus_open_schools_august_30_perc_reduction.jld2") sims_open_schools_august



## Functions --- to be used

include("hospitalisations.jl");
counties = CSV.read(joinpath(homedir(),"Github/KenyaCoV/data/2019_census_age_pyramids_counties.csv"))
names = counties.county


# Change activate simulation, titles and tags depending on data loaded

sims = sims_baseline
file_tag = "baseline"
title_tag = " (Baseline)"

# sims = sims_end_regional_lockdown
# file_tag = "end_regional_lockdown_30_perc_SD"
# title_tag = " (30% SD + end regional lockdown May 16th)"

# sims = sims_open_schools_june
# file_tag = "open_schools_june_30_perc_SD"
# title_tag = " (30% SD + opening schools in June)"

sims = sims_open_schools_august
file_tag = "open_schools_august_30_perc_SD"
title_tag = " (30% SD + opening schools in August)"

### Peak calculations by county
@time a,peak_A,peak_A_value = first_introduction_time_peak_peak_value(sims,1)
b,peak_M,peak_M_value = first_introduction_time_peak_peak_value(sims,2)
c,peak_V,peak_V_value = first_introduction_time_peak_peak_value(sims,3)

first_times = get_first_time(a,b,c)
df = print_onset_report(first_times,peak_A,peak_A_value,peak_M,peak_M_value,peak_V,peak_V_value,names,"plotting/peak_report_"*file_tag*".csv")

### Plot the horizontal ranked bars with PIs

plt_cases,plt_deaths = plot_ranked_bars_cases(sims,title_tag)
display(plt_deaths)
savefig(plt_cases,"plotting/cases_"*file_tag*".png")
savefig(plt_deaths,"plotting/deaths_"*file_tag*".png")

a,b,c,d = get_hosp_forecast(sims)# This bit runs the hospital model

#This bit does the plots
plt1,plt2 = plot_ranked_bars_health_usage(a,b,c,d,title_tag)
plt2
savefig(plt1,"plotting/Hospital_capacity_"*file_tag*".png")
savefig(plt2,"plotting/ICU_capacity_"*file_tag*".png")

#This bit plots the dynamics for Nairobi, Mombasa and the rest of the country
# a= [dayofyear(2020,m,1) for m = 4:12] .- dayofyear(2020,3,13)
# b= [dayofyear(2021,m,1) for m = 1:12] .- dayofyear(2021,1,1)
# b = b.+293
first_of_months = vcat([dayofyear(2020,m,1) for m = 4:12] .- dayofyear(2020,3,13),([dayofyear(2021,m,1) for m = 1:12] .- dayofyear(2021,1,1)).+293)
Nairobi_index = findfirst(counties.county .== "Nairobi")
Mombassa_index = findfirst(counties.county .== "Mombasa")
Kwale_index = findfirst(counties.county .== "Kwale")
Kilifi_index = findfirst(counties.county .== "Kilifi")
Mandera_index  = findfirst(counties.county .== "Mandera")

plt_incidence_nairobi,plt_usage_nairobi = give_plots_for_county(sims,Nairobi_index)
plot!(plt_incidence_nairobi,title = "Nairobi"*title_tag)

plot!(plt_usage_nairobi,[0,657],[spare_capacity_H_by_county[Nairobi_index+1],spare_capacity_H_by_county[Nairobi_index+1]],
        title = "Nairobi"*title_tag,lw = 2,ls = :dash,lab = "spare hosp. capacity",color = :blue)

plot!(plt_usage_nairobi,[0,657],[spare_capacity_ICU_by_county[Nairobi_index+1],spare_capacity_ICU_by_county[Nairobi_index+1]],
        title = "Nairobi"*title_tag,lw = 2,ls = :dash,lab = "spare ICU capacity",color = :green)

savefig(plt_incidence_nairobi,"plotting/nairobi_incidence_"*file_tag*".png")
savefig(plt_usage_nairobi,"plotting/nairobi_health_usage_"*file_tag*".png")


plt_incidence_mombasa,plt_usage_mombasa = give_plots_for_county(sims,Mombassa_index)

plot!(plt_incidence_mombasa,title = "Mombasa"*title_tag)

plot!(plt_usage_mombasa,[0,657],[spare_capacity_H_by_county[Mombassa_index+1],spare_capacity_H_by_county[Mombassa_index+1]],
        title = "Mombasa"*title_tag,lw = 2,ls = :dash,lab = "spare hosp. capacity",color = :blue)

plot!(plt_usage_mombasa,[0,657],[spare_capacity_ICU_by_county[Mombassa_index+1],spare_capacity_ICU_by_county[Mombassa_index+1]],
        title = "Mombasa"*title_tag,lw = 2,ls = :dash,lab = "spare ICU capacity",color = :green)

savefig(plt_incidence_mombasa,"plotting/mombasa_incidence_"*file_tag*".png")
savefig(plt_usage_mombasa,"plotting/mombasa_health_usage_"*file_tag*".png")

plt_incidence_rest,plt_usage_rest = give_plots_for_county(sims,setdiff(1:47,[Mombassa_index,Nairobi_index]))
plot!(plt_incidence_rest,title = "Rest of Kenya"*title_tag)
plot!(plt_usage_rest,title = "Rest of Kenya"*title_tag)
savefig(plt_incidence_rest,"plotting/rest_of_country_incidence_"*file_tag*".png")
savefig(plt_usage_rest,"plotting/rest_of_country_health_usage_"*file_tag*".png")

plt_incidence_whole_country,plt_usage_whole_country = give_plots_for_county(sims,1:47)
plot!(plt_incidence_whole_country,title = "Whole country"*title_tag)
plot!(plt_usage_whole_country,title = "Whole country"*title_tag)
savefig(plt_incidence_whole_country,"plotting/whole_country_incidence_"*file_tag*".png")
savefig(plt_usage_whole_country,"plotting/whole_country_health_usage_"*file_tag*".png")


## Total cases by each SCENARIO
function data_summary(sims,nai_index,mombasa_index)
    median_whole_country= median([sum(sims.u[k][end][1:47,:,3]) for k = 1:1000])
    lb_whole_country = quantile([sum(sims.u[k][end][1:47,:,3]) for k = 1:1000],0.025)
    ub_whole_country = quantile([sum(sims.u[k][end][1:47,:,3]) for k = 1:1000],0.975)

    median_nairobi= median([sum(sims.u[k][end][nai_index,:,3]) for k = 1:1000])
    lb_nairobi = quantile([sum(sims.u[k][end][nai_index,:,3]) for k = 1:1000],0.025)
    ub_nairobi= quantile([sum(sims.u[k][end][nai_index,:,3]) for k = 1:1000],0.975)

    median_mombasa= median([sum(sims.u[k][end][mombasa_index,:,3]) for k = 1:1000])
    lb_mombasa = quantile([sum(sims.u[k][end][mombasa_index,:,3]) for k = 1:1000],0.025)
    ub_mombasa= quantile([sum(sims.u[k][end][mombasa_index,:,3]) for k = 1:1000],0.975)

    median_rest_country= median([sum(sims.u[k][end][setdiff(1:47,[nai_index,mombasa_index]),:,3]) for k = 1:1000])
    lb_rest_country = quantile([sum(sims.u[k][end][setdiff(1:47,[Mombassa_index,Nairobi_index]),:,3]) for k = 1:1000],0.025)
    ub_rest_country = quantile([sum(sims.u[k][end][setdiff(1:47,[Mombassa_index,Nairobi_index]),:,3]) for k = 1:1000],0.975)

    dt = [median_whole_country median_nairobi median_mombasa median_rest_country; lb_whole_country lb_nairobi lb_mombasa lb_rest_country; ub_whole_country ub_nairobi ub_mombasa ub_rest_country]

    return dt

end



@load joinpath(homedir(),"Github/KenyaCoVOutputs/sims_consensus_new_baseline_30_perc_reduction.jld2") sims_baseline
dt_baseline = data_summary(sims_baseline,Nairobi_index,Mombassa_index)
sims_baseline = nothing
GC.gc()

@load joinpath(homedir(),"Github/KenyaCoVOutputs/sims_consensus_end_lockdown_30_perc_reduction.jld2") sims_end_regional_lockdown
dt_end_regional_lockdown = data_summary(sims_end_regional_lockdown,Nairobi_index,Mombassa_index)
sims_end_regional_lockdown = nothing
GC.gc()

@load joinpath(homedir(),"Github/KenyaCoVOutputs/sims_consensus_open_schools_june_30_perc_reduction.jld2") sims_open_schools_june
dt_open_schools_june = data_summary(sims_open_schools_june,Nairobi_index,Mombassa_index)
sims_open_schools_june = nothing
GC.gc()

@load joinpath(homedir(),"Github/KenyaCoVOutputs/sims_consensus_open_schools_august_30_perc_reduction.jld2") sims_open_schools_august
dt_open_schools_august = data_summary(sims_open_schools_august,Nairobi_index,Mombassa_index)
sims_open_schools_august = nothing
GC.gc()


function plot_bar_total_cases_by_scenario(region,title_tag)
    meds_vec = [dt_baseline[1,region],dt_end_regional_lockdown[1,region],dt_open_schools_june[1,region],dt_open_schools_august[1,region]]
    lb_vec =   [dt_baseline[2,region],dt_end_regional_lockdown[2,region],dt_open_schools_june[2,region],dt_open_schools_august[2,region]]
    ub_vec =   [dt_baseline[3,region],dt_end_regional_lockdown[3,region],dt_open_schools_june[3,region],dt_open_schools_august[3,region]]

    plt = bar(meds_vec,
            xticks = (1:4,["Baseline: Full int." "End lockdown May 16th" "Open schools June" "Open schools August"]),
            title = "Total severe cases up to end of 2021 ("*title_tag*")",
            ylabel = "Total severe cases",
            xlabel = "Scenario",
            lab="")

    scatter!(plt,(1:4,meds_vec), yerror = (meds_vec.-lb_vec,ub_vec.-meds_vec),ms = 1.,color = :black,lab ="")

    return plt

end


plt_severe_cases_whole_country = plot_bar_total_cases_by_scenario(1,"Whole country")
savefig(plt_severe_cases_whole_country,"plotting/whole_country_severe_cases_by_scenario.png")

#plt_severe_cases_nairobi = plot_bar_total_cases_by_scenario(2,"Nairobi")
#savefig(plt_severe_cases_nairobi,"plotting/nairobi_severe_cases_by_scenario.png")

#plt_severe_cases_mombasa = plot_bar_total_cases_by_scenario(3,"Mombasa")
#savefig(plt_severe_cases_mombasa,"plotting/mombasa_severe_cases_by_scenario.png")

plt_severe_cases_rest = plot_bar_total_cases_by_scenario(4,"Rest of the country")
savefig(plt_severe_cases_rest,"plotting/rest_of_country_severe_cases_by_scenario.png")
