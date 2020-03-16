push!(LOAD_PATH, joinpath(homedir(),"GitHub/KenyaCoV/src"))
using Plots,Parameters,Distributions,JLD2,DataFrames,StatsPlots,FileIO,MAT
using Statistics: median, quantile

treatment_rates = [(0.,1),(0.,0.5),(1/7.,1.),(1/7,0.5),(1/3.5,1.),(1/3.5,0.5)]
rel_transmission_perc = [0,10,25,50,100]

"""
Plots ---
1) Early growth for both scenarios
2) Peaks
3) Size and age distribution
"""


@load joinpath(homedir(),"Github/KenyaCoVOutputs/results_1.jld2") results_1
@load joinpath(homedir(),"Github/KenyaCoVOutputs/results_2.jld2") results_2
@load joinpath(homedir(),"Github/KenyaCoVOutputs/results_3.jld2") results_3
@load joinpath(homedir(),"Github/KenyaCoVOutputs/results_4.jld2") results_4
@load joinpath(homedir(),"Github/KenyaCoVOutputs/results_5.jld2") results_5

# gr()
include("plotting_functions.jl");
scenario_group = [results_1,results_2,results_3,results_4,results_5]
plt_no_control = plot_total_incidence_group(scenario_group,treatment_rates,1,rel_transmission_perc)
title!(plt_no_control,"Uncontrolled epidemic: 5 generations of undetected transmission")
plot!(plt_no_control,size = (700,400))
savefig(plt_no_control,"plotting/baseline_scenarios.pdf")


plt_controls = plot_total_incidence_group(scenario_group,treatment_rates,6,rel_transmission_perc)
title!(plt_controls,"Target diseased cases (no social distancing): detection period = 3.5 days")
plot!(plt_controls,size = (700,400))
savefig(plt_controls,"plotting/target_scenarios.pdf")


plt_no_control_1 = plot_incidence_spatial(results_1,treatment_rates,1)
title!(plt_no_control_1,"Uncontrolled epidemic: 0% rel. trans. undetecteds ")
plot!(plt_no_control_1,size = (700,500))
savefig(plt_no_control_1,"plotting/spatial_baseline_scenarios_tau_0.pdf")

plt_control_1 = plot_incidence_spatial(results_1,treatment_rates,6)
title!(plt_control_1,"Targetted intervention: 0% rel. trans. undetecteds ")
plot!(plt_control_1,size = (700,500))
savefig(plt_control_1,"plotting/spatial_targetting_scenarios_tau_0.pdf")

plt_no_control_2 = plot_incidence_spatial(results_2,treatment_rates,1)
title!(plt_no_control_1,"Uncontrolled epidemic: 10% rel. trans. undetecteds ")
plot!(plt_no_control_1,size = (700,500))
savefig(plt_no_control_1,"plotting/spatial_baseline_scenarios_tau_01.pdf")

plt_control_2 = plot_incidence_spatial(results_2,treatment_rates,6)
title!(plt_control_2,"Targetted intervention: 10% rel. trans. undetecteds ")
plot!(plt_control_2,size = (700,500))
savefig(plt_control_2,"plotting/spatial_targetting_scenarios_tau_01.pdf")

"""
PEAK TIMING
"""
collect_no_control_peak_timing = hcat(results_1[1][3][:,21],
                                      results_2[1][3][:,21],
                                      results_3[1][3][:,21],
                                      results_4[1][3][:,21],
                                      results_5[1][3][:,21])

peak_plt_no_control = boxplot(collect_no_control_peak_timing,lab="",
                            # lab = ["rel. infect. undetecteds = 0%" "rel. infect. undetecteds = 10%" "rel. infect. undetecteds = 25%" "rel. infect. undetecteds = 50%" "rel. infect. undetecteds = 100%"],
                            xticks = (1:5,["0%","10%","25%","50%","100%"]))
title!(peak_plt_no_control,"Uncontrolled epidemic peak timing")
ylabel!(peak_plt_no_control,"Time to peak in days")
xlabel!(peak_plt_no_control,"Rel. infectiousness of undetected cases")
plot!(size)
savefig(peak_plt_no_control,"plotting/time_to_peak_uncontrolled.pdf")
