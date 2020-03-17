push!(LOAD_PATH, joinpath(homedir(),"GitHub/KenyaCoV/src"))
using Plots,Parameters,Distributions,JLD2,DataFrames,StatsPlots,FileIO,MAT
using Statistics: median, quantile

treatment_rates = [(0.,1),(0.,0.5),(1/7.,1.),(1/7,0.5),(1/3.5,1.),(1/3.5,0.5)]
reducted_treatment_rates = [(0.,1),(1/3.5,0.5)]
reduced_rel_transmission_perc = [10,25]

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
@load joinpath(homedir(),"Github/KenyaCoVOutputs/results_2S.jld2") results_2S
@load joinpath(homedir(),"Github/KenyaCoVOutputs/results_3S.jld2") results_3S
@load joinpath(homedir(),"Github/KenyaCoVOutputs/results_2SI.jld2") results_2SI
@load joinpath(homedir(),"Github/KenyaCoVOutputs/results_3SI.jld2") results_3SI




# gr()
include("plotting_functions.jl");
"""
Country wide incidence
"""
scenario_group = [results_1,results_2,results_3,results_4,results_5]
plt_no_control = plot_total_incidence_group(scenario_group,treatment_rates,1,rel_transmission_perc)
title!(plt_no_control,"Kenya-wide incidence: no intervention")
plot!(plt_no_control,size = (700,500))
xlims!(0.,30.)
savefig(plt_no_control,"plotting/baseline_scenarios.pdf")


plt_controls = plot_total_incidence_group(scenario_group,treatment_rates,6,rel_transmission_perc)
title!(plt_controls,"Kenya-wide incidence: case isolation")
plot!(plt_controls,size = (700,400))
savefig(plt_controls,"plotting/target_scenarios.pdf")

scenario_group_SD = [results_2S,results_3S]
plt_no_control_SD = plot_total_incidence_group(scenario_group_SD,reducted_treatment_rates,1,reduced_rel_transmission_perc)
title!(plt_no_control_SD,"Social distancing only")
plot!(plt_no_control_SD,size = (700,400))
plot!(plt_no_control_SD,legend = :topleft)
# xlims!(0.,30.)
savefig(plt_no_control_SD,"plotting/baseline_scenarios_with_social_distancing.pdf")


scenario_group_SD = [results_2S,results_3S]
plt_control_SD = plot_total_incidence_group(scenario_group_SD,reducted_treatment_rates,2,reduced_rel_transmission_perc)
title!(plt_control_SD,"Social distancing and direct isolation")
plot!(plt_control_SD,size = (700,400))
plot!(plt_control_SD,legend = :topleft)
# xlims!(0.,30.)
savefig(plt_control_SD,"plotting/direct_control_scenarios_with_social_distancing.pdf")

scenario_group_SDI = [results_2SI,results_3SI]
plt_no_control_SD = plot_total_incidence_group(scenario_group_SDI,reducted_treatment_rates,1,reduced_rel_transmission_perc)
title!(plt_no_control_SD,"Social distancing only")
plot!(plt_no_control_SD,size = (700,400))
plot!(plt_no_control_SD,legend = :topright)
# xlims!(0.,30.)
savefig(plt_no_control_SD,"plotting/baseline_scenarios_with_social_distancing.pdf")

plt_control_SD = plot_total_incidence_group(scenario_group_SDI,reducted_treatment_rates,2,reduced_rel_transmission_perc)
title!(plt_control_SD,"Social distancing and case isolation")
plot!(plt_control_SD,size = (700,400))
plot!(plt_control_SD,legend = :topright)
# xlims!(0.,30.)
savefig(plt_control_SD,"plotting/baseline_scenarios_with_social_distancing_and_case_isolation.pdf")

"""
Spatial plots
"""


plt_no_control_1 = plot_incidence_spatial(results_1,treatment_rates,1)
title!(plt_no_control_1,"Uncontrolled epidemic: 0% rel. trans. undetecteds ")
plot!(plt_no_control_1,size = (700,500))
savefig(plt_no_control_1,"plotting/spatial_baseline_scenarios_tau_0.pdf")

plt_control_1 = plot_incidence_spatial(results_1,treatment_rates,6)
title!(plt_control_1,"Targetted intervention: 0% rel. trans. undetecteds ")
plot!(plt_control_1,size = (700,500))
savefig(plt_control_1,"plotting/spatial_targetting_scenarios_tau_0.pdf")

plt_no_control_2 = plot_incidence_spatial(results_2,treatment_rates,1)
title!(plt_no_control_2,"Spatial distribution of incidence: (10% rel. infect.)")
plot!(plt_no_control_2,size = (700,500))
xlims!(0.,60.)
savefig(plt_no_control_2,"plotting/spatial_baseline_scenarios_tau_01.pdf")

plt_control_2 = plot_incidence_spatial(results_2,treatment_rates,6)
title!(plt_control_2,"Targetted intervention: 10% rel. trans. undetecteds ")
plot!(plt_control_2,size = (700,500))
savefig(plt_control_2,"plotting/spatial_targetting_scenarios_tau_01.pdf")


plt_control_SD_2 = plot_incidence_spatial(results_2S,treatment_rates,1)
xlims!(plt_control_SD_2,(0.,150))
title!(plt_control_SD_2,"Social distancing: 10% rel. trans. undetecteds ")
plot!(plt_control_SD_2,size = (700,500))
savefig(plt_control_SD_2,"plotting/spatial_targetting_scenarios_tau_01_social_distancing.pdf")

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
plot!(size = (700,400))
savefig(peak_plt_no_control,"plotting/time_to_peak_uncontrolled.pdf")


collect_control_peak_timing = hcat(results_1[6][3][:,21],
                                      results_2[6][3][:,21],
                                      results_3[6][3][:,21],
                                      results_4[6][3][:,21],
                                      results_5[6][3][:,21])

peak_plt_control = boxplot(collect_control_peak_timing,lab="",
                            # lab = ["rel. infect. undetecteds = 0%" "rel. infect. undetecteds = 10%" "rel. infect. undetecteds = 25%" "rel. infect. undetecteds = 50%" "rel. infect. undetecteds = 100%"],
                            xticks = (1:5,["0%","10%","25%","50%","100%"]))
title!(peak_plt_control,"Epidemic peak timing with case isolation")
ylabel!(peak_plt_control,"Time to peak in days")
xlabel!(peak_plt_control,"Rel. infectiousness of undetected cases")
ylims!((0.,300.))
plot!(size = (700,400))
savefig(peak_plt_control,"plotting/time_to_peak_controlled.pdf")


collect_SDI_peak_timing = hcat(results_2SI[1][3][:,21],
                                      results_3SI[1][3][:,21],
                                      results_2SI[2][3][:,21],
                                      results_3SI[2][3][:,21])

peak_plt_SD = boxplot(collect_SDI_peak_timing,lab="",
                            # lab = ["rel. infect. undetecteds = 0%" "rel. infect. undetecteds = 10%" "rel. infect. undetecteds = 25%" "rel. infect. undetecteds = 50%" "rel. infect. undetecteds = 100%"],
                            xticks = (1:4,["10% - SD","25% - SD","10% - SD+CI","25% - SD+CI"]))
title!(peak_plt_SD,"Epidemic peak timing with social distancing and case isolation")
ylabel!(peak_plt_SD,"Time to peak in days")
xlabel!(peak_plt_SD,"Rel. infectiousness of undetected cases")
plot!(size = (700,400))
savefig(peak_plt_SD,"plotting/time_to_peak_social_distancing_and_case_isolation.pdf")

"""
Final size plots
"""

sum(results_1[1][4],dims= [1,2])[:]

collect_uncontrol_final_size = hcat(sum(results_1[1][4],dims= [1,2])[:],
                                    sum(results_2[1][4],dims= [1,2])[:],
                                    sum(results_3[1][4],dims= [1,2])[:],
                                    sum(results_4[1][4],dims= [1,2])[:],
                                    sum(results_5[1][4],dims= [1,2])[:])

size_plt_no_control = boxplot(collect_uncontrol_final_size./1e6,lab="",
                            # lab = ["rel. infect. undetecteds = 0%" "rel. infect. undetecteds = 10%" "rel. infect. undetecteds = 25%" "rel. infect. undetecteds = 50%" "rel. infect. undetecteds = 100%"],
                            xticks = (1:5,["0%","10%","25%","50%","100%"]))
title!(size_plt_no_control,"Uncontrolled epidemic number of cases")
ylabel!(size_plt_no_control,"Total cases (millions)")
xlabel!(size_plt_no_control,"Rel. infectiousness of undetected cases")
plot!(size = (700,400))
savefig(size_plt_no_control,"plotting/epidemic_size_no_control.pdf")

collect_control_final_size = hcat(sum(results_1[6][4],dims= [1,2])[:],
                                    sum(results_2[6][4],dims= [1,2])[:],
                                    sum(results_3[6][4],dims= [1,2])[:],
                                    sum(results_4[6][4],dims= [1,2])[:],
                                    sum(results_5[6][4],dims= [1,2])[:])

size_plt_control = boxplot(collect_control_final_size./1e6,lab="",
                            # lab = ["rel. infect. undetecteds = 0%" "rel. infect. undetecteds = 10%" "rel. infect. undetecteds = 25%" "rel. infect. undetecteds = 50%" "rel. infect. undetecteds = 100%"],
                            xticks = (1:5,["0%","10%","25%","50%","100%"]))
title!(size_plt_control,"Epidemic number of cases with case isolation")
ylabel!(size_plt_control,"Total cases (millions)")
xlabel!(size_plt_control,"Rel. infectiousness of undetected cases")
ylims!((0.,1.5))
plot!(size = (700,400))
savefig(size_plt_control,"plotting/epidemic_size_with_control.pdf")

collect_SD_size = hcat(sum(results_2S[1][4],dims= [1,2])[:],
                            sum(results_3S[1][4],dims= [1,2])[:],
                            sum(results_2S[2][4],dims= [1,2])[:],
                            sum(results_3S[2][4],dims= [1,2])[:])
size_plt_SD = boxplot(collect_SD_size,lab="",
                            # lab = ["rel. infect. undetecteds = 0%" "rel. infect. undetecteds = 10%" "rel. infect. undetecteds = 25%" "rel. infect. undetecteds = 50%" "rel. infect. undetecteds = 100%"],
                            xticks = (1:4,["10% - no targeting","25% - no targeting","10% - targetting","25% - targetting"]))

collect_SDI_size = hcat(sum(results_2SI[1][4],dims= [1,2])[:],
                            sum(results_3SI[1][4],dims= [1,2])[:],
                            sum(results_2SI[2][4],dims= [1,2])[:],
                            sum(results_3SI[2][4],dims= [1,2])[:])
size_plt_SD = boxplot(collect_SDI_size./1e6,lab="",
                            # lab = ["rel. infect. undetecteds = 0%" "rel. infect. undetecteds = 10%" "rel. infect. undetecteds = 25%" "rel. infect. undetecteds = 50%" "rel. infect. undetecteds = 100%"],
                            xticks = (1:4,["10% - SD","25% - SD","10% - SD+CI","25% - SD+CI"]))
title!(size_plt_SD,"Epidemic number of cases with social distancing and case isolation")
ylabel!(size_plt_SD,"Total cases (millions)")
xlabel!(size_plt_SD,"Rel. infectiousness of undetected cases")
ylims!((0.,1.5))
plot!(size = (700,400))
savefig(size_plt_SD,"plotting/epidemic_size_with_SD_and_case_isolation.pdf")

"""
Age distribution of cases
"""


scenario_group = [results_1,results_2,results_3,results_4,results_5]

plt_age_distrib_uncontrolled = plot_total_incidence_by_age(scenario_group,treatment_rates,rel_transmission_perc,1)
plot!(size = (800,500))
xlabel!("Age of case individual")
ylabel!("Number of cases (millions)")
title!("Age distribution of cases: uncontrolled epidemic")
savefig(plt_age_distrib_uncontrolled,"plotting/age_distrib_cases_uncontrolled.pdf")

plt_age_distrib_controlled = plot_total_incidence_by_age(scenario_group,treatment_rates,rel_transmission_perc,6)
plot!(size = (800,500),legend = :topleft)

xlabel!("Age of case individual")
ylabel!("Number of cases (millions)")
title!("Age distribution of cases: targetted controls")
savefig(plt_age_distrib_controlled,"plotting/age_distrib_cases_targetted.pdf")


scenario_group_SD = [results_1,results_3,results_3SI]
plt_age_distrib_SD = plot_total_incidence_by_age(scenario_group_SD,1,
                    ["0% rel. infectiousness" "25% rel. infectiousness" "25% rel. infectiousness + SD"])
plot!(size = (800,500))
xlabel!("Age of case individual")
ylabel!("Number of cases (thousands)")
title!("Age profile of cases")
savefig(plt_age_distrib_SD,"plotting/age_distrib_cases_selected.pdf")
