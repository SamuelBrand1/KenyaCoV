using JLD2,Plots,StatsPlots,CSV,DelimitedFiles
using Plots.PlotMeasures
include("get_scenario_data.jl");
## Plot files
#generate_death_report(scenariodata_tertiary_june_90perc,"_tertiary_june_90perc");



##


scenarionames = ["Unmitigated",
                "Schools reopen June",
                "Schools reopen August",
                "Candidates only return June",
                "Candidates only return August",
                "Primary age only return June",
                "Secondary age only return June",
                "Tertiary age only return June",
                "Schools remain shut"]

## Deaths

median_deaths_by_scenario_90perc = [scenariodata_unmitigated.total_deaths.med,
                            scenariodata_1a.total_deaths.med,
                            scenariodata_1c.total_deaths.med,
                            scenariodata_candidates_june_50perc.total_deaths.med,
                            scenariodata_candidates_august_50perc.total_deaths.med,
                            scenariodata_primary_june_50perc.total_deaths.med,
                            scenariodata_secondary_june_50perc.total_deaths.med,
                            scenariodata_tertiary_june_50perc.total_deaths.med,
                            scenariodata_full_intervention.total_deaths.med]

lpred_deaths_by_scenario_90perc = [scenariodata_unmitigated.total_deaths.lpred,
                            scenariodata_1a.total_deaths.lpred,
                            scenariodata_1c.total_deaths.lpred,
                            scenariodata_candidates_june_50perc.total_deaths.lpred,
                            scenariodata_candidates_august_50perc.total_deaths.lpred,
                            scenariodata_primary_june_50perc.total_deaths.lpred,
                            scenariodata_secondary_june_50perc.total_deaths.lpred,
                            scenariodata_tertiary_june_50perc.total_deaths.lpred,
                            scenariodata_full_intervention.total_deaths.lpred]
upred_deaths_by_scenario_90perc = [scenariodata_unmitigated.total_deaths.upred,
                            scenariodata_1a.total_deaths.upred,
                            scenariodata_1c.total_deaths.upred,
                            scenariodata_candidates_june_50perc.total_deaths.upred,
                            scenariodata_candidates_august_50perc.total_deaths.upred,
                            scenariodata_primary_june_50perc.total_deaths.upred,
                            scenariodata_secondary_june_50perc.total_deaths.upred,
                            scenariodata_tertiary_june_50perc.total_deaths.upred,
                            scenariodata_full_intervention.total_deaths.upred]

I = sortperm(median_deaths_by_scenario_90perc,rev = true)

plt_deaths_90 = bar(median_deaths_by_scenario_90perc[I],orientations = :horizonal,
            yticks = (1:length(I),scenarionames[I]),size = (800,550),lab ="",
            xlabel = "Deaths",title = "Deaths by scenario - 90% schools contacts after return",left_margin = 10mm,bottom_margin = 5mm)
scatter!(plt_deaths_90,median_deaths_by_scenario_90perc[I],1:length(I),
            xerror = (median_deaths_by_scenario_90perc[I] .- lpred_deaths_by_scenario_90perc[I],
                        upred_deaths_by_scenario_90perc[I] .- median_deaths_by_scenario_90perc[I]),
            ms = 0.,color = :black,lab ="")
savefig(plt,"example_plot_scenarios.png")
