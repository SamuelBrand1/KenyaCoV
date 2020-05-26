using JLD2,Plots,StatsPlots,CSV,DelimitedFiles
using Plots.PlotMeasures
include("get_scenario_data.jl");
## Plot files
generate_death_report(scenariodata_tertiary_june_90perc,"_tertiary_june_90perc");



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

## Deaths 50% contacts on return
# scenariodata_1a.model_str
median_deaths_by_scenario_50perc = [scenariodata_unmitigated.total_deaths.med,
                            scenariodata_1a.total_deaths.med,
                            scenariodata_1c.total_deaths.med,
                            scenariodata_candidates_june_50perc.total_deaths.med,
                            scenariodata_candidates_august_50perc.total_deaths.med,
                            scenariodata_primary_june_50perc.total_deaths.med,
                            scenariodata_secondary_june_50perc.total_deaths.med,
                            scenariodata_tertiary_june_50perc.total_deaths.med,
                            scenariodata_full_intervention.total_deaths.med]

lpred_deaths_by_scenario_50perc = [scenariodata_unmitigated.total_deaths.lpred,
                            scenariodata_1a.total_deaths.lpred,
                            scenariodata_1c.total_deaths.lpred,
                            scenariodata_candidates_june_50perc.total_deaths.lpred,
                            scenariodata_candidates_august_50perc.total_deaths.lpred,
                            scenariodata_primary_june_50perc.total_deaths.lpred,
                            scenariodata_secondary_june_50perc.total_deaths.lpred,
                            scenariodata_tertiary_june_50perc.total_deaths.lpred,
                            scenariodata_full_intervention.total_deaths.lpred]
upred_deaths_by_scenario_50perc = [sum(scenariodata_unmitigated.deaths_by_age.upred),
                            sum(scenariodata_1a.deaths_by_age.upred),
                            sum(scenariodata_1c.deaths_by_age.upred),
                            sum(scenariodata_candidates_june_50perc.deaths_by_age.upred),
                            sum(scenariodata_candidates_august_50perc.deaths_by_age.upred),
                            sum(scenariodata_primary_june_50perc.deaths_by_age.upred),
                            sum(scenariodata_secondary_june_50perc.deaths_by_age.upred),
                            sum(scenariodata_tertiary_june_50perc.deaths_by_age.upred),
                            sum(scenariodata_full_intervention.deaths_by_age.upred)]

I = sortperm(median_deaths_by_scenario_50perc,rev = true)

plt_deaths_50 = bar(median_deaths_by_scenario_50perc[I],orientations = :horizonal,
            yticks = (1:length(I),scenarionames[I]),size = (800,550),lab ="",
            xlabel = "Deaths (thousands)",
            title = "Deaths by scenario - 50% schools contacts after return",
            left_margin = 10mm,bottom_margin = 5mm,
            xticks = ([0,1e4,2e4,3e4,4e4,5e4],[0,10,20,30,40,50]))
scatter!(plt_deaths_50,median_deaths_by_scenario_50perc[I],1:length(I),
            xerror = (median_deaths_by_scenario_50perc[I] .- lpred_deaths_by_scenario_50perc[I],
                        upred_deaths_by_scenario_50perc[I] .- median_deaths_by_scenario_50perc[I]),
            ms = 0.,color = :black,lab ="",lw=3)
savefig(plt_deaths_50,"plotting/deaths_scenarios_50perc.png")

## Deaths 90%

median_deaths_by_scenario_90perc = [scenariodata_unmitigated.total_deaths.med,
                            scenariodata_1b.total_deaths.med,
                            scenariodata_1d.total_deaths.med,
                            scenariodata_candidates_june_90perc.total_deaths.med,
                            scenariodata_candidates_august_90perc.total_deaths.med,
                            scenariodata_primary_june_90perc.total_deaths.med,
                            scenariodata_secondary_june_90perc.total_deaths.med,
                            scenariodata_tertiary_june_90perc.total_deaths.med,
                            scenariodata_full_intervention.total_deaths.med]

lpred_deaths_by_scenario_90perc = [scenariodata_unmitigated.total_deaths.lpred,
                            scenariodata_1b.total_deaths.lpred,
                            scenariodata_1d.total_deaths.lpred,
                            scenariodata_candidates_june_90perc.total_deaths.lpred,
                            scenariodata_candidates_august_90perc.total_deaths.lpred,
                            scenariodata_primary_june_90perc.total_deaths.lpred,
                            scenariodata_secondary_june_90perc.total_deaths.lpred,
                            scenariodata_tertiary_june_90perc.total_deaths.lpred,
                            scenariodata_full_intervention.total_deaths.lpred]

upred_deaths_by_scenario_90perc = [sum(scenariodata_unmitigated.deaths_by_age.upred),
                            sum(scenariodata_1b.deaths_by_age.upred),
                            sum(scenariodata_1d.deaths_by_age.upred),
                            sum(scenariodata_candidates_june_90perc.deaths_by_age.upred),
                            sum(scenariodata_candidates_august_90perc.deaths_by_age.upred),
                            sum(scenariodata_primary_june_90perc.deaths_by_age.upred),
                            sum(scenariodata_secondary_june_90perc.deaths_by_age.upred),
                            sum(scenariodata_tertiary_june_90perc.deaths_by_age.upred),
                            sum(scenariodata_full_intervention.deaths_by_age.upred)]

I = sortperm(median_deaths_by_scenario_90perc,rev = true)

plt_deaths_90 = bar(median_deaths_by_scenario_90perc[I],orientations = :horizonal,
            yticks = (1:length(I),scenarionames[I]),size = (800,550),lab ="",
            xlabel = "Deaths (thousands)",
            title = "Deaths by scenario - 90% schools contacts after return",
            left_margin = 10mm,bottom_margin = 5mm,
            xticks = ([0,1e4,2e4,3e4,4e4,5e4],[0,10,20,30,40,50]))
scatter!(plt_deaths_90,median_deaths_by_scenario_90perc[I],1:length(I),
            xerror = (median_deaths_by_scenario_90perc[I] .- lpred_deaths_by_scenario_90perc[I],
                        upred_deaths_by_scenario_90perc[I] .- median_deaths_by_scenario_90perc[I]),
            ms = 0.,color = :black,lab ="",lw=3)
savefig(plt_deaths_90,"plotting/deaths_scenarios_90perc.png")

## Hospitalised 50%

median_hosp_by_scenario_50perc = [scenariodata_unmitigated.total_severe_cases.med,
                            scenariodata_1a.total_severe_cases.med,
                            scenariodata_1c.total_severe_cases.med,
                            scenariodata_candidates_june_50perc.total_severe_cases.med,
                            scenariodata_candidates_august_50perc.total_severe_cases.med,
                            scenariodata_primary_june_50perc.total_severe_cases.med,
                            scenariodata_secondary_june_50perc.total_severe_cases.med,
                            scenariodata_tertiary_june_50perc.total_severe_cases.med,
                            scenariodata_full_intervention.total_severe_cases.med]

lpred_hosp_by_scenario_50perc = [scenariodata_unmitigated.total_severe_cases.lpred,
                            scenariodata_1a.total_severe_cases.lpred,
                            scenariodata_1c.total_severe_cases.lpred,
                            scenariodata_candidates_june_50perc.total_severe_cases.lpred,
                            scenariodata_candidates_august_50perc.total_severe_cases.lpred,
                            scenariodata_primary_june_50perc.total_severe_cases.lpred,
                            scenariodata_secondary_june_50perc.total_severe_cases.lpred,
                            scenariodata_tertiary_june_50perc.total_severe_cases.lpred,
                            scenariodata_full_intervention.total_severe_cases.lpred]

upred_hosp_by_scenario_50perc = [sum(scenariodata_unmitigated.severe_cases_by_age.upred),
                            sum(scenariodata_1a.severe_cases_by_age.upred),
                            sum(scenariodata_1c.severe_cases_by_age.upred),
                            sum(scenariodata_candidates_june_50perc.severe_cases_by_age.upred),
                            sum(scenariodata_candidates_august_50perc.severe_cases_by_age.upred),
                            sum(scenariodata_primary_june_50perc.severe_cases_by_age.upred),
                            sum(scenariodata_secondary_june_50perc.severe_cases_by_age.upred),
                            sum(scenariodata_tertiary_june_50perc.severe_cases_by_age.upred),
                            sum(scenariodata_full_intervention.severe_cases_by_age.upred)]

I = sortperm(median_hosp_by_scenario_50perc,rev = true)

plt_hosp_50 = bar(median_hosp_by_scenario_50perc[I],orientations = :horizonal,
            yticks = (1:length(I),scenarionames[I]),size = (800,550),lab ="",
            xlabel = "hospitalisations (thousands)",
            title = "Hospitalisations by scenario - 50% schools contacts after return",
            left_margin = 10mm,bottom_margin = 5mm,
            xticks = ([0,0.5e5,1e5,1.5e5,2e5,2.5e5],[0,50,100,150,200,250]))
            # xticks = ([0,1e4,2e4,3e4,4e4,5e4],[0,10,20,30,40,50]))
scatter!(plt_hosp_50,median_hosp_by_scenario_50perc[I],1:length(I),
            xerror = (median_hosp_by_scenario_50perc[I] .- lpred_hosp_by_scenario_50perc[I],
                        upred_hosp_by_scenario_50perc[I] .- median_hosp_by_scenario_50perc[I]),
            ms = 0.,color = :black,lab ="",lw=3)
savefig(plt_hosp_50,"plotting/hosp_scenarios_50perc.png")

## Hospitalised 90%

median_hosp_by_scenario_90perc = [scenariodata_unmitigated.total_severe_cases.med,
                            scenariodata_1b.total_severe_cases.med,
                            scenariodata_1d.total_severe_cases.med,
                            scenariodata_candidates_june_90perc.total_severe_cases.med,
                            scenariodata_candidates_august_90perc.total_severe_cases.med,
                            scenariodata_primary_june_90perc.total_severe_cases.med,
                            scenariodata_secondary_june_90perc.total_severe_cases.med,
                            scenariodata_tertiary_june_90perc.total_severe_cases.med,
                            scenariodata_full_intervention.total_severe_cases.med]

lpred_hosp_by_scenario_90perc = [scenariodata_unmitigated.total_severe_cases.lpred,
                            scenariodata_1b.total_severe_cases.lpred,
                            scenariodata_1d.total_severe_cases.lpred,
                            scenariodata_candidates_june_90perc.total_severe_cases.lpred,
                            scenariodata_candidates_august_90perc.total_severe_cases.lpred,
                            scenariodata_primary_june_90perc.total_severe_cases.lpred,
                            scenariodata_secondary_june_90perc.total_severe_cases.lpred,
                            scenariodata_tertiary_june_90perc.total_severe_cases.lpred,
                            scenariodata_full_intervention.total_severe_cases.lpred]

upred_hosp_by_scenario_90perc = [sum(scenariodata_unmitigated.severe_cases_by_age.upred),
                            sum(scenariodata_1b.severe_cases_by_age.upred),
                            sum(scenariodata_1d.severe_cases_by_age.upred),
                            sum(scenariodata_candidates_june_90perc.severe_cases_by_age.upred),
                            sum(scenariodata_candidates_august_90perc.severe_cases_by_age.upred),
                            sum(scenariodata_primary_june_90perc.severe_cases_by_age.upred),
                            sum(scenariodata_secondary_june_90perc.severe_cases_by_age.upred),
                            sum(scenariodata_tertiary_june_90perc.severe_cases_by_age.upred),
                            sum(scenariodata_full_intervention.severe_cases_by_age.upred)]

I = sortperm(median_hosp_by_scenario_90perc,rev = true)

plt_hosp_90 = bar(median_hosp_by_scenario_90perc[I],orientations = :horizonal,
            yticks = (1:length(I),scenarionames[I]),size = (800,550),lab ="",
            xlabel = "hospitalisations (thousands)",
            title = "Hospitalisations by scenario - 90% schools contacts after return",
            left_margin = 10mm,bottom_margin = 5mm,
            xticks = ([0,0.5e5,1e5,1.5e5,2e5,2.5e5],[0,50,100,150,200,250]))
            # xticks = ([0,1e4,2e4,3e4,4e4,5e4],[0,10,20,30,40,50]))
scatter!(plt_hosp_90,median_hosp_by_scenario_90perc[I],1:length(I),
            xerror = (median_hosp_by_scenario_90perc[I] .- lpred_hosp_by_scenario_90perc[I],
                        upred_hosp_by_scenario_90perc[I] .- median_hosp_by_scenario_90perc[I]),
            ms = 0.,color = :black,lab ="",lw=3)
savefig(plt_hosp_90,"plotting/hosp_scenarios_90perc.png")
