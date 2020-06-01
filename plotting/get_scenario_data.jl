
@load("reports/report_unmitigated/scenario_data_unmitigated.jld2");
scenariodata_unmitigated = deepcopy(scenariodata);
@load("reports/report_full_intervention/scenario_data_full_intervention.jld2");
scenariodata_full_intervention = deepcopy(scenariodata);

@load("reports/report_scenario_1a/scenario_data_scenario_1a.jld2");
scenariodata_1a = deepcopy(scenariodata);
@load("reports/report_scenario_1b/scenario_data_scenario_1b.jld2");
scenariodata_1b = deepcopy(scenariodata);
@load("reports/report_scenario_1c/scenario_data_scenario_1c.jld2");
scenariodata_1c = deepcopy(scenariodata);
@load("reports/report_scenario_1d/scenario_data_scenario_1d.jld2");
scenariodata_1d = deepcopy(scenariodata);

@load("reports/report_candidates_june_90perc/scenario_data_candidates_june_90perc.jld2");
scenariodata_candidates_june_90perc = deepcopy(scenariodata);
@load("reports/report_candidates_june_50perc/scenario_data_candidates_june_50perc.jld2");
scenariodata_candidates_june_50perc = deepcopy(scenariodata);

@load("reports/report_candidates_august_90perc/scenario_data_candidates_august_90perc.jld2");
scenariodata_candidates_august_90perc = deepcopy(scenariodata);

@load("reports/report_candidates_august_50perc/scenario_data_candidates_august_50perc.jld2");
scenariodata_candidates_august_50perc = deepcopy(scenariodata);

@load("reports/report_primary_june_50perc/scenario_data_primary_june_50perc.jld2");
scenariodata_primary_june_50perc = deepcopy(scenariodata);
@load("reports/report_secondary_june_50perc/scenario_data_secondary_june_50perc.jld2");
scenariodata_secondary_june_50perc = deepcopy(scenariodata);

@load("reports/report_tertiary_june_50perc/scenario_data_tertiary_june_50perc.jld2");
scenariodata_tertiary_june_50perc = deepcopy(scenariodata);

@load("reports/report_primary_june_90perc/scenario_data_primary_june_90perc.jld2");
scenariodata_primary_june_90perc = deepcopy(scenariodata);
@load("reports/report_secondary_june_90perc/scenario_data_secondary_june_90perc.jld2");
scenariodata_secondary_june_90perc = deepcopy(scenariodata);
@load("reports/report_tertiary_june_90perc/scenario_data_tertiary_june_90perc.jld2");
scenariodata_tertiary_june_90perc = deepcopy(scenariodata);


##



function generate_death_report(output,simulation_tag)


    writedlm("reports/country_deaths_ts"*simulation_tag*".csv",output.country_incidence_death_ts.med,",")
    writedlm("reports/country_cum_deaths_ts"*simulation_tag*".csv",cumsum(output.country_incidence_death_ts.med),",")

    writedlm("reports/death_incidence_ts_by_county"*simulation_tag*".csv",hcat(areanames,output.incidence_death_ts.med),",")
    writedlm("reports/cum_death_incidence_ts_by_county"*simulation_tag*".csv",hcat(areanames,cumsum(output.incidence_death_ts.med,dims=2)),",")

    return scenariodata
end
