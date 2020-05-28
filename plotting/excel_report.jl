using JLD2,Plots,StatsPlots,CSV,DelimitedFiles, Dates
include("get_scenario_data.jl");




return (total_severe_cases=total_severe_cases,
        severe_cases_by_area=severe_cases_by_area,
        severe_cases_by_age=severe_cases_by_age,
        total_deaths=total_deaths,
        deaths_by_area=deaths_by_area,
        deaths_by_age=deaths_by_age,
        hosp_peak_excess_demand_by_area=hosp_peak_excess_demand_by_area,
        ICU_peak_excess_demand_by_area=ICU_peak_excess_demand_by_area,
        incidence_H_by_area_over_sims=incidence_H_by_area_over_sims,
        incidence_death_by_area_over_sims=incidence_death_by_area_over_sims,
        H_occup_by_area_over_sims=H_occup_by_area_over_sims,
        ICU_occup_by_area_over_sims=ICU_occup_by_area_over_sims,
        country_incidence_A_ts=country_incidence_A_ts,
        incidence_A_ts=incidence_A_ts,
        country_incidence_M_ts=country_incidence_M_ts,
        incidence_M_ts=incidence_M_ts,
        country_incidence_V_ts=country_incidence_V_ts,
        incidence_V_ts=incidence_V_ts,
        country_incidence_death_ts=country_incidence_death_ts,
        incidence_death_ts=incidence_death_ts,
        country_prevalence_H_ts=country_prevalence_H_ts,
        prevalence_H_ts=prevalence_H_ts,
        country_prevalence_ICU_ts=country_prevalence_ICU_ts,
        prevalence_ICU_ts=prevalence_ICU_ts)

# per scenario:
##

function school_effect(scenario_input_data)

    scenariodata = scenario_input_data

    # Number of symptomatic new/incident infections on day
    new_infections = [scenariodata.country_incidence_M_ts.med[t] for t in [30,45,60,90,180,360]] .+
                     [scenariodata.country_incidence_V_ts.med[t] for t in [30,45,60,90,180,360]]  # sum mild, severe and critical incidence
    new_infections_ts = scenariodata.country_incidence_M_ts.med .+ scenariodata.country_incidence_V_ts.med
    peak_value, peak_position = findmax(new_infections_ts)
    new_infections_peak_date = Date(2020,3,13) + Dates.Day(peak_position)
    new_infections = [new_infections;peak_value]  # add peak incidence to vector
    new_infections_peak_date = ["";"";"";"";"";"";new_infections_peak_date]  # add missing values in the first 6 positions
    new_symptomatic_infections_tuple = hcat(new_infections,new_infections_peak_date)

    # Number of incident/new mild cases on day
    new_infections = [scenariodata.country_incidence_M_ts.med[t] for t in [30,45,60,90,180,360]]
    new_infections_ts = scenariodata.country_incidence_M_ts.med
    peak_value, peak_position = findmax(new_infections_ts)
    new_infections_peak_date = Date(2020,3,13) + Dates.Day(peak_position)
    new_infections = [new_infections;peak_value]  # add peak incidence to vector
    new_infections_peak_date = ["";"";"";"";"";"";new_infections_peak_date]  # add missing values in the first 6 positions
    new_mild_infections_tuple = hcat(new_infections,new_infections_peak_date)

    # Number of incident/new severe cases on day
    new_infections = [scenariodata.country_incidence_V_ts.med[t] for t in [30,45,60,90,180,360]]
    new_infections_ts = scenariodata.country_incidence_V_ts.med
    peak_value, peak_position = findmax(new_infections_ts)
    new_infections_peak_date = Date(2020,3,13) + Dates.Day(peak_position)
    new_infections = [new_infections;peak_value]  # add peak incidence to vector
    new_infections_peak_date = ["";"";"";"";"";"";new_infections_peak_date]  # add missing values in the first 6 positions
    new_severe_infections_tuple = hcat(new_infections,new_infections_peak_date)

    # Number of incident/new critical cases on day
            # this wasnt saved! we leave the tuple blanc for now
    new_infections = ["";"";"";"";"";"";""]
    new_infections_peak_date = ["";"";"";"";"";"";""]  # add missing values in the first 6 positions
    new_critical_infections_tuple = hcat(new_infections,new_infections_peak_date)

    # Number of incident/new deaths on day
    new_deaths = [scenariodata.country_incidence_death_ts.med[t] for t in [30,45,60,90,180,360]]
    new_deaths_ts = scenariodata.country_incidence_death_ts.med
    peak_value, peak_position = findmax(new_deaths_ts)
    new_deaths_peak_date = Date(2020,3,13) + Dates.Day(peak_position)
    new_deaths = [new_deaths;peak_value]  # add peak incidence to vector
    new_deaths_peak_date = ["";"";"";"";"";"";new_deaths_peak_date]  # add missing values in the first 6 positions
    new_deaths_tuple = hcat(new_deaths,new_deaths_peak_date)


    # Prevalent severe cases on day (includes incident cases and severe cases not yet recovered or dead)
    prevalent_infections = [scenariodata.country_prevalence_H_ts.med[t] for t in [30,45,60,90,180,360]]
    prevalent_infections_ts = scenariodata.country_prevalence_H_ts.med
    peak_value, peak_position = findmax(prevalent_infections_ts)
    prevalent_infections_peak_date = Date(2020,3,13) + Dates.Day(peak_position)
    prevalent_infections = [prevalent_infections;peak_value]  # add peak prevalence to vector
    prevalent_infections_peak_date = ["";"";"";"";"";"";prevalent_infections_peak_date]  # add missing values in the first 6 positions
    prevalent_severe_infections_tuple = hcat(prevalent_infections,prevalent_infections_peak_date)

    # Prevalent critical cases on day (includes incident cases and those not yet recovered or dead)
    prevalent_infections = [scenariodata.country_prevalence_ICU_ts.med[t] for t in [30,45,60,90,180,360]]
    prevalent_infections_ts = scenariodata.country_prevalence_ICU_ts.med
    peak_value, peak_position = findmax(prevalent_infections_ts)
    prevalent_infections_peak_date = Date(2020,3,13) + Dates.Day(peak_position)
    prevalent_infections = [prevalent_infections;peak_value]  # add peak prevalence to vector
    prevalent_infections_peak_date = ["";"";"";"";"";"";prevalent_infections_peak_date]  # add missing values in the first 6 positions
    prevalent_critical_infections_tuple = hcat(prevalent_infections,prevalent_infections_peak_date)


    # 8. Cumulative infections by day
    new_infections = [cumsum(scenariodata.country_incidence_M_ts.med)[t] for t in [30,45,60,90,180,360]] .+
                     [cumsum(scenariodata.country_incidence_V_ts.med)[t] for t in [30,45,60,90,180,360]]  # sum mild, severe and critical incidence
    new_infections_ts = cumsum(scenariodata.country_incidence_M_ts.med .+ scenariodata.country_incidence_V_ts.med)
    peak_value, peak_position = findmax(new_infections_ts)
    new_infections_peak_date = Date(2020,3,13) + Dates.Day(peak_position)
    new_infections = [new_infections;peak_value]  # add peak incidence to vector
    new_infections_epidemic_end_ts = cumsum(scenariodata.country_incidence_M_ts.med .+ scenariodata.country_incidence_V_ts.med .+ scenariodata.country_incidence_A_ts.med) # consider even asymptomatics infections for end of epidemic judgement
    epidemic_end_time = findmax(new_infections_epidemic_end_ts)[2]
    new_infections = [new_infections;new_infections_ts[epidemic_end_time]]  # add cummulative M and V infections at end of epidemic
    epidemic_end_date = Date(2020,3,13) + Dates.Day(epidemic_end_time) #
    new_infections_peak_date = ["";"";"";"";"";"";new_infections_peak_date;epidemic_end_date]  # add missing values in the first 6 positions
    cum_symptomatic_infections_tuple = hcat(new_infections,new_infections_peak_date)

    # Cumulative mild cases by day
    new_infections = [cumsum(scenariodata.country_incidence_M_ts.med)[t] for t in [30,45,60,90,180,360]]
    new_infections_ts = cumsum(scenariodata.country_incidence_M_ts.med)
    peak_value, peak_position = findmax(new_infections_ts)
    new_infections_peak_date = Date(2020,3,13) + Dates.Day(peak_position)
    new_infections = [new_infections;peak_value]  # add peak incidence to vector
    new_infections_epidemic_end_ts = cumsum(scenariodata.country_incidence_M_ts.med .+ scenariodata.country_incidence_V_ts.med .+ scenariodata.country_incidence_A_ts.med) # consider even asymptomatics infections for end of epidemic judgement
    epidemic_end_time = findmax(new_infections_epidemic_end_ts)[2]
    new_infections = [new_infections;new_infections_ts[epidemic_end_time]]  # add cummulative M infections at end of epidemic
    epidemic_end_date = Date(2020,3,13) + Dates.Day(epidemic_end_time) #
    new_infections_peak_date = ["";"";"";"";"";"";new_infections_peak_date;epidemic_end_date]  # add missing values in the first 6 positions
    cum_mild_infections_tuple = hcat(new_infections,new_infections_peak_date)

    # Cumulative severe cases by day
    new_infections = [cumsum(scenariodata.country_incidence_V_ts.med)[t] for t in [30,45,60,90,180,360]]
    new_infections_ts = cumsum(scenariodata.country_incidence_V_ts.med)
    peak_value, peak_position = findmax(new_infections_ts)
    new_infections_peak_date = Date(2020,3,13) + Dates.Day(peak_position)
    new_infections = [new_infections;peak_value]  # add peak incidence to vector
    new_infections_epidemic_end_ts = cumsum(scenariodata.country_incidence_M_ts.med .+ scenariodata.country_incidence_V_ts.med .+ scenariodata.country_incidence_A_ts.med) # consider even asymptomatics infections for end of epidemic judgement
    epidemic_end_time = findmax(new_infections_epidemic_end_ts)[2]
    new_infections = [new_infections;new_infections_ts[epidemic_end_time]]  # add cummulative M infections at end of epidemic
    epidemic_end_date = Date(2020,3,13) + Dates.Day(epidemic_end_time) #
    new_infections_peak_date = ["";"";"";"";"";"";new_infections_peak_date;epidemic_end_date]  # add missing values in the first 6 positions
    cum_severe_infections_tuple = hcat(new_infections,new_infections_peak_date)

    #  Cumulative critical cases by day
      # this wasnt saved! we leave the tuple blanc for now
    new_infections = ["";"";"";"";"";"";"";""]
    new_infections_peak_date = ["";"";"";"";"";"";"";""]  # add missing values in the first 6 positions
    cum_critical_infections_tuple = hcat(new_infections,new_infections_peak_date)


    # Cumulative deaths by day
    new_deaths = [cumsum(scenariodata.country_incidence_death_ts.med)[t] for t in [30,45,60,90,180,360]]
    new_deaths_ts = cumsum(scenariodata.country_incidence_death_ts.med)
    peak_value, peak_position = findmax(new_deaths_ts)
    new_deaths_peak_date = Date(2020,3,13) + Dates.Day(peak_position)
    new_deaths = [new_deaths;peak_value]  # add peak incidence to vector
    new_deaths_epidemic_end_ts = cumsum(scenariodata.country_incidence_M_ts.med .+ scenariodata.country_incidence_V_ts.med .+ scenariodata.country_incidence_A_ts.med) # consider even asymptomatics infections for end of epidemic judgement
    epidemic_end_time = findmax(new_deaths_epidemic_end_ts)[2]
    new_deaths = [new_deaths;new_deaths_ts[epidemic_end_time]]  # add cummulative M infections at end of epidemic
    epidemic_end_date = Date(2020,3,13) + Dates.Day(epidemic_end_time) #
    new_deaths_peak_date = ["";"";"";"";"";"";new_deaths_peak_date;epidemic_end_date]  # add missing values in the first 6 positions
    cum_deaths_tuple = hcat(new_deaths,new_deaths_peak_date)

    out_tupple = vcat(new_symptomatic_infections_tuple,new_mild_infections_tuple,new_severe_infections_tuple,
                     new_critical_infections_tuple,new_deaths_tuple, prevalent_severe_infections_tuple,
                     prevalent_critical_infections_tuple,cum_symptomatic_infections_tuple,cum_mild_infections_tuple,
                     cum_severe_infections_tuple,cum_critical_infections_tuple,cum_deaths_tuple)

end

scenario_unmitigated =  school_effect(scenariodata_unmitigated)
scenario_current  = school_effect(scenariodata_full_intervention)
scenario_1a = school_effect(scenariodata_1a)
scenario_1b = school_effect(scenariodata_1b)
scenario_2a = school_effect(scenariodata_1c)
scenario_2b = school_effect(scenariodata_1d)
scenario_3a = fill("", 89, 2)  # these were redundant requests
scenario_3b = fill("", 89, 2)  # these were redubdant requests
scenario_4a = school_effect(scenariodata_candidates_june_90perc)
scenario_4b = school_effect(scenariodata_candidates_june_50perc)
scenario_5a = school_effect(scenariodata_candidates_august_90perc)
scenario_5b = school_effect(scenariodata_candidates_august_50perc)
scenario_6a = school_effect(scenariodata_primary_june_90perc)
scenario_6b = school_effect(scenariodata_primary_june_50perc)
scenario_7a = school_effect(scenariodata_secondary_june_90perc)
scenario_7b = school_effect(scenariodata_secondary_june_50perc)
scenario_8a = school_effect(scenariodata_tertiary_june_90perc)
scenario_8b = school_effect(scenariodata_tertiary_june_50perc)

all_school_effects = hcat(scenario_unmitigated,
                          scenario_current,
                          scenario_1a,
                          scenario_1b,
                          scenario_2a,
                          scenario_2b,
                          scenario_3a,
                          scenario_3b,
                          scenario_4a,
                          scenario_4b,
                          scenario_5a,
                          scenario_5b,
                          scenario_6a,
                          scenario_6b,
                          scenario_7a,
                          scenario_7b,
                          scenario_8a,
                          scenario_8b)

writedlm("reports/scenario_numbers.csv",all_school_effects,",")

##

function county_data(scenario_input_data)
    scenariodata = scenario_input_data

    # 1a. Number of new infections at day of peak in new/incident infections

    new_infections = scenariodata.incidence_M_ts.med .+ scenariodata.incidence_V_ts.med
    new_infections_by_county = zeros(Float64,47)
    for i=1:47
        new_infections_by_county[i] = findmax(new_infections[i,:])[1]
    end

    # 1b. Day of peak in new/incident infections
    day_of_peak_new_infections_by_county= zeros(Int64,47)
    for i=1:47
        day_of_peak_new_infections_by_county[i] = findmax(new_infections[i,:])[2]
    end

    # 2a. Number of new mild cases at peak in new/incident mild cases
    new_infections = scenariodata.incidence_M_ts.med
    new_mild_infections_by_county = zeros(Float64,47)
    for i=1:47
        new_mild_infections_by_county[i] = findmax(new_infections[i,:])[1]
    end

    # 2b. Day of peak in new/incident mild cases
    day_of_peak_new_mild_infections_by_county= zeros(Int64,47)
    for i=1:47
        day_of_peak_new_mild_infections_by_county[i] = findmax(new_infections[i,:])[2]
    end

    # 3a. Number of new severe cases at day of peak in new/incident severe cases
    new_infections = scenariodata.incidence_V_ts.med
    new_severe_infections_by_county = zeros(Float64,47)
    for i=1:47
        new_severe_infections_by_county[i] = findmax(new_infections[i,:])[1]
    end

    # 3b. Day of peak in new/incident severe cases
    day_of_peak_new_severe_infections_by_county= zeros(Int64,47)
    for i=1:47
        day_of_peak_new_severe_infections_by_county[i] = findmax(new_infections[i,:])[2]
    end


    # 4a. Number of new critical cases at day of peak in  new/incident critical cases
    new_critical_infections_by_county = fill("",47)  # was not extracted from model

    # 4b. Day of peak in new/incident critical cases
    day_of_peak_new_critical_infections_by_county = fill("",47)


    # 5a. Number of deaths at day of peak in new/incident deaths
    new_infections = scenariodata.incidence_death_ts.med
    new_deaths_by_county = zeros(Float64,47)
    for i=1:47
        new_deaths_by_county[i] = findmax(new_infections[i,:])[1]
    end

    # 5b. Day of peak in new/incident deaths
    day_of_peak_new_deaths_by_county= zeros(Int64,47)
    for i=1:47
        day_of_peak_new_deaths_by_county[i] = findmax(new_infections[i,:])[2]
    end


    # 6a. Number of severe cases at peak of prevalent severe cases  (includes incident cases and severe cases not yet recovered or dead)
    new_infections = scenariodata.prevalence_H_ts.med
    prevalent_severe_cases_by_county = zeros(Float64,47)
    for i=1:47
        prevalent_severe_cases_by_county[i] = findmax(new_infections[i,:])[1]
    end

    # 6b. Day of peak in prevalent severe cases
    day_of_peak_prevalent_severe_cases_by_county= zeros(Int64,47)
    for i=1:47
        day_of_peak_prevalent_severe_cases_by_county[i] = findmax(new_infections[i,:])[2]
    end

    #7a. Number of critical cases at peak of prevalent critical cases  (includes incident cases and those not yet recovered or dead)
    new_infections = scenariodata.prevalence_ICU_ts.med
    prevalent_critical_cases_by_county = zeros(Float64,47)
    for i=1:47
        prevalent_critical_cases_by_county[i] = findmax(new_infections[i,:])[1]
    end

    # 7b. Day of peak in prevalent critical cases
    day_of_peak_prevalent_critical_cases_by_county= zeros(Int64,47)
    for i=1:47
        day_of_peak_prevalent_critical_cases_by_county[i] = findmax(new_infections[i,:])[2]
    end

    # 8. Cumulative infections at 360 days from first reported case
    new_infections = scenariodata.incidence_M_ts.med .+ scenariodata.incidence_V_ts.med
    cum_infections_by_county = zeros(Float64,47)
    for i=1:47
        cum_infections_by_county[i] = cumsum(new_infections[i,:])[360]
    end

    # 9. Cumulative mild cases at 360 days from first reported case
    new_infections = scenariodata.incidence_M_ts.med
    cum_mild_infections_by_county = zeros(Float64,47)
    for i=1:47
        cum_mild_infections_by_county[i] = cumsum(new_infections[i,:])[360]
    end

    # 10. Cumulative severe cases at 360 days from first reported case
    new_infections = scenariodata.incidence_V_ts.med
    cum_severe_infections_by_county = zeros(Float64,47)
    for i=1:47
        cum_severe_infections_by_county[i] = cumsum(new_infections[i,:])[360]
    end

    # 11. Cumulative critical cases at 360 days from first reported case
    cum_critical_infections_by_county = fill("",47) # was not extraced during model return

    # 12. Cumulative deaths at 360 days from first reported case
    new_infections = scenariodata.incidence_death_ts.med
    cum_deaths_by_county = zeros(Float64,47)
    for i=1:47
        cum_deaths_by_county[i] = cumsum(new_infections[i,:])[360]
    end

    # 13b. Day at which  hospital bed capacity for severe cases is first exceeded
    day_hospital_capacity_exceeded  = fill("",47) # not extracted during model run

    # 14b. Day at which ventilator capacity  for critcial cases is first exceeded
    day_ICU_capacity_exceeded  = fill("",47) # not extracted during model run

    #

    all_tuples = fill("",23,47)
    all_tuples[1,:] = string.(new_infections_by_county)
    all_tuples[2,:] = string.(day_of_peak_new_infections_by_county)
    all_tuples[3,:] = string.(new_mild_infections_by_county)
    all_tuples[4,:] = string.(day_of_peak_new_mild_infections_by_county)
    all_tuples[5,:] = string.(new_severe_infections_by_county)
    all_tuples[6,:] = string.(day_of_peak_new_severe_infections_by_county)
    all_tuples[7,:] = new_critical_infections_by_county
    all_tuples[8,:] = day_of_peak_new_critical_infections_by_county
    all_tuples[9,:] = string.(new_deaths_by_county)
    all_tuples[10,:] = string.(day_of_peak_new_deaths_by_county)
    all_tuples[11,:] = string.(prevalent_severe_cases_by_county)
    all_tuples[12,:] = string.(day_of_peak_prevalent_severe_cases_by_county)
    all_tuples[13,:] = string.(prevalent_critical_cases_by_county)
    all_tuples[14,:] = string.(day_of_peak_prevalent_critical_cases_by_county)
    all_tuples[15,:] = string.(cum_infections_by_county)
    all_tuples[16,:] = string.(cum_mild_infections_by_county)
    all_tuples[17,:] = string.(cum_severe_infections_by_county)
    all_tuples[18,:] = cum_critical_infections_by_county
    all_tuples[19,:] = string.(cum_deaths_by_county)
    all_tuples[20,:] = fill("",47) # available hospital capacity: already prefilled in excel
    all_tuples[21,:] = day_hospital_capacity_exceeded
    all_tuples[22,:] = fill("",47) # available ventilator capacity: already prefilled in excel
    all_tuples[23,:] = day_ICU_capacity_exceeded

    return all_tuples
end

scenario_unmitigated =  county_data(scenariodata_unmitigated)
scenario_current  = county_data(scenariodata_full_intervention)
scenario_1a = county_data(scenariodata_1a)
scenario_1b = county_data(scenariodata_1b)
scenario_2a = county_data(scenariodata_1c)
scenario_2b = county_data(scenariodata_1d)
scenario_3a = fill("", 89, 2)  # these were redundant requests
scenario_3b = fill("", 89, 2)  # these were redubdant requests
scenario_4a = county_data(scenariodata_candidates_june_90perc)
scenario_4b = county_data(scenariodata_candidates_june_50perc)
scenario_5a = county_data(scenariodata_candidates_august_90perc)
scenario_5b = county_data(scenariodata_candidates_august_50perc)
scenario_6a = county_data(scenariodata_primary_june_90perc)
scenario_6b = county_data(scenariodata_primary_june_50perc)
scenario_7a = county_data(scenariodata_secondary_june_90perc)
scenario_7b = county_data(scenariodata_secondary_june_50perc)
scenario_8a = county_data(scenariodata_tertiary_june_90perc)
scenario_8b = county_data(scenariodata_tertiary_june_50perc)

all_scenarios = ["scenario_unmitigated",
                          "scenario_current",
                          "scenario_1a",
                          "scenario_1b",
                          "scenario_2a",
                          "scenario_2b",
                          "scenario_3a",
                          "scenario_3b",
                          "scenario_4a",
                          "scenario_4b",
                          "scenario_5a",
                          "scenario_5b",
                          "scenario_6a",
                          "scenario_6b",
                          "scenario_7a",
                          "scenario_7b",
                          "scenario_8a",
                          "scenario_8b"]


for i = 1:size(all_scenarios)[1]
    writedlm("reports/county_numbers_"*all_scenarios[i]*".csv",Meta.parse(all_scenarios[i)],",")
end
