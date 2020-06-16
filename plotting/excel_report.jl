using JLD2,Plots,StatsPlots,CSV,DelimitedFiles, Dates
include("get_scenario_data.jl");

hosp_capacity = CSV.read(joinpath(homedir(),"Documents/Covid-19/jl_models/KenyaCoVOutputs/Health_system_capacity_data_Kenya.csv"))
spare_capacity_H_by_county = (hosp_capacity[:,2].*hosp_capacity[:,5])[2:end]
spare_capacity_ICU_by_county = (hosp_capacity[:,3].*hosp_capacity[:,6])[2:end]
govt_hosp_capacity = CSV.read(joinpath(homedir(),"Documents/Covid-19/jl_models/data/Beds in Public Hospitals-04062020.csv"))
govt_spare_capacity_H_by_county =  govt_hosp_capacity[:,2].*(hosp_capacity[:,5][2:end])


# per scenario:
##

function school_effect(scenario_input_data)

    scenariodata = scenario_input_data

    # Number of ALL new/incident infections on day
    new_infections = [scenariodata.country_incidence_M_ts.med[t] for t in [30,45,60,90,180,360]] .+
                     [scenariodata.country_incidence_V_ts.med[t] for t in [30,45,60,90,180,360]] .+ # sum ALL incidence
                     [scenariodata.country_incidence_A_ts.med[t] for t in [30,45,60,90,180,360]]
    new_infections_ts = scenariodata.country_incidence_M_ts.med .+ scenariodata.country_incidence_V_ts.med .+ scenariodata.country_incidence_A_ts.med
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
                     [cumsum(scenariodata.country_incidence_V_ts.med)[t] for t in [30,45,60,90,180,360]] .+
                     [cumsum(scenariodata.country_incidence_A_ts.med)[t] for t in [30,45,60,90,180,360]]  # ALLincidence
    new_infections_ts = cumsum(scenariodata.country_incidence_M_ts.med .+ scenariodata.country_incidence_V_ts.med .+ scenariodata.country_incidence_A_ts.med)
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

    new_infections = scenariodata.incidence_M_ts.med .+ scenariodata.incidence_V_ts.med .+ scenariodata.incidence_A_ts.med
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
    new_infections = scenariodata.incidence_M_ts.med .+ scenariodata.incidence_V_ts.med .+ scenariodata.incidence_A_ts.med
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
    writedlm("reports/county_numbers_"*all_scenarios[i]*".csv",eval(Meta.parse(all_scenarios[i])),",")
end

# Correct plotting x-axis

function give_plots_for_county_temp(output,areas_to_plot,scenario_tag,areanames)

    monthdates = [Date(2020,3,1) + Month(i) for i = 1:21 ]
    monthnames = [monthname(d)[1:3]*"-$(year(d)-2000)" for d in monthdates]
    tick_times = [(d - Date(2020,3,13)).value for d in monthdates]

    median_H = output.incidence_V_ts.med[areas_to_plot,:]
    lb_H = output.incidence_V_ts.lpred[areas_to_plot,:]
    ub_H = output.incidence_V_ts.upred[areas_to_plot,:]
    median_D = output.incidence_death_ts.med[areas_to_plot,:]
    lb_D = output.incidence_death_ts.lpred[areas_to_plot,:]
    ub_D = output.incidence_death_ts.upred[areas_to_plot,:]
    median_H_occup = output.prevalence_H_ts.med[areas_to_plot,:]
    lb_H_occup = output.prevalence_H_ts.lpred[areas_to_plot,:]
    ub_H_occup = output.prevalence_H_ts.upred[areas_to_plot,:]
    median_ICU_occup = output.prevalence_ICU_ts.med[areas_to_plot,:]
    lb_ICU_occup = output.prevalence_ICU_ts.lpred[areas_to_plot,:]
    ub_ICU_occup = output.prevalence_ICU_ts.upred[areas_to_plot,:]

    area_name = ["Rest of Kenya"]
    if length(areas_to_plot) == 1
        area_name = areanames[areas_to_plot]
    end

    T= size(median_H)[1]

    plt_incidence = plot(0:(T-1),median_H,
                        lab = "Hospitalisations",
                        lw = 3,color = :red,
                        ribbon = (median_H.-lb_H,ub_H .- median_H),
                        fillalpha = 0.25,
                        xticks = (tick_times,monthnames),
                        title = "Incidence of disease - "*area_name*scenario_tag,
                        size = (1000,600))

    plot!(plt_incidence,0:(T-1),median_D,
                        lab = "Deaths",
                        lw = 3,color = :black,
                        ribbon = (median_D.-lb_D,ub_D .- median_D),
                        fillalpha = 0.5 )

    plt_health_usage = plot(0:(T-1),median_H_occup,
                            lab = "hospital beds",
                            lw = 3,color = :blue,
                            ribbon = (median_H_occup.-lb_H_occup,ub_H_occup .- median_H_occup),
                            fillalpha = 0.25,
                            xticks = (tick_times,monthnames),
                            ylabel = "Daily Occupancy",
                            title = "Health system usage - "*area_name*scenario_tag,
                            size = (1000,600))


    plot!(plt_health_usage,0:(T-1),median_ICU_occup,
                            lab = "ICU beds",
                            lw = 3,color = :green,
                            ribbon = (median_ICU_occup.-lb_ICU_occup,ub_ICU_occup .- median_ICU_occup),
                            fillalpha = 0.25)
    if length(areas_to_plot) == 1
        plot!(plt_health_usage,[0,T-1],[spare_capacity_H_by_county[areas_to_plot][1],spare_capacity_H_by_county[areas_to_plot][1]],
                lw = 2,ls = :dash,lab = "spare hosp. capacity",color = :blue)
        plot!(plt_health_usage,[0,T-1],[spare_capacity_H_by_county[areas_to_plot][1],spare_capacity_H_by_county[areas_to_plot][1]].*1.2,
                lw = 2,ls = :dash,lab = "20% surge capacity",color = :red)
        plot!(plt_health_usage,[0,T-1],[govt_spare_capacity_H_by_county[areas_to_plot][1],govt_spare_capacity_H_by_county[areas_to_plot][1]],
                lw = 2,ls = :dash,lab = "spare govt. hosp. capacity",color = :black)
        plot!(plt_health_usage,[0,T-1],[spare_capacity_ICU_by_county[areas_to_plot][1],spare_capacity_ICU_by_county[areas_to_plot][1]],
                lw = 2,ls = :dash,lab = "spare ICU capacity",color = :green)
    end

    return plt_incidence,plt_health_usage
end


function correct_county_plots(output,scenario_tag,areanames,simulation_tag)
    plt_incidence_mombasa,plt_health_usage_mombasa=give_plots_for_county_temp(output,28,scenario_tag,areanames)
    plt_incidence_nairobi,plt_health_usage_nairobi=give_plots_for_county_temp(output,30,scenario_tag,areanames)
    savefig(plt_incidence_mombasa,"reports/report"*simulation_tag*"/incidence_mombasa"*simulation_tag*".png")
    savefig(plt_health_usage_mombasa,"reports/report"*simulation_tag*"/healthsystem_mombasa"*simulation_tag*".png")
    savefig(plt_incidence_nairobi,"reports/report"*simulation_tag*"/incidence_nairobi"*simulation_tag*".png")
    savefig(plt_health_usage_nairobi,"reports/report"*simulation_tag*"/healthsystem_nairobi"*simulation_tag*".png")

    return nothing
end


correct_county_plots(scenariodata_unmitigated," (Unmitigated)",counties.county,"_unmitigated")
correct_county_plots(scenariodata_tertiary_june_90perc," (Tertiary only 90%)",counties.county,"_tertiary_june_90perc")
correct_county_plots(scenariodata_secondary_june_90perc," (Secondary only 90%)",counties.county,"_secondary_june_90perc")
correct_county_plots(scenariodata_primary_june_90perc," (Primary only 90%)",counties.county,"_primary_june_90perc")
correct_county_plots(scenariodata_tertiary_june_50perc," (Tertiary only 50%)",counties.county,"_tertiary_june_50perc")
correct_county_plots(scenariodata_secondary_june_50perc," (Secondary only 50%)",counties.county,"_secondary_june_50perc")
correct_county_plots(scenariodata_primary_june_50perc," (Primary only 50%)",counties.county,"_primary_june_50perc")
correct_county_plots(scenariodata_full_intervention," (Full intervention)",counties.county,"_full_intervention")
correct_county_plots(scenariodata_1a," (June opening, contacts at 50%)",counties.county,"_scenario_1a")
correct_county_plots(scenariodata_1b," (June opening, contacts at 90%)",counties.county,"_scenario_1b")
correct_county_plots(scenariodata_1c," (August opening, contacts at 50%)",counties.county,"_scenario_1c")
correct_county_plots(scenariodata_1d," (August opening, contacts at 90%)",counties.county,"_scenario_1d")
correct_county_plots(scenariodata_end_regional_lockdown," (End movement restrictions 15th June)",counties.county,"_end_regional_lockdown")


## Plot vulnarability index


function plot_ranked_bars_health_usage_and_vulnarability(output,scenario_tag,areanames)
    vulnarability_index = CSV.read(joinpath(homedir(),"Documents/Covid-19/jl_models/data/Ranked vulnerability indices_county (0406120).csv"))[:,7]
    median_hospital_exceed = output.hosp_peak_excess_demand_by_area.med
    median_ICU_exceed = output.ICU_peak_excess_demand_by_area.med
    nc = length(areanames)

    I = sortperm(median_hospital_exceed)
    plt_HU = bar(median_hospital_exceed[I]*100,orientations = :horizonal,
                yticks = (1:47,areanames[I]),size = (700,550),lab ="",
                xlabel = "% Hospital available bed demand",title = "Hospital demand at peak"*scenario_tag,
                fill_z=vulnarability_index[I],c=ColorGradient(:heat))
    scatter!(plt_HU,median_hospital_exceed[I]*100,1:nc,
                xerror = ((median_hospital_exceed[I] .- output.hosp_peak_excess_demand_by_area.lpred[I])*100
                            ,(output.hosp_peak_excess_demand_by_area.upred[I] .- median_hospital_exceed[I])*100 ),
                ms = 0.,color = :black,lab ="")

    I = sortperm(median_ICU_exceed)
    plt_ICU = bar(median_ICU_exceed[I,1],orientations = :horizonal,
                yticks = (1:nc,areanames[I]),size = (700,550),lab ="",
                xlabel = "ICU bed excess demand (numbers)",title = "ICU demand at peak"*scenario_tag,
                fill_z=vulnarability_index[I],c=ColorGradient(:heat) )
    scatter!(plt_ICU,median_ICU_exceed[I],1:nc,
                xerror = ((median_ICU_exceed[I] .- output.ICU_peak_excess_demand_by_area.lpred[I]),
                            (output.ICU_peak_excess_demand_by_area.upred[I] .- median_ICU_exceed[I]) ),
                ms = 0.,color = :black,lab ="")



        return plt_HU,plt_ICU
end

function peak_H_ICU_vulnarability_plots(output,scenario_tag,areanames,simulation_tag)
    plt_ranked_HU,plt_ranked_ICU = plot_ranked_bars_health_usage_and_vulnarability(output,scenario_tag,areanames)
    savefig(plt_ranked_HU,"reports/report"*simulation_tag*"/peak_hospital_usage_and_vulnarability_by_county"*simulation_tag*".png")
    savefig(plt_ranked_ICU,"reports/report"*simulation_tag*"/peak_ICU_usage_and_vulnarability_by_county"*simulation_tag*".png")
end

peak_H_ICU_vulnarability_plots(scenariodata_unmitigated," (Unmitigated)",counties.county,"_unmitigated")
peak_H_ICU_vulnarability_plots(scenariodata_tertiary_june_90perc," (Tertiary only 90%)",counties.county,"_tertiary_june_90perc")
peak_H_ICU_vulnarability_plots(scenariodata_secondary_june_90perc," (Secondary only 90%)",counties.county,"_secondary_june_90perc")
peak_H_ICU_vulnarability_plots(scenariodata_primary_june_90perc," (Primary only 90%)",counties.county,"_primary_june_90perc")
peak_H_ICU_vulnarability_plots(scenariodata_tertiary_june_50perc," (Tertiary only 50%)",counties.county,"_tertiary_june_50perc")
peak_H_ICU_vulnarability_plots(scenariodata_secondary_june_50perc," (Secondary only 50%)",counties.county,"_secondary_june_50perc")
peak_H_ICU_vulnarability_plots(scenariodata_primary_june_50perc," (Primary only 50%)",counties.county,"_primary_june_50perc")
peak_H_ICU_vulnarability_plots(scenariodata_full_intervention," (Full intervention)",counties.county,"_full_intervention")
peak_H_ICU_vulnarability_plots(scenariodata_1a," (June opening, contacts at 50%)",counties.county,"_scenario_1a")
peak_H_ICU_vulnarability_plots(scenariodata_1b," (June opening, contacts at 90%)",counties.county,"_scenario_1b")
peak_H_ICU_vulnarability_plots(scenariodata_1c," (August opening, contacts at 50%)",counties.county,"_scenario_1c")
peak_H_ICU_vulnarability_plots(scenariodata_1d," (August opening, contacts at 90%)",counties.county,"_scenario_1d")
peak_H_ICU_vulnarability_plots(scenariodata_end_regional_lockdown," (End movement restrictions 15th June)",counties.county,"_end_regional_lockdown")
