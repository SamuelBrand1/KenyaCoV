
function extract_information_from_simulations(sims)
    n = length(sims)
    nc,T = size(sims[1].hosp_occup_by_area_ts)
    nc,na = size(sims[1].total_ICU_by_area_and_age)

    hosp_by_area_over_sims = VectorOfArray([sum(sims[k].total_hosp_by_area_and_age,dims=2) for k = 1:n])[:,1,:]
    deaths_by_area_over_sims = VectorOfArray([sum(sims[k].total_deaths_by_area_and_age,dims=2) for k = 1:n])[:,1,:]
    hosp_by_age_over_sims = VectorOfArray([sum(sims[k].total_hosp_by_area_and_age,dims=1) for k = 1:n])[1,:,:]
    deaths_by_age_over_sims = VectorOfArray([sum(sims[k].total_deaths_by_area_and_age,dims=1) for k = 1:n])[1,:,:]
    incidence_H_by_area_over_sims = VectorOfArray([sims[k].incidence_H for k = 1:n])[:,:,:]
    incidence_death_by_area_over_sims = VectorOfArray([sims[k].death_incidence_by_area_ts for k = 1:n])[:,:,:]
    H_occup_by_area_over_sims = VectorOfArray([sims[k].hosp_occup_by_area_ts for k = 1:n])[:,:,:]
    ICU_occup_by_area_over_sims = VectorOfArray([sims[k].ICU_occup_by_area_ts for k = 1:n])[:,:,:]

    country_incidence_A = VectorOfArray([sum(sims[k].incidence_A,dims=1) for k = 1:n])[1,:,:]
    incidence_A_by_area_and_sims = VectorOfArray([sims[k].incidence_A for k = 1:n])[:,:,:]

    incidence_A_by_area_med = zeros(nc,T)
    incidence_A_by_area_lpred = zeros(nc,T)
    incidence_A_by_area_upred = zeros(nc,T)
    for cn = 1:nc,t = 1:T
        incidence_A_by_area_med[cn,t] = median(incidence_A_by_area_and_sims[cn,t,:])
        incidence_A_by_area_lpred[cn,t] = quantile(incidence_A_by_area_and_sims[cn,t,:],0.025)
        incidence_A_by_area_upred[cn,t] = quantile(incidence_A_by_area_and_sims[cn,t,:],0.975)
    end

    country_incidence_M = VectorOfArray([sum(sims[k].incidence_M,dims=1) for k = 1:n])[1,:,:]
    incidence_M_by_area_and_sims = VectorOfArray([sims[k].incidence_M for k = 1:n])[:,:,:]

    incidence_M_by_area_med = zeros(nc,T)
    incidence_M_by_area_lpred = zeros(nc,T)
    incidence_M_by_area_upred = zeros(nc,T)
    for cn = 1:nc,t = 1:T
        incidence_M_by_area_med[cn,t] = median(incidence_M_by_area_and_sims[cn,t,:])
        incidence_M_by_area_lpred[cn,t] = quantile(incidence_M_by_area_and_sims[cn,t,:],0.025)
        incidence_M_by_area_upred[cn,t] = quantile(incidence_M_by_area_and_sims[cn,t,:],0.975)
    end


    country_incidence_V = VectorOfArray([sum(sims[k].incidence_V,dims=1) for k = 1:n])[1,:,:]
    incidence_V_by_area_and_sims = VectorOfArray([sims[k].incidence_V for k = 1:n])[:,:,:]

    incidence_V_by_area_med = zeros(nc,T)
    incidence_V_by_area_lpred = zeros(nc,T)
    incidence_V_by_area_upred = zeros(nc,T)
    for cn = 1:nc,t = 1:T
        incidence_V_by_area_med[cn,t] = median(incidence_V_by_area_and_sims[cn,t,:])
        incidence_V_by_area_lpred[cn,t] = quantile(incidence_V_by_area_and_sims[cn,t,:],0.025)
        incidence_V_by_area_upred[cn,t] = quantile(incidence_V_by_area_and_sims[cn,t,:],0.975)
    end

    country_incidence_death = VectorOfArray([sum(sims[k].death_incidence_by_area_ts,dims=1) for k = 1:n])[1,:,:]
    incidence_death_by_area_and_sims = VectorOfArray([sims[k].death_incidence_by_area_ts for k = 1:n])[:,:,:]

    incidence_death_by_area_med = zeros(nc,T)
    incidence_death_by_area_lpred = zeros(nc,T)
    incidence_death_by_area_upred = zeros(nc,T)
    for cn = 1:nc,t = 1:T
        incidence_death_by_area_med[cn,t] = median(incidence_death_by_area_and_sims[cn,t,:])
        incidence_death_by_area_lpred[cn,t] = quantile(incidence_death_by_area_and_sims[cn,t,:],0.025)
        incidence_death_by_area_upred[cn,t] = quantile(incidence_death_by_area_and_sims[cn,t,:],0.975)
    end

    country_prevalence_H = VectorOfArray([sum(sims[k].hosp_occup_by_area_ts,dims=1) for k = 1:n])[1,:,:]
    prevalence_H_by_area_and_sims = VectorOfArray([sims[k].hosp_occup_by_area_ts for k = 1:n])[:,:,:]

    prevalence_H_by_area_med = zeros(nc,T)
    prevalence_H_by_area_lpred = zeros(nc,T)
    prevalence_H_by_area_upred = zeros(nc,T)
    for cn = 1:nc,t = 1:T
        prevalence_H_by_area_med[cn,t] = median(prevalence_H_by_area_and_sims[cn,t,:])
        prevalence_H_by_area_lpred[cn,t] = quantile(prevalence_H_by_area_and_sims[cn,t,:],0.025)
        prevalence_H_by_area_upred[cn,t] = quantile(prevalence_H_by_area_and_sims[cn,t,:],0.975)
    end

    country_prevalence_ICU = VectorOfArray([sum(sims[k].ICU_occup_by_area_ts,dims=1) for k = 1:n])[1,:,:]
    prevalence_ICU_by_area_and_sims = VectorOfArray([sims[k].ICU_occup_by_area_ts for k = 1:n])[:,:,:]

    prevalence_ICU_by_area_med = zeros(nc,T)
    prevalence_ICU_by_area_lpred = zeros(nc,T)
    prevalence_ICU_by_area_upred = zeros(nc,T)
    for cn = 1:nc,t = 1:T
        prevalence_ICU_by_area_med[cn,t] = median(prevalence_ICU_by_area_and_sims[cn,t,:])
        prevalence_ICU_by_area_lpred[cn,t] = quantile(prevalence_ICU_by_area_and_sims[cn,t,:],0.025)
        prevalence_ICU_by_area_upred[cn,t] = quantile(prevalence_ICU_by_area_and_sims[cn,t,:],0.975)
    end

    peak_hosp_by_area_by_sim = zeros(nc,n)
    peak_ICU_by_area_by_sim = zeros(nc,n)

    for cn = 1:nc,k = 1:n
        peak_hosp_by_area_by_sim[cn,k] = maximum(sims[k].hosp_occup_by_area_ts[cn,:])
        peak_ICU_by_area_by_sim[cn,k] = maximum(sims[k].ICU_occup_by_area_ts[cn,:])
    end

    country_incidence_A_ts = (med = [median(country_incidence_A[t,:]) for t = 1:T],
                            lpred = [quantile(country_incidence_A[t,:],0.025) for t = 1:T],
                            upred = [quantile(country_incidence_A[t,:],0.025) for t = 1:T])
    incidence_A_ts = (med = incidence_A_by_area_med,
                    lpred = incidence_A_by_area_lpred,
                    upred = incidence_A_by_area_upred)

    country_incidence_M_ts = (med = [median(country_incidence_M[t,:]) for t = 1:T],
                        lpred = [quantile(country_incidence_M[t,:],0.025) for t = 1:T],
                        upred = [quantile(country_incidence_M[t,:],0.025) for t = 1:T])
    incidence_M_ts = (med = incidence_M_by_area_med,
                lpred = incidence_M_by_area_lpred,
                upred = incidence_M_by_area_upred)

    country_incidence_V_ts = (med = [median(country_incidence_V[t,:]) for t = 1:T],
                        lpred = [quantile(country_incidence_V[t,:],0.025) for t = 1:T],
                        upred = [quantile(country_incidence_V[t,:],0.025) for t = 1:T])
    incidence_V_ts = (med = incidence_V_by_area_med,
                lpred = incidence_V_by_area_lpred,
                upred = incidence_V_by_area_upred)

    country_incidence_death_ts = (med = [median(country_incidence_death[t,:]) for t = 1:T],
                        lpred = [quantile(country_incidence_death[t,:],0.025) for t = 1:T],
                        upred = [quantile(country_incidence_death[t,:],0.025) for t = 1:T])

    incidence_death_ts = (med = incidence_death_by_area_med,
                lpred = incidence_death_by_area_lpred,
                upred = incidence_death_by_area_upred)

    country_prevalence_H_ts = (med = [median(country_prevalence_H[t,:]) for t = 1:T],
                        lpred = [quantile(country_prevalence_H[t,:],0.025) for t = 1:T],
                        upred = [quantile(country_prevalence_H[t,:],0.025) for t = 1:T])

    prevalence_H_ts = (med = prevalence_H_by_area_med,
                lpred = prevalence_H_by_area_lpred,
                upred = prevalence_H_by_area_upred)


    country_prevalence_ICU_ts = (med = [median(country_prevalence_ICU[t,:]) for t = 1:T],
                        lpred = [quantile(country_prevalence_ICU[t,:],0.025) for t = 1:T],
                        upred = [quantile(country_prevalence_ICU[t,:],0.025) for t = 1:T])

    prevalence_ICU_ts = (med = prevalence_ICU_by_area_med,
                lpred = prevalence_ICU_by_area_lpred,
                upred = prevalence_ICU_by_area_upred)

    total_severe_cases = (med = median(sum(hosp_by_area_over_sims,dims=1)),
                            lpred = quantile(sum(hosp_by_area_over_sims,dims=1)[:],0.025),
                            upred = quantile(sum(hosp_by_area_over_sims[:],dims=1),0.975))

    severe_cases_by_area = (med = [median(hosp_by_area_over_sims[cn,:]) for cn = 1:nc],
                            lpred = [quantile(hosp_by_area_over_sims[cn,:],0.025) for cn = 1:nc],
                            upred = [quantile(hosp_by_area_over_sims[cn,:],0.975) for cn = 1:nc])
    severe_cases_by_age = (med = [median(hosp_by_age_over_sims[a,:]) for a = 1:na],
                            lpred = [quantile(hosp_by_age_over_sims[a,:],0.025) for a = 1:na],
                            upred = [quantile(hosp_by_age_over_sims[a,:],0.975) for a = 1:na])
    total_deaths = (med = median(sum(deaths_by_area_over_sims,dims=1)),
                            lpred = quantile(sum(deaths_by_area_over_sims,dims=1)[:],0.025),
                            upred = quantile(sum(deaths_by_area_over_sims[:],dims=1),0.975))
    deaths_by_area = (med = [median(deaths_by_area_over_sims[cn,:]) for cn = 1:nc],
                            lpred = [quantile(deaths_by_area_over_sims[cn,:],0.025) for cn = 1:nc],
                            upred = [quantile(deaths_by_area_over_sims[cn,:],0.975) for cn = 1:nc])
    deaths_by_age = (med = [median(deaths_by_age_over_sims[a,:]) for a = 1:na],
                            lpred = [quantile(deaths_by_age_over_sims[a,:],0.025) for a = 1:na],
                            upred = [quantile(deaths_by_age_over_sims[a,:],0.975) for a = 1:na])
    hosp_peak_excess_demand_by_area = (med = [median(peak_hosp_by_area_by_sim[cn,:]/KenyaCoV.spare_capacity_H_by_county[cn]) for cn = 1:nc],
                            lpred = [quantile(peak_hosp_by_area_by_sim[cn,:]/KenyaCoV.spare_capacity_H_by_county[cn],0.025) for cn = 1:nc],
                            upred = [quantile(peak_hosp_by_area_by_sim[cn,:]/KenyaCoV.spare_capacity_H_by_county[cn],0.975) for cn = 1:nc])

    ICU_peak_excess_demand_by_area = (med = [median(peak_ICU_by_area_by_sim[cn,:] .- KenyaCoV.spare_capacity_ICU_by_county[cn]) for cn = 1:nc],
                            lpred = [quantile(peak_ICU_by_area_by_sim[cn,:] .- KenyaCoV.spare_capacity_ICU_by_county[cn],0.025) for cn = 1:nc],
                            upred = [quantile(peak_ICU_by_area_by_sim[cn,:] .- KenyaCoV.spare_capacity_ICU_by_county[cn],0.975) for cn = 1:nc])

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
end


function generate_report(output,model_str,simulation_tag,scenario_tag,areanames;make_new_directory::Bool = false)
    if make_new_directory
        mkdir("reports/report"*simulation_tag)
    end
    scenariodata = (total_severe_cases=output.total_severe_cases,
            severe_cases_by_area=output.severe_cases_by_area,
            severe_cases_by_age=output.severe_cases_by_age,
            total_deaths=output.total_deaths,
            deaths_by_area=output.deaths_by_area,
            deaths_by_age=output.deaths_by_age,
            hosp_peak_excess_demand_by_area=output.hosp_peak_excess_demand_by_area,
            ICU_peak_excess_demand_by_area=output.ICU_peak_excess_demand_by_area,
            country_incidence_A_ts=output.country_incidence_A_ts,
            incidence_A_ts=output.incidence_A_ts,
            country_incidence_M_ts=output.country_incidence_M_ts,
            incidence_M_ts=output.incidence_M_ts,
            country_incidence_V_ts=output.country_incidence_V_ts,
            incidence_V_ts=output.incidence_V_ts,
            country_incidence_death_ts=output.country_incidence_death_ts,
            incidence_death_ts=output.incidence_death_ts,
            country_prevalence_H_ts=output.country_prevalence_H_ts,
            prevalence_H_ts=output.prevalence_H_ts,
            country_prevalence_ICU_ts=output.country_prevalence_ICU_ts,
            prevalence_ICU_ts=output.prevalence_ICU_ts,
            model_str=model_str)
    @save("reports/report"*simulation_tag*"/scenario_data"*simulation_tag*".jld2",scenariodata)

    plt_incidence_mombasa,plt_health_usage_mombasa = KenyaCoV.give_plots_for_county(output,[28],scenario_tag,areanames)
    savefig(plt_incidence_mombasa,"reports/report"*simulation_tag*"/incidence_mombasa"*simulation_tag*".png")
    savefig(plt_health_usage_mombasa,"reports/report"*simulation_tag*"/healthsystem_mombasa"*simulation_tag*".png")
    plt_incidence_nairobi,plt_health_usage_nairobi = KenyaCoV.give_plots_for_county(output,[30],scenario_tag,areanames)
    savefig(plt_incidence_nairobi,"reports/report"*simulation_tag*"/incidence_nairobi"*simulation_tag*".png")
    savefig(plt_health_usage_nairobi,"reports/report"*simulation_tag*"/healthsystem_nairobi"*simulation_tag*".png")
    plt_incidence_restofcountry,plt_health_usage_restofcountry = KenyaCoV.give_plots_for_county(output,setdiff(1:47,[28,30]),scenario_tag,areanames)
    savefig(plt_incidence_restofcountry,"reports/report"*simulation_tag*"/plt_incidence_restofcountry"*simulation_tag*".png")
    savefig(plt_health_usage_restofcountry,"reports/report"*simulation_tag*"/healthsystem_restofcountry"*simulation_tag*".png")

    plt_ranked_HU,plt_ranked_ICU = KenyaCoV.plot_ranked_bars_health_usage(output,scenario_tag,areanames)
    savefig(plt_ranked_HU,"reports/report"*simulation_tag*"/peak_hospital_usage_by_county"*simulation_tag*".png")
    savefig(plt_ranked_ICU,"reports/report"*simulation_tag*"/peak_ICU_usage_by_county"*simulation_tag*".png")

    plt_severe,plt_deaths = KenyaCoV.plot_ranked_bars_cases(output,scenario_tag,areanames)
    savefig(plt_severe,"reports/report"*simulation_tag*"/total_severe_cases_by_county"*simulation_tag*".png")
    savefig(plt_deaths,"reports/report"*simulation_tag*"/total_deaths_by_county"*simulation_tag*".png")

    #Total country output
    writedlm("reports/report"*simulation_tag*"/total_A_incidence_ts"*simulation_tag*".csv",output.country_incidence_A_ts.med,",")
    writedlm("reports/report"*simulation_tag*"/total_cum_A_incidence_ts"*simulation_tag*".csv",cumsum(output.country_incidence_A_ts.med),",")

    writedlm("reports/report"*simulation_tag*"/total_M_incidence_ts"*simulation_tag*".csv",output.country_incidence_M_ts.med,",")
    writedlm("reports/report"*simulation_tag*"/total_V_incidence_ts"*simulation_tag*".csv",output.country_incidence_V_ts.med,",")
    writedlm("reports/report"*simulation_tag*"/total_all_incidence_ts"*simulation_tag*".csv",output.country_incidence_A_ts.med .+ output.country_incidence_M_ts.med .+ output.country_incidence_V_ts.med,",")

    writedlm("reports/report"*simulation_tag*"/A_incidence_ts_by_county"*simulation_tag*".csv",hcat(areanames,output.incidence_A_ts.med),",")
    writedlm("reports/report"*simulation_tag*"/M_incidence_ts_by_county"*simulation_tag*".csv",hcat(areanames,output.incidence_M_ts.med),",")
    writedlm("reports/report"*simulation_tag*"/V_incidence_ts_by_county"*simulation_tag*".csv",hcat(areanames,output.incidence_V_ts.med),",")
    writedlm("reports/report"*simulation_tag*"/all_incidence_ts_by_county"*simulation_tag*".csv",hcat(areanames,output.incidence_A_ts.med .+ output.incidence_M_ts.med .+ output.incidence_V_ts.med),",")


    return scenariodata
end
