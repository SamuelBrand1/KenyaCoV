


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

    peak_hosp_by_area_by_sim = zeros(nc,n)
    peak_ICU_by_area_by_sim = zeros(nc,n)

    for cn = 1:nc,k = 1:n
        peak_hosp_by_area_by_sim[cn,k] = maximum(sims[k].hosp_occup_by_area_ts[cn,:])
        peak_ICU_by_area_by_sim[cn,k] = maximum(sims[k].ICU_occup_by_area_ts[cn,:])
    end


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
            ICU_occup_by_area_over_sims=ICU_occup_by_area_over_sims)
end

function give_plots_for_county(output,areas_to_plot,scenario_tag,areanames)
    incidence_H = sum(output.incidence_H_by_area_over_sims[areas_to_plot,:,:],dims=1)[1,:,:]
    incidence_D = sum(output.incidence_death_by_area_over_sims[areas_to_plot,:,:],dims=1)[1,:,:]
    H_occup = sum(output.H_occup_by_area_over_sims[areas_to_plot,:,:],dims=1)[1,:,:]
    ICU_occup = sum(output.ICU_occup_by_area_over_sims[areas_to_plot,:,:],dims=1)[1,:,:]

    T,n = size(incidence_H)
    monthdates = [Date(2020,3,1) + Month(i) for i = 1:21 ]
    monthnames = [monthname(d)[1:3]*"-$(year(d)-2000)" for d in monthdates]
    tick_times = [(d - Date(2020,3,13)).value for d in monthdates]

    median_H = [median(incidence_H[t,:]) for t = 1:T]
    lb_H = [quantile(incidence_H[t,:],0.025) for t = 1:T]
    ub_H = [quantile(incidence_H[t,:],0.975) for t = 1:T]
    median_D = [median(incidence_D[t,:]) for t = 1:T]
    lb_D = [quantile(incidence_D[t,:],0.025) for t = 1:T]
    ub_D = [quantile(incidence_D[t,:],0.975) for t = 1:T]
    median_H_occup = [median(H_occup[t,:]) for t = 1:T]
    lb_H_occup = [quantile(H_occup[t,:],0.025) for t = 1:T]
    ub_H_occup = [quantile(H_occup[t,:],0.975) for t = 1:T]
    median_ICU_occup = [median(ICU_occup[t,:]) for t = 1:T]
    lb_ICU_occup = [quantile(ICU_occup[t,:],0.025) for t = 1:T]
    ub_ICU_occup = [quantile(ICU_occup[t,:],0.975) for t = 1:T]

    area_name = ["Rest of Kenya"]
    if length(areas_to_plot) == 1
        area_name = areanames[areas_to_plot]
    end

    plt_incidence = plot(0:(T-1),median_H,
                        lab = "Hospitalisations",
                        lw = 3,color = :red,
                        ribbon = (median_H.-lb_H,ub_H .- median_H),
                        fillalpha = 0.25,
                        xticks = (tick_times,monthnames),
                        title = "Incidence of disease - "*area_name[1]*scenario_tag)
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
                            title = "Health system usage - "*area_name[1]*scenario_tag)

    plot!(plt_health_usage,0:(T-1),median_ICU_occup,
                            lab = "ICU beds",
                            lw = 3,color = :green,
                            ribbon = (median_ICU_occup.-lb_ICU_occup,ub_ICU_occup .- median_ICU_occup),
                            fillalpha = 0.25)
    

    return plt_incidence,plt_health_usage
end




function plot_ranked_bars_health_usage(output,scenario_tag,areanames)
    median_hospital_exceed = output.hosp_peak_excess_demand_by_area.med
    median_ICU_exceed = output.ICU_peak_excess_demand_by_area.med
    nc = length(areanames)

    I = sortperm(median_hospital_exceed)
    plt_HU = bar(median_hospital_exceed[I]*100,orientations = :horizonal,
                yticks = (1:47,areanames[I]),size = (700,550),lab ="",
                xlabel = "% Hospital available bed demand",title = "Hospital demand at peak"*scenario_tag)
    scatter!(plt_HU,median_hospital_exceed[I]*100,1:nc,
                xerror = ((median_hospital_exceed[I] .- output.hosp_peak_excess_demand_by_area.lpred[I])*100
                            ,(output.hosp_peak_excess_demand_by_area.upred[I] .- median_hospital_exceed[I])*100 ),
                ms = 0.,color = :black,lab ="")

    I = sortperm(median_ICU_exceed)
    plt_ICU = bar(median_ICU_exceed[I,1],orientations = :horizonal,
                yticks = (1:nc,areanames[I]),size = (700,550),lab ="",
                xlabel = "ICU bed excess demand (numbers)",title = "ICU demand at peak"*scenario_tag)
    scatter!(plt_ICU,median_ICU_exceed[I],1:nc,
                xerror = ((median_ICU_exceed[I] .- output.ICU_peak_excess_demand_by_area.lpred[I]),
                            (output.ICU_peak_excess_demand_by_area.upred[I] .- median_ICU_exceed[I]) ),
                ms = 0.,color = :black,lab ="")



    return plt_HU,plt_ICU
end


function plot_ranked_bars_cases(output,scenario_tag,areanames)
    median_severe_cases = output.severe_cases_by_area.med
    median_deaths = output.deaths_by_area.med
    nc = length(areanames)
    I = sortperm(median_severe_cases)
    plt = bar(median_severe_cases[I],orientations = :horizonal,
                yticks = (1:nc,areanames[I]),size = (700,600),lab ="",
                xlabel = "Total severe cases",title = "Severe cases by county"*scenario_tag)
    scatter!(plt,median_severe_cases[I],1:nc,
            xerror = (median_severe_cases[I] .- output.severe_cases_by_area.lpred[I],
                        output.severe_cases_by_area.upred[I] .- median_severe_cases[I] ),
             ms = 0.,color = :black,lab ="")

    I = sortperm(median_deaths)

    plt_deaths = bar(median_deaths[I],orientations = :horizonal,
                yticks = (1:nc,areanames[I]),size = (700,600),lab ="",
                xlabel = "Total deaths",title = "Deaths by county"*scenario_tag)
    scatter!(plt_deaths,median_deaths[I],1:nc,
            xerror = (median_deaths[I] .- output.deaths_by_area.lpred[I],
                        output.deaths_by_area.upred[I] .- median_deaths[I] ),
             ms = 0.,color = :black,lab ="")

    return plt,plt_deaths
end
