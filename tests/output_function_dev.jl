
"""
function output_simulation_data
    This function calculates daily incidence (asymptomatic, mild and severe) from the time-stepping of the underlying dynamic model.
    It also runs the hospital model as a post-processing layer after simulation
"""
function output_simulation_data(sol,i)
    times = sol.prob.tspan[1]:1:sol.prob.tspan[end]
    incidence_A = diff([sum(sol(t)[:,:,9],dims=2)[:]  for t in times])
    incidence_M = diff([sum(sol(t)[:,:,10],dims=2)[:]  for t in times])
    incidence_V = diff([sum(sol(t)[:,:,11],dims=2)[:]  for t in times])
    incidence_H = diff([sol(t)[:,:,12]  for t in times])
    T = length(incidence_H)
    nc,na = size(incidence_H[1])
    total_hosp_occup = zeros(nc,T)
    total_ICU_occup = zeros(nc,T)
    total_new_ICU = zeros(nc,T)
    total_death_incidence = zeros(nc,T)
    hosp_by_area_and_age = zeros(nc,na)
    ICU_by_area_and_age = zeros(nc,na)
    deaths_by_area_and_age = zeros(nc,na)

    for cn in 1:nc, a in 1:na
        hosp_by_area_and_age[cn,a] = sum(sol(sol.prob.tspan[end])[cn,a,12])
    end

    for cn in 1:nc, a in 1:na
        hosp_occup, ICU_occup, new_ICU,death_incidence = KenyaCoV.generate_hospitalisation_outcomes([inc_h[cn,a] for inc_h in incidence_H],a,cn)
        total_hosp_occup[cn,:] .+= hosp_occup
        total_ICU_occup[cn,:] .+= ICU_occup
        total_new_ICU[cn,:] .+= new_ICU
        total_death_incidence[cn,:] .+= death_incidence
        ICU_by_area_and_age[cn,a] = sum(new_ICU)
        deaths_by_area_and_age[cn,a] = sum(death_incidence)
    end

    return (incidence_A=VectorOfArray(incidence_A)[:,:],
            incidence_M=VectorOfArray(incidence_M)[:,:],
            incidence_V=VectorOfArray(incidence_V)[:,:],
            incidence_H=VectorOfArray(incidence_H)[:,:],
            hosp_occup_by_area_ts = total_hosp_occup,
            ICU_occup_by_area_ts = total_ICU_occup,
            incidence_ICU_by_area_ts = total_new_ICU,
            death_incidence_by_area_ts = total_death_incidence,
            total_hosp_by_area_and_age = hosp_by_area_and_age,
            total_ICU_by_area_and_age = ICU_by_area_and_age,
            total_deaths_by_area_and_age = deaths_by_area_and_age),false
end

function extract_information_from_simulations(sims)
    n = length(sims)
    nc,T = size(sims[1].hosp_occup_by_area_ts)
    nc,na = size(sims[1].total_ICU_by_area_and_age)

    hosp_by_area_over_sims = VectorOfArray([sum(sims[k].total_hosp_by_area_and_age,dims=2) for k = 1:n])[:,1,:]
    deaths_by_area_over_sims = VectorOfArray([sum(sims[k].total_deaths_by_area_and_age,dims=2) for k = 1:n])[:,1,:]
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
    total_deaths = (med = median(sum(deaths_by_area_over_sims,dims=1)),
                            lpred = quantile(sum(deaths_by_area_over_sims,dims=1)[:],0.025),
                            upred = quantile(sum(deaths_by_area_over_sims[:],dims=1),0.975))
    deaths_by_area = (med = [median(deaths_by_area_over_sims[cn,:]) for cn = 1:nc],
                            lpred = [quantile(deaths_by_area_over_sims[cn,:],0.025) for cn = 1:nc],
                            upred = [quantile(deaths_by_area_over_sims[cn,:],0.975) for cn = 1:nc])

    hosp_peak_excess_demand_by_area = (med = [median(peak_hosp_by_area_by_sim[cn,:]/KenyaCoV.spare_capacity_H_by_county[cn]) for cn = 1:nc],
                            lpred = [quantile(peak_hosp_by_area_by_sim[cn,:]/KenyaCoV.spare_capacity_H_by_county[cn],0.025) for cn = 1:nc],
                            upred = [quantile(peak_hosp_by_area_by_sim[cn,:]/KenyaCoV.spare_capacity_H_by_county[cn],0.975) for cn = 1:nc])

    ICU_peak_excess_demand_by_area = (med = [median(peak_ICU_by_area_by_sim[cn,:] .- KenyaCoV.spare_capacity_ICU_by_county[cn]) for cn = 1:nc],
                            lpred = [quantile(peak_ICU_by_area_by_sim[cn,:] .- KenyaCoV.spare_capacity_ICU_by_county[cn],0.025) for cn = 1:nc],
                            upred = [quantile(peak_ICU_by_area_by_sim[cn,:] .- KenyaCoV.spare_capacity_ICU_by_county[cn],0.975) for cn = 1:nc])

    return (total_severe_cases=total_severe_cases,
            severe_cases_by_area=severe_cases_by_area,
            total_deaths=total_deaths,
            deaths_by_area=deaths_by_area,
            hosp_peak_excess_demand_by_area=hosp_peak_excess_demand_by_area,
            ICU_peak_excess_demand_by_area=ICU_peak_excess_demand_by_area)
end

function give_plots_for_county(output,areas_to_plot,scenario_tag)
    cum_A,cum_M,cum_V,cum_H = cum_incidence_for_each_sim_by_type(sims,cn)
    _hosp_occup_per_sim,_ICU_occup_per_sim,_new_ICU_per_sim,_death_incidence_per_sim = total_hospital_outcomes_per_sim(sims,cn)
    _incidence_H = diff(cum_H,dims=2)
    median__incidence_H = [median(_incidence_H[:,t]) for t in 1:365]
    lb__incidence_H = [quantile(_incidence_H[:,t],0.025) for t in 1:365]
    ub__incidence_H = [quantile(_incidence_H[:,t],0.975) for t in 1:365]
    median__incidence_death = [median(_death_incidence_per_sim[:,t]) for t in 1:365]
    lb__incidence_death = [quantile(_death_incidence_per_sim[:,t],0.025) for t in 1:365]
    ub__incidence_death = [quantile(_death_incidence_per_sim[:,t],0.975) for t in 1:365]
    median__H_occ = [median(_hosp_occup_per_sim[:,t]) for t in 1:365]
    lb__H_occ = [quantile(_hosp_occup_per_sim[:,t],0.025) for t in 1:365]
    ub__H_occ = [quantile(_hosp_occup_per_sim[:,t],0.975) for t in 1:365]
    median__ICU_occ = [median(_ICU_occup_per_sim[:,t]) for t in 1:365]
    lb__ICU_occ = [quantile(_ICU_occup_per_sim[:,t],0.025) for t in 1:365]
    ub__ICU_occ = [quantile(_ICU_occup_per_sim[:,t],0.975) for t in 1:365]

    plt_incidence_ = plot(0:364,median__incidence_H,
                        lab = "Hospitalisations",
                        lw = 3,color = :red,
                        ribbon = (median__incidence_H.-lb__incidence_H,ub__incidence_H .- median__incidence_H),
                        fillalpha = 0.25,
                        xticks = (first_of_months,["Apr","May","June","July","Aug","Sept","Oct","Nov","Dec"]),
                        xlims = (0.,275.),
                        title = "Incidence of disease - cn $(cn)")
    plot!(plt_incidence_,0:364,median__incidence_death,
                        lab = "Deaths",
                        lw = 3,color = :black,
                        ribbon = (median__incidence_death.-lb__incidence_death,ub__incidence_death .- median__incidence_death),
                        fillalpha = 0.5,
                        xticks = (first_of_months,["Apr","May","June","July","Aug","Sept","Oct","Nov","Dec"]),
                        xlims = (0.,275.),
                        ylabel = "Daily Incidence" )

    plt_health_usage = plot(0:364,median__H_occ,
                            lab = "hospital beds",
                            lw = 3,color = :blue,
                            ribbon = (median__H_occ.-lb__H_occ,ub__H_occ .- median__H_occ),
                            fillalpha = 0.25,
                            xticks = (first_of_months,["Apr","May","June","July","Aug","Sept","Oct","Nov","Dec"]),
                            xlims = (0.,275.),
                            ylabel = "Daily Occupancy",
                            title = "Health system usage")

    plot!(plt_health_usage,0:364,median__ICU_occ,
                            lab = "ICU beds",
                            lw = 3,color = :green,
                            ribbon = (median__ICU_occ.-lb__ICU_occ,ub__ICU_occ .- median__ICU_occ),
                            fillalpha = 0.25,
                            xticks = (first_of_months,["Apr","May","June","July","Aug","Sept","Oct","Nov","Dec"]),
                            xlims = (0.,275.),
                            ylabel = "Daily Occupancy",
                            title = "Health system usage")
            return plt_incidence_,plt_health_usage


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
