
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

    hosp_by_area_over_sims =  VectorOfArray([sum(sims[k].total_hosp_by_area_and_age,dims=2) for k = 1:n])[:,1,:]
    nc,na,T
end

function plot_ranked_bars_cases(sims,scenario)
    median_final_number = [median([sum(sims.u[k][end][cn,:,3]) for k = 1:1000] ) for cn = 1:47]
    lb_final_number = [quantile([sum(sims.u[k][end][cn,:,3]) for k = 1:1000],0.025) for cn = 1:47]
    ub_final_number = [quantile([sum(sims.u[k][end][cn,:,3]) for k = 1:1000],0.975) for cn = 1:47]
    deaths_matrix = zeros(1000,47)
    for k = 1:1000,cn = 1:47
        deaths_matrix[k,cn] = sum(sims.u[k][end][cn,:,4].*ICU_rate_by_age_cond_hosp.*0.625)
    end
    median_deaths = [median(deaths_matrix[:,cn] ) for cn = 1:47]
    lb_deaths = [quantile(deaths_matrix[:,cn],0.025) for cn = 1:47]
    ub_deaths = [quantile(deaths_matrix[:,cn],0.975) for cn = 1:47]


    I = sortperm(median_final_number)
    plt = bar(median_final_number[I],orientations = :horizonal,
                yticks = (1:47,names[I]),size = (700,600),lab ="",
                xlabel = "Total severe cases",title = "Severe cases by county"*scenario)
    scatter!(plt,median_final_number[I],1:47,
                xerror = (median_final_number[I] .-lb_final_number[I],ub_final_number[I] .- median_final_number[I] ),ms = 0.,color = :black,lab ="")

    I = sortperm(median_deaths)
    plt_deaths = bar(median_deaths[I],orientations = :horizonal,
                yticks = (1:47,names[I]),size = (700,500),lab ="",
                xlabel = "Total deaths",title = "Deaths by county"*scenario)
    scatter!(plt_deaths,median_deaths[I],1:47,
                xerror = (median_deaths[I] .-lb_deaths[I],ub_deaths[I] .- median_deaths[I] ),ms = 0.,color = :black,lab ="")

    return plt,plt_deaths
end
