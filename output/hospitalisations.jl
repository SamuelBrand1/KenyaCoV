#In this file I prototype a simple hospital model
#Only input is number and age of each new hospitalised person on each day
#Only looking at Hosp -> ICU and ICU -> death
include("hospitalisation_data.jl");


function plot_ranked_bars_health_usage(overall_hosp_exceed,overall_ICU_exceed,chance_exceeds_hosp,chance_exceeds_ICU,scenario)

    I = sortperm(overall_hosp_exceed[:,1])
    plt_HU = bar(overall_hosp_exceed[I,1]*100,orientations = :horizonal,
                yticks = (1:47,names[I]),size = (700,550),lab ="",
                xlabel = "% Hospital available bed demand",title = "Hospital demand at peak"*scenario)
    scatter!(plt_HU,overall_hosp_exceed[I,1]*100,1:47,
                xerror = ((overall_hosp_exceed[I,1] .-overall_hosp_exceed[I,2])*100,(overall_hosp_exceed[I,3] .- overall_hosp_exceed[I,1])*100 ),ms = 0.,color = :black,lab ="")

    I = sortperm(overall_ICU_exceed[:,1])
    plt_ICU = bar(overall_ICU_exceed[I,1],orientations = :horizonal,
                yticks = (1:47,names[I]),size = (700,550),lab ="",
                xlabel = "ICU bed excess demand (numbers)",title = "ICU demand at peak"*scenario)
    scatter!(plt_ICU,overall_ICU_exceed[I,1],1:47,
                xerror = ((overall_ICU_exceed[I,1] .-overall_ICU_exceed[I,2]),(overall_ICU_exceed[I,3] .- overall_ICU_exceed[I,1]) ),ms = 0.,color = :black,lab ="")



    return plt_HU,plt_ICU
end

function get_hosp_forecast(sims)
    overall_hosp_exceed = zeros(47,3)
    overall_ICU_exceed = zeros(47,3)
    chance_exceeds_hosp = zeros(47)
    chance_exceeds_ICU = zeros(47)
    for cn = 1:47
        hosp_occup_per_sim,ICU_occup_per_sim,new_ICU_per_sim,death_incidence_per_sim = total_hospital_outcomes_per_sim(sims,cn)
        hosp_maximums = [(hosp_max,hosp_max_time) = findmax(hosp_occup_per_sim[k,:]) for k = 1:1000]
        ICU_maximums = [(ICU_max,ICU_max_time) = findmax(ICU_occup_per_sim[k,:]) for k = 1:1000]
        perc_hosp_exceeds = [h[1]/spare_capacity_H_by_county[cn] for h in hosp_maximums]
        perc_ICU_exceeds = [h[1] - spare_capacity_ICU_by_county[cn] for h in ICU_maximums]
        chance_exceeds_hosp[cn] = mean([p >= 1 for p in perc_hosp_exceeds])
        chance_exceeds_ICU[cn] = mean([p >= 1 for p in perc_ICU_exceeds])
        overall_hosp_exceed[cn,1] = median(perc_hosp_exceeds)
        overall_hosp_exceed[cn,2] = quantile(perc_hosp_exceeds,0.025)
        overall_hosp_exceed[cn,3] = quantile(perc_hosp_exceeds,0.975)
        overall_ICU_exceed[cn,1] = median(perc_ICU_exceeds)
        overall_ICU_exceed[cn,2] = quantile(perc_ICU_exceeds,0.025)
        overall_ICU_exceed[cn,3] = quantile(perc_ICU_exceeds,0.975)
    end
    return overall_hosp_exceed,overall_ICU_exceed,chance_exceeds_hosp,chance_exceeds_ICU
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

function first_introduction_time_peak_peak_value(sims,type)
    n = length(sims.u)
    T = length(sims.u[1])
    nc,na,ns = size(sims.u[1][1])
    first_times_matrix = zeros(n,nc)
    peak_times_matrix = zeros(n,nc)
    peak_value_matrix = zeros(n,nc)
    total_incidence = zeros(T-1)

    for k = 1:n
        for cn in 1:47
            cum_incidence = @view VectorOfArray(sims.u[k])[cn,:,type,:]
            for i = 1:(T-1)
                total_incidence[i] = cum_incidence[1,i+1] - cum_incidence[1,i]
                for a = 2:na
                    @inbounds total_incidence[i] += cum_incidence[a,i+1] - cum_incidence[a,i]
                end
            end
             first_time = findfirst(total_incidence .> 0)
             (peak_value,peak_time) = findmax(total_incidence)
             peak_value_matrix[k,cn] = peak_value
             peak_times_matrix[k,cn] = peak_time
             if !isnothing(first_time)
                 first_times_matrix[k,cn] = first_time
             else
                 first_times_matrix[k,cn] = 9999
             end
        end
    end
    return first_times_matrix,peak_times_matrix,peak_value_matrix
end

function get_first_time(a,b,c)
    n,nc = size(a)
    first_times  = similar(a)
    for k = 1:n,cn = 1:nc
        first_times[k,cn] = min(a[k,cn],b[k,cn],c[k,cn])
    end
    return first_times
end

function print_onset_report(first_times,peak_A,peak_A_value,peak_M,peak_M_value,peak_V,peak_V_value,names,filename)
    n,nc = size(first_times)
    df = DataFrame(county = names)
    median_time_to_peak_A = [median(peak_A[:,cn] ) for cn in 1:47]
    lb_time_to_peak_A = [quantile(peak_A[:,cn],0.025) for cn in 1:47]
    ub_time_to_peak_A = [quantile(peak_A[:,cn] ,0.975) for cn in 1:47]
    df.median_time_to_peak_asymptomatics = median_time_to_peak_A
    df.lower_estimate_time_to_peak_asymptomatics = lb_time_to_peak_A
    df.upper_estimate_time_to_peak_asymptomatics = ub_time_to_peak_A

    median_time_to_peak_M = [median(peak_M[:,cn] ) for cn in 1:47]
    lb_time_to_peak_M = [quantile(peak_M[:,cn] ,0.025) for cn in 1:47]
    ub_time_to_peak_M = [quantile(peak_M[:,cn] ,0.975) for cn in 1:47]
    df.median_time_to_peak_Mild_cases = median_time_to_peak_M
    df.lower_estimate_time_to_peak_Mild_cases = lb_time_to_peak_M
    df.upper_estimate_time_to_peak_Mild_cases = ub_time_to_peak_M

    median_time_to_peak_V = [median(peak_V[:,cn] ) for cn in 1:47]
    lb_time_to_peak_V = [quantile(peak_V[:,cn] ,0.025) for cn in 1:47]
    ub_time_to_peak_V = [quantile(peak_V[:,cn] ,0.975) for cn in 1:47]
    df.median_time_to_peak_Severe_cases = median_time_to_peak_V
    df.lower_estimate_time_to_peak_Severe_cases  = lb_time_to_peak_V
    df.upper_estimate_time_to_peak_Severe_cases  = ub_time_to_peak_V

    median_peak_size_A = [median(peak_A_value[:,cn] ) for cn in 1:47]
    lb_peak_size_A = [quantile(peak_A_value[:,cn] ,0.025) for cn in 1:47]
    ub_peak_size_A = [quantile(peak_A_value[:,cn] ,0.975) for cn in 1:47]
    df.median_peak_size_A = median_peak_size_A
    df.lower_estimate_peak_size_A  = lb_peak_size_A
    df.upper_estimate_peak_size_A  = ub_peak_size_A

    median_peak_size_M = [median(peak_M_value[:,cn] ) for cn in 1:47]
    lb_peak_size_M = [quantile(peak_M_value[:,cn] ,0.025) for cn in 1:47]
    ub_peak_size_M = [quantile(peak_M_value[:,cn] ,0.975) for cn in 1:47]
    df.median_peak_size_Mild_cases = median_peak_size_M
    df.lower_estimate_peak_size_Mild_cases  = lb_peak_size_M
    df.upper_estimate_peak_size_Mild_cases  = ub_peak_size_M

    median_peak_size_V = [median(peak_V_value[:,cn] ) for cn in 1:47]
    lb_peak_size_V = [quantile(peak_V_value[:,cn] ,0.025) for cn in 1:47]
    ub_peak_size_V = [quantile(peak_V_value[:,cn] ,0.975) for cn in 1:47]
    df.median_peak_size_Severe_cases = median_peak_size_V
    df.lower_estimate_peak_size_Severe_cases  = lb_peak_size_V
    df.upper_estimate_peak_size_Severe_cases  = ub_peak_size_V

    CSV.write(filename,df)
    return df
end


function generate_checklist(n)
    checklist = Any[]
    push!(checklist,1:n)
    for i = 1:n
        push!(checklist,i)
    end
    return checklist
end

function total_incidence_and_cum_incidence_for_each_sim(sims)
    n = length(sims.u)
    T = length(sims.u[1])
    cum_incidence_for_each_sim = zeros(n,T)
    for k = 1:n,t = 1:T
        cum_incidence_for_each_sim[k,t] = sum(sims.u[k][t][:,:,1:3])
    end
    return cum_incidence_for_each_sim,diff(cum_incidence_for_each_sim,dims = 2)
end

function total_incidence_and_cum_incidence_for_each_sim_for_a_county(sims,cn)
    n = length(sims.u)
    T = length(sims.u[1])
    cum_incidence_for_each_sim = zeros(n,T)
    for k = 1:n,t = 1:T
        cum_incidence_for_each_sim[k,t] = sum(sims.u[k][t][cn,:,1:3])
    end
    return cum_incidence_for_each_sim,diff(cum_incidence_for_each_sim,dims = 2)
end

function find_peak_time_by_sim(sims)
    n = length(sims.u)
    T = length(sims.u[1])
    cum_incidence_for_each_sim = zeros(n,T)
    for k = 1:n,t = 1:T
        cum_incidence_for_each_sim[k,t] = sum(sims.u[k][t][:,:,1:3])
    end
    incidence = diff(cum_incidence_for_each_sim,dims = 2)
    return [findmax(incidence[k,:]) for k = 1:n]
end

function find_peaks(incidence_array)
    n,T = size(incidence_array)
    return [findmax(incidence_array[k,:])[1] for k = 1:n]
end

function cum_incidence_for_each_sim_by_type(sims)
    n = length(sims.u)
    T = length(sims.u[1])
    cum_incidence_A = zeros(n,T)
    cum_incidence_M = zeros(n,T)
    cum_incidence_V = zeros(n,T)
    cum_incidence_H = zeros(n,T)

    for k = 1:n,t = 1:T
        cum_incidence_A[k,t] = sum(sims.u[k][t][:,:,1])
        cum_incidence_M[k,t] = sum(sims.u[k][t][:,:,2])
        cum_incidence_V[k,t] = sum(sims.u[k][t][:,:,3])
        cum_incidence_H[k,t] = sum(sims.u[k][t][:,:,4])
    end
    return cum_incidence_A,cum_incidence_M,cum_incidence_V,cum_incidence_H
end

function cum_incidence_for_each_sim_by_type(sims,cn)
    n = length(sims.u)
    T = length(sims.u[1])
    cum_incidence_A = zeros(n,T)
    cum_incidence_M = zeros(n,T)
    cum_incidence_V = zeros(n,T)
    cum_incidence_H = zeros(n,T)

    for k = 1:n,t = 1:T
        cum_incidence_A[k,t] = sum(sims.u[k][t][cn,:,1])
        cum_incidence_M[k,t] = sum(sims.u[k][t][cn,:,2])
        cum_incidence_V[k,t] = sum(sims.u[k][t][cn,:,3])
        cum_incidence_H[k,t] = sum(sims.u[k][t][cn,:,4])
    end
    return cum_incidence_A,cum_incidence_M,cum_incidence_V,cum_incidence_H
end
"""
generate_hospitalisation_outcomes(H)
"""
death_rate_critical_cases = mean([0.5,0.75])

function split_incoming_hospitalisations_into_outcomes(H,age,area)
    T = length(H)
    #randomly split population into hospital and ICU
    new_hosp = zeros(Int64,T);new_ICU = zeros(Int64,T); new_to_die = zeros(Int64,T);
    for t = 1:T
        outcomes = rand(Multinomial(Int64.(H[t]),[1 - ICU_rate_by_age_cond_hosp[age],
                                            ICU_rate_by_age_cond_hosp[age]*(1-death_rate_critical_cases),
                                            ICU_rate_by_age_cond_hosp[age]*death_rate_critical_cases]) )
        new_hosp[t] = outcomes[1]
        new_ICU[t] = outcomes[2]
        new_to_die[t] = outcomes[3]
    end
    return new_hosp,new_ICU,new_to_die
end

function generate_hospitalisation_outcomes(H,age,area)
    T = length(H)
    hosp_occup = zeros(Int64,T)
    ICU_occup = zeros(Int64,T)
    death_incidence = zeros(Int64,T)
    #randomly split population into hospital and ICU
    new_hosp,new_ICU,new_to_die = split_incoming_hospitalisations_into_outcomes(H,age,area)
    #Assign population outcomes
    for t = 1:T
        for k = 1:new_hosp[t]
            tₑ = rand(6:10)
            timesinhosp = t:(min(T,t+tₑ))#6-10 days in Hospital then discharge
            hosp_occup[timesinhosp] .+= 1
        end
        for k = 1:new_ICU[t]
            tₑ = rand(7:14)
            timesinICU = t:(min(T,t+tₑ))#7-14 days in ICU then discharge
            ICU_occup[timesinICU] .+= 1
        end
        for k = 1:new_to_die[t]
            tₑ = rand(7:14)
            timesinICU = t:(min(T,t+tₑ))#7-14 days in ICU then death
            ICU_occup[timesinICU] .+= 1
            death_incidence[min(T,t+tₑ)] += 1
        end
    end

    return hosp_occup, ICU_occup, new_ICU,death_incidence
end



function total_hospital_outcomes_per_sim(sims)
    n = length(sims.u)
    T = length(sims.u[1])
    hosp_occup_per_sim = zeros(Int64,n,T-1)
    ICU_occup_per_sim = zeros(Int64,n,T-1)
    new_ICU_per_sim = zeros(Int64,n,T-1)
    death_incidence_per_sim = zeros(Int64,n,T-1)
    for k = 1:n
        for a = 1:17
            H = diff([sum(sims.u[k][t][:,a,4]) for t = 1:659])
            hosp_occup, ICU_occup, new_ICU,death_incidence = generate_hospitalisation_outcomes(H,a,1)
            hosp_occup_per_sim[k,:] .+= hosp_occup
            ICU_occup_per_sim[k,:] .+= ICU_occup
            new_ICU_per_sim[k,:] .+= new_ICU
            death_incidence_per_sim[k,:] .+= death_incidence
        end
    end
    return hosp_occup_per_sim,ICU_occup_per_sim,new_ICU_per_sim,death_incidence_per_sim
end

function total_hospital_outcomes_per_sim(sims,cn)
    n = length(sims.u)
    T = length(sims.u[1])
    hosp_occup_per_sim = zeros(Int64,n,T-1)
    ICU_occup_per_sim = zeros(Int64,n,T-1)
    new_ICU_per_sim = zeros(Int64,n,T-1)
    death_incidence_per_sim = zeros(Int64,n,T-1)
    for k = 1:n
        for a = 1:17
            H = diff([sum(sims.u[k][t][cn,a,4]) for t = 1:659])
            hosp_occup, ICU_occup, new_ICU,death_incidence = generate_hospitalisation_outcomes(H,a,1)
            hosp_occup_per_sim[k,:] .+= hosp_occup
            ICU_occup_per_sim[k,:] .+= ICU_occup
            new_ICU_per_sim[k,:] .+= new_ICU
            death_incidence_per_sim[k,:] .+= death_incidence
        end
    end
    return hosp_occup_per_sim,ICU_occup_per_sim,new_ICU_per_sim,death_incidence_per_sim
end

function give_plots_for_county(sims,cn)
    cum_A,cum_M,cum_V,cum_H = cum_incidence_for_each_sim_by_type(sims,cn)
    _hosp_occup_per_sim,_ICU_occup_per_sim,_new_ICU_per_sim,_death_incidence_per_sim = total_hospital_outcomes_per_sim(sims,cn)
    _incidence_H = diff(cum_H,dims=2)
    median__incidence_H = [median(_incidence_H[:,t]) for t in 1:658]
    lb__incidence_H = [quantile(_incidence_H[:,t],0.025) for t in 1:658]
    ub__incidence_H = [quantile(_incidence_H[:,t],0.975) for t in 1:658]
    median__incidence_death = [median(_death_incidence_per_sim[:,t]) for t in 1:658]
    lb__incidence_death = [quantile(_death_incidence_per_sim[:,t],0.025) for t in 1:658]
    ub__incidence_death = [quantile(_death_incidence_per_sim[:,t],0.975) for t in 1:658]
    median__H_occ = [median(_hosp_occup_per_sim[:,t]) for t in 1:658]
    lb__H_occ = [quantile(_hosp_occup_per_sim[:,t],0.025) for t in 1:658]
    ub__H_occ = [quantile(_hosp_occup_per_sim[:,t],0.975) for t in 1:658]
    median__ICU_occ = [median(_ICU_occup_per_sim[:,t]) for t in 1:658]
    lb__ICU_occ = [quantile(_ICU_occup_per_sim[:,t],0.025) for t in 1:658]
    ub__ICU_occ = [quantile(_ICU_occup_per_sim[:,t],0.975) for t in 1:658]

    plt_incidence_ = plot(0:657,median__incidence_H,
                        lab = "Hospitalisations",
                        lw = 3,color = :red,
                        ribbon = (median__incidence_H.-lb__incidence_H,ub__incidence_H .- median__incidence_H),
                        fillalpha = 0.25,
                        xticks = (first_of_months,["Apr","May","June","July","Aug","Sept","Oct","Nov","Dec","Jan","Feb","Mar","Apr","May","June","July","Aug","Sept","Oct","Nov","Dec"]),
                        xlims = (0.,657.),
                        title = "Incidence of disease - cn $(cn)")
    plot!(plt_incidence_,0:657,median__incidence_death,
                        lab = "Deaths",
                        lw = 3,color = :black,
                        ribbon = (median__incidence_death.-lb__incidence_death,ub__incidence_death .- median__incidence_death),
                        fillalpha = 0.5,
                        xticks = (first_of_months,["Apr","May","June","July","Aug","Sept","Oct","Nov","Dec","Jan","Feb","Mar","Apr","May","June","July","Aug","Sept","Oct","Nov","Dec"]),
                        xlims = (0.,657.),
                        ylabel = "Daily Incidence" )

    plt_health_usage = plot(0:657,median__H_occ,
                            lab = "hospital beds",
                            lw = 3,color = :blue,
                            ribbon = (median__H_occ.-lb__H_occ,ub__H_occ .- median__H_occ),
                            fillalpha = 0.25,
                            xticks = (first_of_months,["Apr","May","June","July","Aug","Sept","Oct","Nov","Dec","Jan","Feb","Mar","Apr","May","June","July","Aug","Sept","Oct","Nov","Dec"]),
                            xlims = (0.,657.),
                            ylabel = "Daily Occupancy",
                            title = "Health system usage")

    plot!(plt_health_usage,0:657,median__ICU_occ,
                            lab = "ICU beds",
                            lw = 3,color = :green,
                            ribbon = (median__ICU_occ.-lb__ICU_occ,ub__ICU_occ .- median__ICU_occ),
                            fillalpha = 0.25,
                            xticks = (first_of_months,["Apr","May","June","July","Aug","Sept","Oct","Nov","Dec","Jan","Feb","Mar","Apr","May","June","July","Aug","Sept","Oct","Nov","Dec"]),
                            xlims = (0.,657.),
                            ylabel = "Daily Occupancy",
                            title = "Health system usage")
            return plt_incidence_,plt_health_usage


end
