#In this file I prototype a simple hospital model
#Only input is number and age of each new hospitalised person on each day
#Only looking at Hosp -> ICU and ICU -> death
include("hospitalisation_data.jl");
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
            H = diff([sum(sims.u[k][t][:,a,4]) for t = 1:366])
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
            H = diff([sum(sims.u[k][t][cn,a,4]) for t = 1:366])
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
