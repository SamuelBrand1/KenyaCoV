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
