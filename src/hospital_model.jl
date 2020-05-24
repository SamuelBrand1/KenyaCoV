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
