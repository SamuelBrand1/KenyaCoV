"""
Functions for analysing simulations

Output functions: Target the peak timing for each county, cases by each county,
timing of total peak and total cases
"""

"""
Define uncertainty of parameter estimates
"""
d_incubation = LogNormal(log(5.),0.25) #Liu et al
mean(d_incubation)
(quantile(d_incubation,0.025),median(d_incubation),quantile(d_incubation,0.975))
d_R₀ = Gamma(100,2.5/100) ##Liu et al
mean(d_R₀)
(quantile(d_R₀,0.025),median(d_R₀),quantile(d_R₀,0.975))

########
"""
Simulation functions
"""

#=function consensus_randomise_params_β25(prob,i,repeat) #Remember to rescale susceptibility by the inverse leading eigenvalue
    _P = deepcopy(prob.p)
    #_P.β = rand(d_R₀)
    return remake(prob,p=_P)
end=#
function randomise_R₀(prob,i,repeat) #Remember to rescale susceptibility by the inverse leading eigenvalue
    _P = deepcopy(prob.p)
    _P.β = rand(d_R₀)
    return remake(prob,p=_P)
end
function output_daily_incidence_and_hosp_cums(sol,i)
    times = sol.prob.tspan[1]:1:sol.prob.tspan[end]
    cumA=sol[end][:,:,9]
    cumM=sol[end][:,:,10]
    cumV=sol[end][:,:,11]
    cumH=sol[end][:,:,12]
    cumQ=sol[end][:,:,16]
    I = [sum(sol(t)[:,:,4:6])  for t in times]
    return (cumA,cumM,cumV,cumH,cumQ,I),false
end
function output_simulation_data(sol,i)
    times = sol.prob.tspan[1]:1:sol.prob.tspan[end]

    total_cases_A_by_area_and_age=sol(sol.prob.tspan[end])[:,:,9]   ##
    total_cases_M_by_area_and_age=sol(sol.prob.tspan[end])[:,:,10]  ##
    total_cases_V_by_area_and_age=sol(sol.prob.tspan[end])[:,:,11]  ##
    total_q_by_area_and_age=sol(sol.prob.tspan[end])[:,:,16]        ##
    #total_cases=sum(sol(sol.prob.tspan[end])[:,:,9])+sum(sol(sol.prob.tspan[end])[:,:,10])+sum(sol(sol.prob.tspan[end])[:,:,11])  ###
    #total_q=sum(sol(sol.prob.tspan[end])[:,:,16])                   ###

    incidence_A = diff([sum(sol(t)[:,:,9],dims=2)[:]  for t in times])
    incidence_M = diff([sum(sol(t)[:,:,10],dims=2)[:]  for t in times])
    incidence_V = diff([sum(sol(t)[:,:,11],dims=2)[:]  for t in times])
    incidence_H_by_area = diff([sum(sol(t)[:,:,12],dims=2)[:]  for t in times])
    incidence_H_by_area_and_age = diff([sol(t)[:,:,12]  for t in times])
    T = length(incidence_H_by_area_and_age)
    nc,na = size(incidence_H_by_area_and_age[1])
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
        hosp_occup, ICU_occup, new_ICU,death_incidence = generate_hospitalisation_outcomes([inc_h[cn,a] for inc_h in incidence_H_by_area_and_age],a,cn)
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
            incidence_H=VectorOfArray(incidence_H_by_area)[:,:],
            hosp_occup_by_area_ts = total_hosp_occup,
            ICU_occup_by_area_ts = total_ICU_occup,
            incidence_ICU_by_area_ts = total_new_ICU,
            death_incidence_by_area_ts = total_death_incidence,
            total_hosp_by_area_and_age = hosp_by_area_and_age,
            total_ICU_by_area_and_age = ICU_by_area_and_age,
            total_deaths_by_area_and_age = deaths_by_area_and_age,
            total_cases_A_by_area_and_age=total_cases_A_by_area_and_age,        ##
            total_cases_M_by_area_and_age=total_cases_M_by_area_and_age,        ##
            total_cases_V_by_area_and_age=total_cases_V_by_area_and_age,         ##
            total_q_by_area_and_age=total_q_by_area_and_age),false              ##
end

#=function run_simulations_β25(P::CoVParameters_Screening,prob,n_traj,cb)
    ensemble_prob = EnsembleProblem(prob,
                                    prob_func = consensus_randomise_params_β25,
                                    output_func = output_daily_incidence_and_hosp_cums)
    return solve(ensemble_prob,FunctionMap(),dt = P.dt,callback = cb,trajectories = n_traj)
end
function run_simulations_βrand(P::CoVParameters_Screening,prob,n_traj,cb)
    ensemble_prob = EnsembleProblem(prob,
                                    prob_func = consensus_randomise_params_βrand,
                                    output_func = output_daily_incidence_and_hosp_cums)
    return solve(ensemble_prob,FunctionMap(),dt = P.dt,callback = cb,trajectories = n_traj)
end=#

function run_simulations_βrand_withH(P::CoVParameters_Screening,prob,n_traj;interventions = CallbackSet())
    ensemble_prob = EnsembleProblem(prob,
                                    prob_func = randomise_R₀,
                                    output_func = output_simulation_data)
    return solve(ensemble_prob,FunctionMap(),dt = P.dt,callback = interventions,trajectories = n_traj)
end

######## Function to initialize dynamics parameters
interventions=[#=1=#"Contact tracing of hospitalized only",#=2=#"Symptomatic screening",#=3=#"Symptomatic screening with contact tracing",
                #=4=#"Mass screening",#=5=#"Mass screening with contact tracing"]

function run_scenarios()
    @load "data/detection_rates_for_different_taus.jld2" d_0 d_01 d_025 d_05 d_1
    #u0,P,P_dest = model_ingredients_from_data_screening("data/data_for_age_structuredmodel.jld2","data/flight_numbers.csv","data/projected_global_prevelance.csv")
    u0,P,P_dest = model_ingredients_from_data_screening("data/data_for_age_structuredmodel_with_counties.jld2","data/flight_numbers.csv","data/projected_global_prevelance.csv")
        @load "data/detection_rates_for_different_taus.jld2" d_0 d_01 d_025 d_05 d_1
        #Redistribute susceptibility JUST to rescale infectiousness so we get the correct R₀/r
        P.χ = ones(n_a)
        P.rel_detection_rate = d_1
        P.dt = 0.25
        P.ext_inf_rate = 0.
        P.ϵ = 1.
        #Set the susceptibility vector --- just to specify the correct R₀
        sus_matrix = repeat(P.χ,1,17)
        R_A = P.ϵ*((1/P.σ₂) + (1/P.γ) ) #effective duration of asymptomatic
        R_M = (P.ϵ/P.σ₂) + (P.ϵ_D/P.γ) #effective duration of mild
        R_V = (P.ϵ/P.σ₂) + (P.ϵ_V/P.τ) #effective duration of severe
        R_vector = [(1-P.rel_detection_rate[a])*R_A + P.rel_detection_rate[a]*(1-P.hₐ[a])*R_M + P.rel_detection_rate[a]*P.hₐ[a]*R_V for a = 1:17]
        inf_matrix = repeat(R_vector',17,1)
        eigs, = eigen(sus_matrix.*P.M.*inf_matrix)
        max_eigval = Real(eigs[end])
        P.χ = ones(n_a)/max_eigval #This rescales everything so β is the same as R₀

    #P.β = β#rand(d_R₀) #Choose R₀ randomly from 2-3 range
    u0[4,8,3] = 5 #10 initial Asymptomatics in Nairobi;    #u0[12,8,3] = 10 #10 initial pre-symptomatics in Mombasa
    return u0,P
end




########
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
                            upred = quantile(sum(hosp_by_area_over_sims,dims=1)[:],0.975))

    severe_cases_by_area = (med = [median(hosp_by_area_over_sims[cn,:]) for cn = 1:nc],
                            lpred = [quantile(hosp_by_area_over_sims[cn,:],0.025) for cn = 1:nc],
                            upred = [quantile(hosp_by_area_over_sims[cn,:],0.975) for cn = 1:nc])
    severe_cases_by_age = (med = [median(hosp_by_age_over_sims[a,:]) for a = 1:na],
                            lpred = [quantile(hosp_by_age_over_sims[a,:],0.025) for a = 1:na],
                            upred = [quantile(hosp_by_age_over_sims[a,:],0.975) for a = 1:na])
    total_deaths = (med = median(sum(deaths_by_area_over_sims,dims=1)),
                            lpred = quantile(sum(deaths_by_area_over_sims,dims=1)[:],0.025),
                            upred = quantile(sum(deaths_by_area_over_sims,dims=1)[:],0.975))
    deaths_by_area = (med = [median(deaths_by_area_over_sims[cn,:]) for cn = 1:nc],
                            lpred = [quantile(deaths_by_area_over_sims[cn,:],0.025) for cn = 1:nc],
                            upred = [quantile(deaths_by_area_over_sims[cn,:],0.975) for cn = 1:nc])
    deaths_by_age = (med = [median(deaths_by_age_over_sims[a,:]) for a = 1:na],
                            lpred = [quantile(deaths_by_age_over_sims[a,:],0.025) for a = 1:na],
                            upred = [quantile(deaths_by_age_over_sims[a,:],0.975) for a = 1:na])
    hosp_peak_excess_demand_by_area = (med = [median(peak_hosp_by_area_by_sim[cn,:]/KenyaCoV_screening.spare_capacity_H_by_county[cn]) for cn = 1:nc],
                            lpred = [quantile(peak_hosp_by_area_by_sim[cn,:]/KenyaCoV_screening.spare_capacity_H_by_county[cn],0.025) for cn = 1:nc],
                            upred = [quantile(peak_hosp_by_area_by_sim[cn,:]/KenyaCoV_screening.spare_capacity_H_by_county[cn],0.975) for cn = 1:nc])

    ICU_peak_excess_demand_by_area = (med = [median(peak_ICU_by_area_by_sim[cn,:] .- KenyaCoV_screening.spare_capacity_ICU_by_county[cn]) for cn = 1:nc],
                            lpred = [quantile(peak_ICU_by_area_by_sim[cn,:] .- KenyaCoV_screening.spare_capacity_ICU_by_county[cn],0.025) for cn = 1:nc],
                            upred = [quantile(peak_ICU_by_area_by_sim[cn,:] .- KenyaCoV_screening.spare_capacity_ICU_by_county[cn],0.975) for cn = 1:nc])



    ##Added for final cumulative cases
    total_cases_A_by_area_and_age_med = zeros(nc,na)
    total_cases_A_by_area_and_age_lpred = zeros(nc,na)
    total_cases_A_by_area_and_age_upred = zeros(nc,na)
    for cn = 1:nc,an = 1:na
        total_cases_A_by_area_and_age_med[cn,an] = median([sims[k].total_cases_A_by_area_and_age[cn,an] for k=1:n])
        total_cases_A_by_area_and_age_lpred[cn,an] = quantile([sims[k].total_cases_A_by_area_and_age[cn,an] for k=1:n],0.025)
        total_cases_A_by_area_and_age_upred[cn,an] = quantile([sims[k].total_cases_A_by_area_and_age[cn,an] for k=1:n],0.975)
    end
    total_cases_A_by_area_and_age=(med=total_cases_A_by_area_and_age_med,
                                    lpred=total_cases_A_by_area_and_age_lpred,
                                    upred=total_cases_A_by_area_and_age_upred)

    total_cases_M_by_area_and_age_med = zeros(nc,na)
    total_cases_M_by_area_and_age_lpred = zeros(nc,na)
    total_cases_M_by_area_and_age_upred = zeros(nc,na)
    for cn = 1:nc,an = 1:na
        total_cases_M_by_area_and_age_med[cn,an] = median([sims[k].total_cases_M_by_area_and_age[cn,an] for k=1:n])
        total_cases_M_by_area_and_age_lpred[cn,an] = quantile([sims[k].total_cases_M_by_area_and_age[cn,an] for k=1:n],0.025)
        total_cases_M_by_area_and_age_upred[cn,an] = quantile([sims[k].total_cases_M_by_area_and_age[cn,an] for k=1:n],0.975)
    end
    total_cases_M_by_area_and_age=(med=total_cases_M_by_area_and_age_med,
                                    lpred=total_cases_M_by_area_and_age_lpred,
                                    upred=total_cases_M_by_area_and_age_upred)
    total_cases_V_by_area_and_age_med = zeros(nc,na)
    total_cases_V_by_area_and_age_lpred = zeros(nc,na)
    total_cases_V_by_area_and_age_upred = zeros(nc,na)
    for cn = 1:nc,an = 1:na
        total_cases_V_by_area_and_age_med[cn,an] = median([sims[k].total_cases_V_by_area_and_age[cn,an] for k=1:n])
        total_cases_V_by_area_and_age_lpred[cn,an] = quantile([sims[k].total_cases_V_by_area_and_age[cn,an] for k=1:n],0.025)
        total_cases_V_by_area_and_age_upred[cn,an] = quantile([sims[k].total_cases_V_by_area_and_age[cn,an] for k=1:n],0.975)
    end
    total_cases_V_by_area_and_age=(med=total_cases_V_by_area_and_age_med,
                                    lpred=total_cases_V_by_area_and_age_lpred,
                                    upred=total_cases_V_by_area_and_age_upred)

    total_q_by_area_and_age_med = zeros(nc,na)
    total_q_by_area_and_age_lpred = zeros(nc,na)
    total_q_by_area_and_age_upred = zeros(nc,na)
    for cn = 1:nc,an = 1:na
        total_q_by_area_and_age_med[cn,an] = median([sims[k].total_q_by_area_and_age[cn,an] for k=1:n])
        total_q_by_area_and_age_lpred[cn,an] = quantile([sims[k].total_q_by_area_and_age[cn,an] for k=1:n],0.025)
        total_q_by_area_and_age_upred[cn,an] = quantile([sims[k].total_q_by_area_and_age[cn,an] for k=1:n],0.975)
    end
    total_q_by_area_and_age=(med=total_q_by_area_and_age_med,
                                    lpred=total_q_by_area_and_age_lpred,
                                    upred=total_q_by_area_and_age_upred)

    total_cases_A=(med=median([sum(sims[k].total_cases_A_by_area_and_age[:,:]) for k=1:n]),
                    lpred=quantile([sum(sims[k].total_cases_A_by_area_and_age[:,:]) for k=1:n],0.025),
                    upred=quantile([sum(sims[k].total_cases_A_by_area_and_age[:,:]) for k=1:n],0.975))
    total_cases_M=(med=median([sum(sims[k].total_cases_M_by_area_and_age[:,:]) for k=1:n]),
                    lpred=quantile([sum(sims[k].total_cases_M_by_area_and_age[:,:]) for k=1:n],0.025),
                    upred=quantile([sum(sims[k].total_cases_M_by_area_and_age[:,:]) for k=1:n],0.975))
    total_cases=(med=median([sum(sims[k].total_cases_A_by_area_and_age[:,:])+sum(sims[k].total_cases_M_by_area_and_age[:,:])+sum(sims[k].total_cases_V_by_area_and_age[:,:]) for k=1:n]),
                    lpred=quantile([sum(sims[k].total_cases_A_by_area_and_age[:,:])+sum(sims[k].total_cases_M_by_area_and_age[:,:])+sum(sims[k].total_cases_V_by_area_and_age[:,:]) for k=1:n],0.025),
                    upred=quantile([sum(sims[k].total_cases_A_by_area_and_age[:,:])+sum(sims[k].total_cases_M_by_area_and_age[:,:])+sum(sims[k].total_cases_V_by_area_and_age[:,:]) for k=1:n],0.975))
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
            prevalence_ICU_ts=prevalence_ICU_ts,
            total_cases_A_by_area_and_age=total_cases_A_by_area_and_age,    ##
            total_cases_M_by_area_and_age=total_cases_M_by_area_and_age,    ##
            total_cases_V_by_area_and_age=total_cases_V_by_area_and_age,    ##
            total_cases_A=total_cases_A,                                    ##
            total_cases_M=total_cases_M,                                    ##
            total_cases=total_cases,                                        ##
            total_q_by_area_and_age=total_q_by_area_and_age)                ##
end

function generate_report_screening(folder,output,#=model_str,=#simulation_tag,areanames)
    report_folder=folder*"report_"*simulation_tag*"/"
    if !isdir(report_folder)   mkdir(report_folder)   end

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
            prevalence_ICU_ts=output.prevalence_ICU_ts#=,
            model_str=model_str=#)
    @save(report_folder*"scenario_data"*simulation_tag*".jld2",scenariodata)

    plt_incidence_mombasa,plt_health_usage_mombasa = KenyaCoV_screening.give_plots_for_county(output,[28],simulation_tag,areanames)
    savefig(plt_incidence_mombasa,report_folder*"incidence_mombasa"*simulation_tag*".png")
    savefig(plt_health_usage_mombasa,report_folder*"healthsystem_mombasa"*simulation_tag*".png")
    plt_incidence_nairobi,plt_health_usage_nairobi = KenyaCoV_screening.give_plots_for_county(output,[30],simulation_tag,areanames)
    savefig(plt_incidence_nairobi,report_folder*"incidence_nairobi"*simulation_tag*".png")
    savefig(plt_health_usage_nairobi,report_folder*"healthsystem_nairobi"*simulation_tag*".png")
    plt_incidence_restofcountry,plt_health_usage_restofcountry = KenyaCoV_screening.give_plots_for_county(output,setdiff(1:47,[28,30]),simulation_tag,areanames)
    savefig(plt_incidence_restofcountry,report_folder*"plt_incidence_restofcountry"*simulation_tag*".png")
    savefig(plt_health_usage_restofcountry,report_folder*"healthsystem_restofcountry"*simulation_tag*".png")

    plt_ranked_HU,plt_ranked_ICU = KenyaCoV_screening.plot_ranked_bars_health_usage(output,simulation_tag,areanames)
    savefig(plt_ranked_HU,report_folder*"peak_hospital_usage_by_county"*simulation_tag*".png")
    savefig(plt_ranked_ICU,report_folder*"peak_ICU_usage_by_county"*simulation_tag*".png")

    plt_severe,plt_deaths = KenyaCoV_screening.plot_ranked_bars_cases(output,simulation_tag,areanames)
    savefig(plt_severe,report_folder*"total_severe_cases_by_county"*simulation_tag*".png")
    savefig(plt_deaths,report_folder*"total_deaths_by_county"*simulation_tag*".png")

    #Total country output
    writedlm(report_folder*"country_A_incidence_ts"*simulation_tag*".csv",output.country_incidence_A_ts.med,",")
    writedlm(report_folder*"country_cum_A_incidence_ts"*simulation_tag*".csv",cumsum(output.country_incidence_A_ts.med),",")
    writedlm(report_folder*"country_M_incidence_ts"*simulation_tag*".csv",output.country_incidence_M_ts.med,",")
    writedlm(report_folder*"country_cum_M_incidence_ts"*simulation_tag*".csv",cumsum(output.country_incidence_M_ts.med),",")
    writedlm(report_folder*"country_V_incidence_ts"*simulation_tag*".csv",output.country_incidence_V_ts.med,",")
    writedlm(report_folder*"country_cum_V_incidence_ts"*simulation_tag*".csv",cumsum(output.country_incidence_V_ts.med),",")
    writedlm(report_folder*"country_all_incidence_ts"*simulation_tag*".csv",output.country_incidence_A_ts.med .+ output.country_incidence_M_ts.med .+ output.country_incidence_V_ts.med,",")
    writedlm(report_folder*"country_cum_all_incidence_ts"*simulation_tag*".csv",cumsum(output.country_incidence_A_ts.med .+ output.country_incidence_M_ts.med .+ output.country_incidence_V_ts.med),",")
    writedlm(report_folder*"country_deaths_ts"*simulation_tag*".csv",output.country_incidence_death_ts.med,",")
    writedlm(report_folder*"country_cum_deaths_ts"*simulation_tag*".csv",cumsum(output.country_incidence_death_ts.med),",")

    writedlm(report_folder*"country_hosp_prev"*simulation_tag*".csv",output.country_prevalence_H_ts,",")
    writedlm(report_folder*"country_ICU_prev"*simulation_tag*".csv",output.country_prevalence_ICU_ts,",")


    writedlm(report_folder*"A_incidence_ts_by_county"*simulation_tag*".csv",hcat(areanames,output.incidence_A_ts.med),",")
    writedlm(report_folder*"cum_A_incidence_ts_by_county"*simulation_tag*".csv",hcat(areanames,cumsum(output.incidence_A_ts.med,dims=2)),",")
    writedlm(report_folder*"M_incidence_ts_by_county"*simulation_tag*".csv",hcat(areanames,output.incidence_M_ts.med),",")
    writedlm(report_folder*"cum_M_incidence_ts_by_county"*simulation_tag*".csv",hcat(areanames,cumsum(output.incidence_M_ts.med,dims=2)),",")
    writedlm(report_folder*"V_incidence_ts_by_county"*simulation_tag*".csv",hcat(areanames,output.incidence_V_ts.med),",")
    writedlm(report_folder*"cum_V_incidence_ts_by_county"*simulation_tag*".csv",hcat(areanames,cumsum(output.incidence_V_ts.med,dims=2)),",")
    writedlm(report_folder*"all_incidence_ts_by_county"*simulation_tag*".csv",hcat(areanames,output.incidence_A_ts.med .+ output.incidence_M_ts.med .+ output.incidence_V_ts.med),",")
    writedlm(report_folder*"cum_all_incidence_ts_by_county"*simulation_tag*".csv",hcat(areanames,cumsum(output.incidence_A_ts.med .+ output.incidence_M_ts.med .+ output.incidence_V_ts.med,dims=2)),",")

    writedlm(report_folder*"death_incidence_ts_by_county"*simulation_tag*".csv",hcat(areanames,output.incidence_death_ts.med),",")
    writedlm(report_folder*"cum_death_incidence_ts_by_county"*simulation_tag*".csv",hcat(areanames,cumsum(output.incidence_death_ts.med,dims=2)),",")

    writedlm(report_folder*"hosp_prev_ts_by_county"*simulation_tag*".csv",hcat(areanames,output.prevalence_H_ts.med),",")
    writedlm(report_folder*"ICU_prev_ts_by_county"*simulation_tag*".csv",hcat(areanames,output.prevalence_ICU_ts.med),",")

    return scenariodata
end

function give_plots_for_county(output,areas_to_plot,simulation_tag,areanames)
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
                        title = "Incidence of disease - "*area_name[1]*simulation_tag)
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
                            title = "Health system usage - "*area_name[1]*simulation_tag)

    plot!(plt_health_usage,0:(T-1),median_ICU_occup,
                            lab = "ICU beds",
                            lw = 3,color = :green,
                            ribbon = (median_ICU_occup.-lb_ICU_occup,ub_ICU_occup .- median_ICU_occup),
                            fillalpha = 0.25)
    if length(areas_to_plot) == 1
        plot!(plt_health_usage,[0,T-1],[spare_capacity_H_by_county[areas_to_plot][1],spare_capacity_H_by_county[areas_to_plot][1]],
                lw = 2,ls = :dash,lab = "spare hosp. capacity",color = :blue)

        plot!(plt_health_usage,[0,T-1],[spare_capacity_ICU_by_county[areas_to_plot][1],spare_capacity_ICU_by_county[areas_to_plot][1]],
                lw = 2,ls = :dash,lab = "spare ICU capacity",color = :green)
    end

    return plt_incidence,plt_health_usage
end

function plot_ranked_bars_health_usage(output,simulation_tag,areanames)
    median_hospital_exceed = output.hosp_peak_excess_demand_by_area.med
    median_ICU_exceed = output.ICU_peak_excess_demand_by_area.med
    nc = length(areanames)

    I = sortperm(median_hospital_exceed)
    plt_HU = bar(median_hospital_exceed[I]*100,orientations = :horizonal,
                yticks = (1:47,areanames[I]),size = (700,550),lab ="",
                xlabel = "% Hospital available bed demand",title = "Hospital demand at peak"*simulation_tag)
    scatter!(plt_HU,median_hospital_exceed[I]*100,1:nc,
                xerror = ((median_hospital_exceed[I] .- output.hosp_peak_excess_demand_by_area.lpred[I])*100
                            ,(output.hosp_peak_excess_demand_by_area.upred[I] .- median_hospital_exceed[I])*100 ),
                ms = 0.,color = :black,lab ="")

    I = sortperm(median_ICU_exceed)
    plt_ICU = bar(median_ICU_exceed[I,1],orientations = :horizonal,
                yticks = (1:nc,areanames[I]),size = (700,550),lab ="",
                xlabel = "ICU bed excess demand (numbers)",title = "ICU demand at peak"*simulation_tag)
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

function run_βrand_withHreports(P::KenyaCoV_screening.CoVParameters_Screening,prob,n_traj,s,sc,folder,simulation_tag,areanames;interventions = CallbackSet())
    #sims = run_simulations(P,prob,n_traj;interventions=interventions)
    sims = run_simulations_βrand_withH(P,prob,n_traj;interventions=interventions)
    output = extract_information_from_simulations(sims);
    scenariodata = generate_report_screening(folder,output,simulation_tag,areanames,n_traj,s,sc)
    return scenariodata
end

function run_βrand(P::KenyaCoV_screening.CoVParameters_Screening,prob,n_traj;interventions = CallbackSet())
    sims = run_simulations_βrand_withH(P,prob,n_traj;interventions=interventions)
    output = extract_information_from_simulations(sims);
    return output
end

counties = CSV.read("data/2019_census_age_pyramids_counties.csv")
