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
d_R₀ = Gamma(100,2/100) ##Liu et al
mean(d_R₀)
(quantile(d_R₀,0.025),median(d_R₀),quantile(d_R₀,0.975))


######## Function to initialize dynamics parameters
interventions=[#=1=#"Contact tracing of hospitalized only",#=2=#"Symptomatic screening",#=3=#"Symptomatic screening with contact tracing",
                #=4=#"Mass screening",#=5=#"Mass screening with contact tracing"]

##############
function baseline_config(ϵ)
    @load "data/detection_rates_for_different_taus.jld2" d_0 d_01 d_025 d_05 d_1
    #u0,P,P_dest = model_ingredients_from_data_screening("data/data_for_age_structuredmodel.jld2","data/flight_numbers.csv","data/projected_global_prevelance.csv")
    u0,P,P_dest = model_ingredients_from_data_screening("data/data_for_age_structuredmodel_with_counties.jld2","data/flight_numbers.csv","data/projected_global_prevelance.csv")
        @load "data/detection_rates_for_different_taus.jld2" d_0 d_01 d_025 d_05 d_1
        #Redistribute susceptibility JUST to rescale infectiousness so we get the correct R₀/r
        P.χ = ones(n_a)
        P.rel_detection_rate = d_1
        P.dt = 0.25
        P.ext_inf_rate = 0.
        P.ϵ = ϵ
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
        for r=1:n,a=1:n_a,c_a=1:n_a
            P.M_rescaled[r,c_a,a]=P.β*P.M[c_a,a]/P.N̂[r,a]
        end

    #P.β = rand(d_R₀) #Choose R₀ randomly from 1.5-2.5 range
    u0[4,8,3] = 5 #10 initial Asymptomatics in Nairobi;    #u0[12,8,3] = 10 #10 initial pre-symptomatics in Mombasa
    return u0,P
end
#######
function run_simulations(P::KenyaCoV_screening.CoVParameters_Screening,prob,n_traj,R0_is_random,is_distributed;interventions = CallbackSet())
    ensemble_prob=0;sims=0;
    if R0_is_random
        ensemble_prob = EnsembleProblem(prob,
                                        prob_func = randomise_R₀,
                                        output_func = output_simulation_data)
    else
        ensemble_prob = EnsembleProblem(prob,
                                        prob_func = remake_prob,
                                        output_func = output_simulation_data)
    end
    if is_distributed
        sims=solve(ensemble_prob,EnsembleDistributed(),dt = P.dt,callback = interventions,trajectories = n_traj)
    else
        sims=solve(ensemble_prob,FunctionMap(),dt = P.dt,callback = interventions,trajectories = n_traj)
    end
    output = extract_information_from_simulations(sims)
    return output
end

#=function run_βrand(P::KenyaCoV_screening.CoVParameters_Screening,prob,n_traj;interventions = CallbackSet())
    sims = run_simulations_βrand(P,prob,n_traj;interventions=interventions)
    output = extract_information_from_simulations(sims);
    return output
end

function run_βfixed(P::KenyaCoV_screening.CoVParameters_Screening,prob,n_traj;interventions = CallbackSet())
    sims = run_simulations_βfixed(P,prob,n_traj;interventions=interventions)
    output = extract_information_from_simulations(sims);
    return output
end

function run_simulations_βrand(P::CoVParameters_Screening,prob,n_traj;interventions = CallbackSet())
    ensemble_prob = EnsembleProblem(prob,
                                    prob_func = randomise_R₀,
                                    output_func = output_simulation_data)
    return solve(ensemble_prob,FunctionMap(),dt = P.dt,callback = interventions,trajectories = n_traj)
end

function run_simulations_βfixed(P::CoVParameters_Screening,prob,n_traj;interventions = CallbackSet())
    ensemble_prob = EnsembleProblem(prob,
                                    prob_func = remake_prob,
                                    output_func = output_simulation_data)
    return solve(ensemble_prob,FunctionMap(),dt = P.dt,callback = interventions,trajectories = n_traj)
end

function run_βfixed_distributed(P::KenyaCoV_screening.CoVParameters_Screening,prob,n_traj;interventions = CallbackSet())
    sims = run_simulations_βfixed_distributed(P,prob,n_traj;interventions=interventions)
    output = extract_information_from_simulations(sims);
    return output
end

function run_simulations_βfixed_distributed(P::CoVParameters_Screening,prob,n_traj;interventions = CallbackSet())
    ensemble_prob = EnsembleProblem(prob,
                                    prob_func = remake_prob,
                                    output_func = output_simulation_data)
    return solve(ensemble_prob,EnsembleDistributed(),dt = P.dt,callback = interventions,trajectories = n_traj)
end

function run_βrand_distributed(P::KenyaCoV_screening.CoVParameters_Screening,prob,n_traj;interventions = CallbackSet())
    sims = run_simulations_βrand_distributed(P,prob,n_traj;interventions=interventions)
    output = extract_information_from_simulations(sims);
    return output
end

function run_simulations_βrand_distributed(P::CoVParameters_Screening,prob,n_traj;interventions = CallbackSet())
    ensemble_prob = EnsembleProblem(prob,
                                    prob_func = randomise_R₀,
                                    output_func = output_simulation_data)
    return solve(ensemble_prob,EnsembleDistributed(),dt = P.dt,callback = interventions,trajectories = n_traj)
end=#

counties = CSV.read("data/2019_census_age_pyramids_counties.csv")

########
function randomise_R₀(prob,i,repeat) #Remember to rescale susceptibility by the inverse leading eigenvalue
    _P = deepcopy(prob.p)
    _P.β = rand(d_R₀)
    return remake(prob,p=_P)
end
function remake_prob(prob,i,repeat) #Remember to rescale susceptibility by the inverse leading eigenvalue
    _P = deepcopy(prob.p)
    #_P.β = rand(d_R₀)
    return remake(prob,p=_P)
end
#=function output_daily_incidence_and_hosp_cums(sol,i)
    times = sol.prob.tspan[1]:1:sol.prob.tspan[end]
    cumA=sol[end][:,:,9]
    cumM=sol[end][:,:,10]
    cumV=sol[end][:,:,11]
    cumH=sol[end][:,:,12]
    cumQ=sol[end][:,:,16]
    I = [sum(sol(t)[:,:,4:6])  for t in times]
    return (cumA,cumM,cumV,cumH,cumQ,I),false
end=#
function output_simulation_data(sol,i)
    times = sol.prob.tspan[1]:1:sol.prob.tspan[end]

    total_cases_A_by_area_and_age=sol(sol.prob.tspan[end])[:,:,9]   ##
    total_cases_M_by_area_and_age=sol(sol.prob.tspan[end])[:,:,10]  ##
    total_cases_V_by_area_and_age=sol(sol.prob.tspan[end])[:,:,11]  ##
    total_q_by_area_and_age=sol(sol.prob.tspan[end])[:,:,16]        ##
    total_qs_by_area_and_age=sol(sol.prob.tspan[end])[:,:,17]       ##

    incidence_A = diff([sum(sol(t)[:,:,9],dims=2)[:]  for t in times])
    incidence_M = diff([sum(sol(t)[:,:,10],dims=2)[:]  for t in times])
    incidence_V = diff([sum(sol(t)[:,:,11],dims=2)[:]  for t in times])
    incidence_I = incidence_A + incidence_M + incidence_V           ##
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
            incidence_I=VectorOfArray(incidence_I)[:,:],                        ##
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
            total_cases_V_by_area_and_age=total_cases_V_by_area_and_age,        ##
            total_q_by_area_and_age=total_q_by_area_and_age,
            total_qs_by_area_and_age=total_qs_by_area_and_age),false            ##
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

    country_incidence_I = VectorOfArray([sum(sims[k].incidence_I,dims=1) for k = 1:n])[1,:,:]       ##
    incidence_I_by_area_and_sims = VectorOfArray([sims[k].incidence_I for k = 1:n])[:,:,:]          ##
    incidence_I_by_area_med = zeros(nc,T)                                                           ##
    incidence_I_by_area_lpred = zeros(nc,T)                                                         ##
    incidence_I_by_area_upred = zeros(nc,T)                                                         ##
    for cn = 1:nc,t = 1:T                                                                           ##
        incidence_I_by_area_med[cn,t] = median(incidence_I_by_area_and_sims[cn,t,:])                ##
        incidence_I_by_area_lpred[cn,t] = quantile(incidence_I_by_area_and_sims[cn,t,:],0.025)      ##
        incidence_I_by_area_upred[cn,t] = quantile(incidence_I_by_area_and_sims[cn,t,:],0.975)      ##
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

    country_incidence_I_ts = (med = [median(country_incidence_I[t,:]) for t = 1:T],             ##
                            lpred = [quantile(country_incidence_I[t,:],0.025) for t = 1:T],     ##
                            upred = [quantile(country_incidence_I[t,:],0.025) for t = 1:T])     ##
    incidence_I_ts = (med = incidence_I_by_area_med,                                            ##
                    lpred = incidence_I_by_area_lpred,                                          ##
                    upred = incidence_I_by_area_upred)                                          ##

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

    total_qs_by_area_and_age_med = zeros(nc,na)
    total_qs_by_area_and_age_lpred = zeros(nc,na)
    total_qs_by_area_and_age_upred = zeros(nc,na)
    for cn = 1:nc,an = 1:na
        total_qs_by_area_and_age_med[cn,an] = median([sims[k].total_qs_by_area_and_age[cn,an] for k=1:n])
        total_qs_by_area_and_age_lpred[cn,an] = quantile([sims[k].total_qs_by_area_and_age[cn,an] for k=1:n],0.025)
        total_qs_by_area_and_age_upred[cn,an] = quantile([sims[k].total_qs_by_area_and_age[cn,an] for k=1:n],0.975)
    end
    total_qs_by_area_and_age=(med=total_qs_by_area_and_age_med,
                                    lpred=total_qs_by_area_and_age_lpred,
                                    upred=total_qs_by_area_and_age_upred)

    total_cases_A=(med=median([sum(sims[k].total_cases_A_by_area_and_age[:,:]) for k=1:n]),
                    lpred=quantile([sum(sims[k].total_cases_A_by_area_and_age[:,:]) for k=1:n],0.025),
                    upred=quantile([sum(sims[k].total_cases_A_by_area_and_age[:,:]) for k=1:n],0.975))
    total_cases_M=(med=median([sum(sims[k].total_cases_M_by_area_and_age[:,:]) for k=1:n]),
                    lpred=quantile([sum(sims[k].total_cases_M_by_area_and_age[:,:]) for k=1:n],0.025),
                    upred=quantile([sum(sims[k].total_cases_M_by_area_and_age[:,:]) for k=1:n],0.975))
    total_cases_MV=(med=median([sum(sims[k].total_cases_M_by_area_and_age[:,:])+sum(sims[k].total_hosp_by_area_and_age) for k=1:n]),     ###
                    lpred=quantile([sum(sims[k].total_cases_M_by_area_and_age[:,:])+sum(sims[k].total_hosp_by_area_and_age) for k=1:n],0.025),
                    upred=quantile([sum(sims[k].total_cases_M_by_area_and_age[:,:])+sum(sims[k].total_hosp_by_area_and_age) for k=1:n],0.975))
    total_cases=(med=median([sum(sims[k].total_cases_A_by_area_and_age[:,:])+sum(sims[k].total_cases_M_by_area_and_age[:,:])+sum(sims[k].total_cases_V_by_area_and_age[:,:]) for k=1:n]),
                    lpred=quantile([sum(sims[k].total_cases_A_by_area_and_age[:,:])+sum(sims[k].total_cases_M_by_area_and_age[:,:])+sum(sims[k].total_cases_V_by_area_and_age[:,:]) for k=1:n #=if sum(sims[k].total_cases_A_by_area_and_age[:,:])+sum(sims[k].total_cases_M_by_area_and_age[:,:])+sum(sims[k].total_cases_V_by_area_and_age[:,:])>1e5=#],0.025),
                    upred=quantile([sum(sims[k].total_cases_A_by_area_and_age[:,:])+sum(sims[k].total_cases_M_by_area_and_age[:,:])+sum(sims[k].total_cases_V_by_area_and_age[:,:]) for k=1:n],0.975))
    total_q=(med=median([sum(sims[k].total_q_by_area_and_age[:,:]) for k=1:n]),
                    lpred=quantile([sum(sims[k].total_q_by_area_and_age[:,:]) for k=1:n],0.025),
                    upred=quantile([sum(sims[k].total_q_by_area_and_age[:,:]) for k=1:n],0.975))
    total_qs=(med=median([sum(sims[k].total_qs_by_area_and_age[:,:]) for k=1:n]),
                    lpred=quantile([sum(sims[k].total_qs_by_area_and_age[:,:]) for k=1:n],0.025),
                    upred=quantile([sum(sims[k].total_qs_by_area_and_age[:,:]) for k=1:n],0.975))
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
            country_incidence_I_ts=country_incidence_I_ts,                  ##
            incidence_I_ts=incidence_I_ts,                                  ##
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
            total_cases_MV=total_cases_MV,                                   ###
            total_cases=total_cases,                                        ##
            total_q_by_area_and_age=total_q_by_area_and_age,                ##
            total_qs_by_area_and_age=total_qs_by_area_and_age,              ##
            total_q=total_q,                ##
            total_qs=total_qs)              ##
end
