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

"""
Simulation functions
"""


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
            total_deaths_by_area_and_age = deaths_by_area_and_age),false
end




function randomise_params(prob,i,repeat) #Remember to rescale susceptibility by the inverse leading eigenvalue
    _P = deepcopy(prob.p)
    _P.isolating_detecteds = true
    _P.τ = _P.τ_initial
    _P.σ = 1/rand(d_incubation)
    _P.β = rand(d_R₀)*_P.γ
    return remake(prob,p=_P)
end

function consensus_randomise_params(prob,i,repeat) #Remember to rescale susceptibility by the inverse leading eigenvalue
    _P = deepcopy(prob.p)
    _P.β = rand(d_R₀)
    return remake(prob,p=_P)
end

function randomise_params_and_infectiousness(prob,i,repeat)
    _P = deepcopy(prob.p)
    _P.isolating_detecteds = true
    _P.τ = _P.τ_initial
    _P.σ = 1/rand(d_incubation)
    _P.β = rand(d_R₀)*_P.γ
    _P.ϵ = rand(Uniform(0.,0.5))
    return remake(prob,p=_P)
end

"""
Scenario functions
"""
# function run_simulations(P::KenyaCoV.CoVParameters_AS,prob,n_traj,τ,ϵ_D)
#     P.τ_initial = τ
#     P.ϵ_D = ϵ_D
#     ensemble_prob = EnsembleProblem(prob,
#                                     prob_func = randomise_params,
#                                     output_func = output_daily_and_final_incidence)
#     return solve(ensemble_prob,FunctionMap(),dt = P.dt,trajectories = n_traj)
# end

function run_simulations(P::KenyaCoV.CoVParameters_AS,prob,n_traj,τ,ϵ_D,cb)
    P.τ_initial = τ
    P.ϵ_D = ϵ_D
    ensemble_prob = EnsembleProblem(prob,
                                    prob_func = randomise_params,
                                    output_func = output_daily_and_final_incidence)
    return solve(ensemble_prob,FunctionMap(),dt = P.dt,callback = cb,trajectories = n_traj)
end

function run_consensus_simulations(P::KenyaCoV.CoVParameters_AS,prob,n_traj,cb)
    ensemble_prob = EnsembleProblem(prob,
                                    prob_func = consensus_randomise_params,
                                    output_func = output_simulation_data)
    return solve(ensemble_prob,FunctionMap(),dt = P.dt,callback = cb,trajectories = n_traj)
end

# function run_simulations(P::KenyaCoV.CoVParameters_AS,prob,n_traj,τ,ϵ_D,prob_func)
#     P.τ_initial = τ
#     P.ϵ_D = ϵ_D
#     ensemble_prob = EnsembleProblem(prob,
#                                     prob_func = prob_func,
#                                     output_func = output_daily_and_final_incidence)
#     return solve(ensemble_prob,FunctionMap(),dt = P.dt,trajectories = n_traj)
# end

# function run_scenario(P::KenyaCoV.CoVParameters_AS,prob,n_traj,treatment_rates)
#     results = []
#     for (τ,ϵ_D) in treatment_rates
#         sims = run_simulations(P,prob,n_traj,τ,ϵ_D)
#         analysisdata = incidence_from_sims(sims)
#         push!(results,analysisdata)
#     end
#     return results
# end

function run_scenario(P::KenyaCoV.CoVParameters_AS,prob,n_traj,treatment_rates,cb)
    results = []
    for (τ,ϵ_D) in treatment_rates
        sims = run_simulations(P,prob,n_traj,τ,ϵ_D,cb)
        analysisdata = incidence_from_sims(sims)
        push!(results,analysisdata)
    end
    return results
end



"""
Analysis functions
"""
function incidence_from_sims(sims)
    n = length(sims)
    m = length(sims[1][1])
    inc_A_data = zeros(n_wa+1,m-1,n)
    inc_D_data = zeros(n_wa+1,m-1,n)
    final_D = zeros(n_wa,n_a,n)
    final_A = zeros(n_wa,n_a,n)
    for (k,sim) in enumerate(sims)
        for i = 1:n_wa,t = 1:(m-1)
            inc_A_data[i,t,k] = sim[1][t+1][i,1] - sim[1][t][i,1]#Daily incidence rather than cumulative
            inc_D_data[i,t,k] = sim[1][t+1][i,2] - sim[1][t][i,2]
        end
        for t = 1:(m-1)
            inc_A_data[n_wa+1,t,k] = sum(sim[1][t+1][:,1] .- sim[1][t][:,1])
            inc_D_data[n_wa+1,t,k] = sum(sim[1][t+1][:,2] .- sim[1][t][:,2])
        end
        final_A[:,:,k] .= sim[2][:,:,1]
        final_D[:,:,k] .= sim[2][:,:,2]
    end
    peak_times = zeros(n,21)
    inc_A_conf_intvs = zeros(21,m-1,3)
    inc_D_conf_intvs = zeros(21,m-1,3)
    for i = 1:(n_wa+1),t = 1:(m-1)
        inc_A_conf_intvs[i,t,1] = quantile(inc_A_data[i,t,:],0.5)
        inc_A_conf_intvs[i,t,2] = inc_A_conf_intvs[i,t,1] - quantile(inc_A_data[i,t,:],0.025) #Lower conf. int.
        inc_A_conf_intvs[i,t,3] =  quantile(inc_A_data[i,t,:],0.975) - inc_A_conf_intvs[i,t,1] #Upper conf. int.
        inc_D_conf_intvs[i,t,1] = quantile(inc_D_data[i,t,:],0.5)
        inc_D_conf_intvs[i,t,2] = inc_D_conf_intvs[i,t,1] - quantile(inc_D_data[i,t,:],0.025) #Lower conf. int.
        inc_D_conf_intvs[i,t,3] =  quantile(inc_D_data[i,t,:],0.975) - inc_D_conf_intvs[i,t,1] #Upper conf. int.
    end

    for i = 1:(n_wa+1), k =1:n
        peak_times[k,i] = argmax(inc_A_data[i,:,k])
    end

    return inc_D_conf_intvs,inc_A_conf_intvs,peak_times,final_D,final_A
end
