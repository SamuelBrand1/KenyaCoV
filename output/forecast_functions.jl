"""
Functions for analysing simulations

Output functions: Target the peak timing for each county, cases by each county,
timing of total peak and total cases
"""

function output_daily_and_final_incidence(sol,i)
    times = sol.prob.tspan[1]:1:sol.prob.tspan[end]
    z = [sum(sol(t)[:,:,7:8],dims = 2)[:,1,:]  for t in times]
    return (z,sol[end][:,:,7:8]),false
end

function randomise_params(prob,i,repeat)
    _P = deepcopy(prob.p)
    _P.σ = 1/rand(d_incubation)
    _P.β = rand(d_R₀)*_P.γ/(_P.δ + _P.ϵ*(1-_P.δ))
    return remake(prob,p=_P)
end

"""
Simulation functions
"""
function run_simulations(P::KenyaCoV.CoVParameters_AS,prob,n_traj,τ)
    P.τ = τ
    ensemble_prob = EnsembleProblem(prob,
                                    prob_func = randomise_params,
                                    output_func = output_daily_and_final_incidence)
    return solve(ensemble_prob,FunctionMap(),dt = P.dt,trajectories = n_traj)
end

function run_scenario(P::KenyaCoV.CoVParameters_AS,prob,n_traj,treatment_rates)
    results = []
    for τ in treatment_rates
        sims = run_simulations(P,prob,n_traj,τ)
        analysisdata = incidence_from_sims(sims)
        push!(results,analysisdata)
    end
    return results
end

"""
Analysis functions:
"""
function incidence_from_sims(sims)
    n = length(sims)
    m = length(sims[1][1])
    inc_A_data = zeros(21,m-1,n)
    inc_D_data = zeros(21,m-1,n)
    final_case_data = zeros(21,16,n)
    for (k,sim) in enumerate(sims)
        for i = 1:20,t = 1:(m-1)
            inc_A_data[i,t,k] = sim[1][t+1][i,1] - sim[1][t][i,1]#Daily incidence rather than cumulative
            inc_D_data[i,t,k] = sim[1][t+1][i,2] - sim[1][t][i,2]
        end
        for t = 1:(m-1)
            inc_A_data[21,t,k] = sum(sim[1][t+1][:,1] .- sim[1][t][:,1])
            inc_D_data[21,t,k] = sum(sim[1][t+1][:,2] .- sim[1][t][:,2])
        end
        for i =1:20,a=1:16
            final_case_data[i,a,k] = sum(sim[2][i,a,:])
        end
        for a = 1:16
            final_case_data[21,a,k] = sum(sim[2][:,a,:])
        end
    end
    peak_times = zeros(n,21)
    inc_A_conf_intvs = zeros(21,m-1,3)
    inc_D_conf_intvs = zeros(21,m-1,3)
    final_case_intvs = zeros(21,16,3)
    for i = 1:21,t = 1:(m-1)
        inc_A_conf_intvs[i,t,1] = quantile(inc_A_data[i,t,:],0.5)
        inc_A_conf_intvs[i,t,2] = inc_A_conf_intvs[i,t,1] - quantile(inc_A_data[i,t,:],0.025) #Lower conf. int.
        inc_A_conf_intvs[i,t,3] =  quantile(inc_A_data[i,t,:],0.975) - inc_A_conf_intvs[i,t,1] #Upper conf. int.
        inc_D_conf_intvs[i,t,1] = quantile(inc_D_data[i,t,:],0.5)
        inc_D_conf_intvs[i,t,2] = inc_D_conf_intvs[i,t,1] - quantile(inc_D_data[i,t,:],0.025) #Lower conf. int.
        inc_D_conf_intvs[i,t,3] =  quantile(inc_D_data[i,t,:],0.975) - inc_D_conf_intvs[i,t,1] #Upper conf. int.
    end

    for i = 1:21, k =1:n
        peak_times[k,i] = argmax(inc_A_data[i,:,k])
    end

    for i = 1:21,a=1:16
        final_case_intvs[i,a,1] = quantile(final_case_data[i,a,:],0.5)
        final_case_intvs[i,a,2] = quantile(final_case_data[i,a,:],0.5) - quantile(final_case_data[i,a,:],0.025)
        final_case_intvs[i,a,3] = quantile(final_case_data[i,a,:],0.975) - quantile(final_case_data[i,a,:],0.5)
    end



    return inc_A_conf_intvs,inc_D_conf_intvs,peak_times,final_case_intvs
end
