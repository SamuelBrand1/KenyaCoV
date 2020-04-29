push!(LOAD_PATH, joinpath(homedir(),"GitHub/KenyaCoV/src"))
using Plots,Parameters,Distributions,DifferentialEquations,JLD2,DataFrames,StatsPlots,FileIO,MAT,RecursiveArrayTools
using Revise
import KenyaCoV
using LinearAlgebra:eigen
using Statistics: median, quantile

gr()#Plotting frontend
"""
Consensus modelling
"""

"""
Load age structured data, define initial state and declare the KenyaCoV problem
"""

u0,P,P_dest = KenyaCoV.model_ingredients_from_data("data/data_for_age_structuredmodel.jld2",
                                            "data/flight_numbers.csv",
                                            "data/projected_global_prevelance.csv")
#Put in the lockdown effects
T_normal = deepcopy(P.T)
T_regional_lockdown = deepcopy(P.T)
#Outgoing travel
#Nairobi
T_regional_lockdown[:,4] *= 0.1;T_regional_lockdown[4,4] += 1 - sum(T_regional_lockdown[:,4])
#Mombasa + Kwale
T_regional_lockdown[:,12] *= 0.1;T_regional_lockdown[12,12] += 1 - sum(T_regional_lockdown[:,12])
#Kilifi + Kwale
T_regional_lockdown[:,20] *= 0.1;T_regional_lockdown[20,20] += 1 - sum(T_regional_lockdown[:,20])
#Incoming travel
for leaving_area in 1:20,arriving_area in 1:20
    if !(leaving_area in [4,12,20]) && arriving_area in [4,12,20]
        T_regional_lockdown[arriving_area,leaving_area] *= 0.1
    end
end





@load "data/detection_rates_for_different_taus.jld2" d_0 d_01 d_025 d_05 d_1
@load "data/susceptibility_rates.jld2" σ

#Initial infecteds
"""
Load age mixing matrices (these are all in to (row) from (col) format)
"""

@load "data/agemixingmatrix_Kenya_all_types.jld2" M_Kenya M_Kenya_ho M_Kenya_other M_Kenya_school M_Kenya_work


#Function for chaning contacts from 100% to 45% over 14 days
function ramp_down(t)
    if t < 0.
        return 1.
    elseif t>= 0. || t <= 14.
        return (1-t/14.) + 0.45*t/14
    elseif t > 14.
        return 0.45
    end
end

#Put in the regional lockdown on day 25 (April 7th)
function regional_lockdown_timing(u,t,integrator)
  integrator.p.lockdown && t > 25.
end
function affect_regional_lockdown!(integrator)
  integrator.p.T = T_regional_lockdown
  integrator.p.lockdown = true
end
cb_regional_lockdown = DiscreteCallback(regional_lockdown_timing,affect_regional_lockdown!)


# both_cbs = CallbackSet(cb_iso_limit,cb_SD_limit)


"""
SCENARIO no intervention baseline
"""



u0,P,P_dest = KenyaCoV.model_ingredients_from_data("data/data_for_age_structuredmodel.jld2",
                                            "data/flight_numbers.csv",
                                            "data/projected_global_prevelance.csv")
#Redistribute susceptibility JUST to rescale infectiousness so we get the correct R₀/r
P.χ = ones(KenyaCoV.n_a)
P.rel_detection_rate = d_1
P.dt = 0.25;
P.ext_inf_rate = 0.;
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
P.χ = ones(KenyaCoV.n_a)/max_eigval

P.β = rand(KenyaCoV.d_R₀) #Choose R₀ randomly from 2-3 range

u0[4,8,3] = 30 #10 initial Asymptomatics in Nairobi
u0[12,8,3] = 10 #10 initial pre-symptomatics in Mombasa
prob = KenyaCoV.create_KenyaCoV_non_neg_prob(u0,(0.,1*365.),P)

sims_baseline = KenyaCoV.run_consensus_simulations(P::KenyaCoV.CoVParameters_AS,prob,1,CallbackSet())
P.χ .*= 1.384
sims_baseline_scaled = KenyaCoV.run_consensus_simulations(P::KenyaCoV.CoVParameters_AS,prob,1,CallbackSet())

@save joinpath(homedir(),"Github/KenyaCoVOutputs/sims_consensus_baseline.jld2") sims_baseline
@save joinpath(homedir(),"Github/KenyaCoVOutputs/sims_consensus_baseline_scaled.jld2") sims_baseline_scaled


println("Finished baseline sims for consensus modelling ")



"""
SCENARIO 2 --- Imperfect social distancing + spatial lockdown:
* No age-specific susceptibilties --- its disease difference
* Asymptomatics are 10% as infectious as symptomatics
"""



u0,P,P_dest = KenyaCoV.model_ingredients_from_data("data/data_for_age_structuredmodel.jld2",
                                            "data/flight_numbers.csv",
                                            "data/projected_global_prevelance.csv")
#Redistribute susceptibility JUST to rescale infectiousness so we get the correct R₀/r
P.χ = ones(KenyaCoV.n_a)
P.rel_detection_rate = d_1
P.dt = 0.25;
P.ext_inf_rate = 0.;
P.ϵ = 1.
P.c_t = ramp_d
#Set the susceptibility vector --- just to specify the correct R₀
sus_matrix = repeat(P.χ,1,17)
R_A = P.ϵ*((1/P.σ₂) + (1/P.γ) ) #effective duration of asymptomatic
R_M = (P.ϵ/P.σ₂) + (P.ϵ_D/P.γ) #effective duration of mild
R_V = (P.ϵ/P.σ₂) + (P.ϵ_V/P.τ) #effective duration of severe
R_vector = [(1-P.rel_detection_rate[a])*R_A + P.rel_detection_rate[a]*(1-P.hₐ[a])*R_M + P.rel_detection_rate[a]*P.hₐ[a]*R_V for a = 1:17]
inf_matrix = repeat(R_vector',17,1)

eigs, = eigen(sus_matrix.*P.M.*inf_matrix)
max_eigval = Real(eigs[end])
P.χ = ones(KenyaCoV.n_a)/max_eigval

P.β = 5.#rand(KenyaCoV.d_R₀) #Choose R₀ randomly from mean 2.5 (2-3) range
P.c_t = ramp_down

u0[4,8,3] = 30 #30 initial Asymptomatics in Nairobi
u0[12,8,3] = 10 #10 initial pre-symptomatics in Mombasa
prob = KenyaCoV.create_KenyaCoV_non_neg_prob(u0,(0.,1*365.),P)

sims_controls = KenyaCoV.run_consensus_simulations(P::KenyaCoV.CoVParameters_AS,prob,1,cb_regional_lockdown)
P.χ .*= 1.384
sims_controls_scaled = KenyaCoV.run_consensus_simulations(P::KenyaCoV.CoVParameters_AS,prob,1,cb_regional_lockdown)

@save joinpath(homedir(),"Github/KenyaCoVOutputs/sims_consensus_control.jld2") sims_controls
@save joinpath(homedir(),"Github/KenyaCoVOutputs/sims_consensus_control.jld2") sims_controls_scaled
