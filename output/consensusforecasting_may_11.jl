push!(LOAD_PATH, joinpath(homedir(),"/Documents/Covid-19/jl_models/src"))
using Plots,Parameters,Distributions,DifferentialEquations,JLD2,DataFrames,StatsPlots,FileIO,MAT,RecursiveArrayTools,CSV
using Revise
include("/users/Ojal/Documents/Covid-19/jl_models/src/KenyaCoV.jl")
using LinearAlgebra:eigen
using Statistics: median, quantile

"""
Consensus modelling --- May week 1

1) estimate epidemic spread by county
2) estimate the effect of reopening schools on either 2 June 2020 or delaying opening to 31 August 2020
3) the effect of lifting travel restrictions between counties in [16th] May 2020. 
"""

"""
Load age structured data, define callback control measures, and effect of regional lockdown
"""

u0,P,P_dest = KenyaCoV.model_ingredients_from_data("data/data_for_age_structuredmodel_with_counties.jld2",
                                            "data/flight_numbers.csv",
                                            "data/projected_global_prevelance.csv")
counties = CSV.read("data/2019_census_age_pyramids_counties.csv")
Nairobi_index = findfirst(counties.county .== "Nairobi")
Mombassa_index = findfirst(counties.county .== "Mombasa")
Kwale_index = findfirst(counties.county .== "Kwale")
Kilifi_index = findfirst(counties.county .== "Kilifi")
Mandera_index  = findfirst(counties.county .== "Mandera")


#Put in the lockdown effects
T_normal = deepcopy(P.T)
T_regional_lockdown = deepcopy(P.T)

#Outgoing travel
#Nairobi
T_regional_lockdown[:,Nairobi_index] *= 0.1;T_regional_lockdown[Nairobi_index,Nairobi_index] += 1 - sum(T_regional_lockdown[:,Nairobi_index])
#Mombasa
T_regional_lockdown[:,Mombassa_index] *= 0.1;T_regional_lockdown[Mombassa_index,Mombassa_index] += 1 - sum(T_regional_lockdown[:,Mombassa_index])
#Kilifi
T_regional_lockdown[:,Kilifi_index] *= 0.1;T_regional_lockdown[Kilifi_index,Kilifi_index] += 1 - sum(T_regional_lockdown[:,Kilifi_index])
#Kwale
T_regional_lockdown[:,Kwale_index] *= 0.1;T_regional_lockdown[Kwale_index,Kwale_index] += 1 - sum(T_regional_lockdown[:,Kwale_index])


#Incoming travel
for leaving_area in 1:47,arriving_area in 1:47
    if !(leaving_area in [Nairobi_index,Mombassa_index,Kwale_index,Kilifi_index]) && arriving_area in [Nairobi_index,Mombassa_index,Kwale_index,Kilifi_index]
        amount_reduced = 0.9*T_regional_lockdown[arriving_area,leaving_area]
        T_regional_lockdown[arriving_area,leaving_area] -= amount_reduced
        T_regional_lockdown[leaving_area,leaving_area] += amount_reduced #All avoided trips to locked down areas lead to staying at home
    end
end


@load "data/inference_for_age_dependent_symptomatic_rates/detection_rates_for_different_epsilons_model2.jld2" d_0 d_01 d_025 d_05 d_1
χ_zhang = vcat(0.34*ones(3),ones(10),1.47*ones(4))

#Initial infecteds
"""
Load age mixing matrices (these are all in to (row) from (col) format)
"""

@load "data/agemixingmatrix_Kenya_all_types.jld2" M_Kenya M_Kenya_ho M_Kenya_other M_Kenya_school M_Kenya_work
@load "data/agemixingmatrix_china.jld2" M_China

#Function for changing contacts so as to have -45% over 14 days
function ramp_down(t)
    if t < 0.
        return 1.
    elseif t>= 0. && t <= 14.
        return (1-t/14.) + 0.5*t/14.
    elseif t > 14.
        return 0.5
    end
end
using Dates
Date(2020,4,7) - Date(2020,3,13) # lockdown 7th April 2020
Date(2020,6,2) - Date(2020,3,13) # School open 2nd June 2020
Date(2020,8,14) - Date(2020,3,13) # Schools closed 14th Augus 2020

Date(2020,8,31) - Date(2020,3,13) # Schools  open 31st august 2020
Date(2020,5,16) - Date(2020,3,13) # Regional lockdown end 16th May
Date(2020,10,30) - Date(2020,3,13) # Schools closed 30th October 2020

Date(2021,1,4) - Date(2020,3,13)  # schools open 4th Jan 2021
Date(2021,4,9) - Date(2020,3,13)  # Schools closed 9th April 2021
Date(2021,5,3) - Date(2020,3,13)  # Schools open 3rd may 2021
Date(2021,8,6) - Date(2020,3,13)  # Schools closed 6th August 2021
Date(2021,8,30) - Date(2020,3,13)  # Schools open 30th August 2021
Date(2021,10,22) - Date(2020,3,13)  # Schools closed 22nd October 2021


Date(2021,12,31) - Date(2020,3,13)  # we now run intil Dec 2021


#Put in the regional lockdown on day 25 (April 7th)
function regional_lockdown_timing(u,t,integrator)
  !integrator.p.lockdown && t > 25.
end
function affect_regional_lockdown!(integrator)
  integrator.p.T = T_regional_lockdown
  integrator.p.lockdown = true
end
cb_regional_lockdown = DiscreteCallback(regional_lockdown_timing,affect_regional_lockdown!)

#End the regional lockdown on day 64 (May 16th)
function regional_lockdown_ending(u,t,integrator)
  integrator.p.lockdown && t > 64.
end
function affect_regional_lockdown_end!(integrator)
  integrator.p.T = T_normal
  integrator.p.lockdown = false
end
cb_regional_lockdown_end = DiscreteCallback(regional_lockdown_ending,affect_regional_lockdown_end!)

regional_lockdown_starts = CallbackSet(cb_regional_lockdown)
regional_lockdown_starts_and_finishes = CallbackSet(cb_regional_lockdown,cb_regional_lockdown_end)

#Closure and opening of schools
function open_schools_june(u,t,integrator)
  integrator.p.schools_closed && t > 81.
end

function close_schools_august(u,t,integrator)
  !integrator.p.schools_closed && t > 154.
end

function open_schools_august(u,t,integrator)
  integrator.p.schools_closed && t > 171.
end

function close_schools_october(u,t,integrator)
  !integrator.p.schools_closed && t > 231.
end


function open_schools_jan2021(u,t,integrator)
  integrator.p.schools_closed && t > 297.
end

function close_schools_apr2021(u,t,integrator)
  !integrator.p.schools_closed && t > 392.
end

function open_schools_may2021(u,t,integrator)
  integrator.p.schools_closed && t > 416.
end

function close_schools_aug2021(u,t,integrator)
  !integrator.p.schools_closed && t > 511.
end

function open_schools_aug2021(u,t,integrator)
  integrator.p.schools_closed && t > 535.
end

function close_schools_oct2021(u,t,integrator)
  !integrator.p.schools_closed && t > 588.
end

function affect_open_schools!(integrator)
  integrator.p.M = 1.1*M_Kenya_ho .+ M_Kenya_other .+ M_Kenya_work .+ 0.7*M_Kenya_school
  integrator.p.schools_closed = false
end
function affect_close_schools!(integrator)
  integrator.p.M = 1.2*M_Kenya_ho .+ M_Kenya_other .+ M_Kenya_work
  integrator.p.schools_closed = true
end

cb_open_schools_june = DiscreteCallback(open_schools_june,affect_open_schools!)
cb_open_schools_august = DiscreteCallback(open_schools_august,affect_open_schools!)
cb_open_schools_jan2021 = DiscreteCallback(open_schools_jan2021,affect_open_schools!)
cb_open_schools_may2021 = DiscreteCallback(open_schools_may2021,affect_open_schools!)
cb_open_schools_aug2021 = DiscreteCallback(open_schools_aug2021,affect_open_schools!)

cb_close_schools_august = DiscreteCallback(close_schools_august,affect_close_schools!)
cb_close_schools_october = DiscreteCallback(close_schools_october,affect_close_schools!)
cb_close_schools_apr2021 = DiscreteCallback(close_schools_apr2021,affect_close_schools!)
cb_close_schools_aug2021 = DiscreteCallback(close_schools_aug2021,affect_close_schools!)
cb_close_schools_oct2021 = DiscreteCallback(close_schools_oct2021,affect_close_schools!)

measures_schools_open_june_2020 = CallbackSet(cb_regional_lockdown,
                                            cb_open_schools_june,
                                            cb_close_schools_august,
                                            cb_open_schools_august,
                                            cb_close_schools_october,
                                            cb_open_schools_jan2021,
                                            cb_close_schools_apr2021,
                                            cb_open_schools_may2021,
                                            cb_close_schools_aug2021,
                                            cb_open_schools_aug2021,
                                            cb_close_schools_oct2021)
measures_schools_open_august_2020 = CallbackSet(cb_regional_lockdown,
                                            cb_open_schools_august,
                                            cb_close_schools_october,
                                            cb_open_schools_jan2021,
                                            cb_close_schools_apr2021,
                                            cb_open_schools_may2021,
                                            cb_close_schools_aug2021,
                                            cb_open_schools_aug2021,
                                            cb_close_schools_oct2021)

"""
Set up parameters
"""



u0,P,P_dest = KenyaCoV.model_ingredients_from_data("data/data_for_age_structuredmodel_with_counties.jld2",
                                            "data/flight_numbers.csv",
                                            "data/projected_global_prevelance.csv")
#Redistribute susceptibility JUST to rescale infectiousness so we get the correct R₀/r
P.χ = copy(χ_zhang)
P.rel_detection_rate = d_1
P.dt = 0.25;
P.ext_inf_rate = 0.;
P.ϵ = 1.
#Set the susceptibility vector --- just to specify the correct R₀
sus_matrix = repeat(χ_zhang,1,17)
R_A = P.ϵ*((1/P.σ₂) + (1/P.γ) ) #effective duration of asymptomatic
R_M = (P.ϵ/P.σ₂) + (P.ϵ_D/P.γ) #effective duration of mild
R_V = (P.ϵ/P.σ₂) + (P.ϵ_V/P.τ) #effective duration of severe
R_vector = [(1-P.rel_detection_rate[a])*R_A + P.rel_detection_rate[a]*(1-P.hₐ[a])*R_M + P.rel_detection_rate[a]*P.hₐ[a]*R_V for a = 1:17]
inf_matrix = repeat(R_vector',17,1)

eigs_china, = eigen(sus_matrix.*M_China.*inf_matrix)
max_eigval_china = Real(eigs_china[end])
eigs_kenya, = eigen(sus_matrix.*M_Kenya.*inf_matrix)
max_eigval_Kenya = Real(eigs_kenya[end])
multiplier_for_kenya = max_eigval_Kenya/max_eigval_china
P.χ .= χ_zhang ./max_eigval_china #This rescales everything so β is the same as R₀ for China

u0[Nairobi_index,8,3] = 30 #10 initial pre-symptomatics in Nairobi
u0[Mombassa_index,8,3] = 10 #10 initial pre-symptomatics in Mombasa
u0[Mandera_index,8,3] = 5 #5 initial pre-symptomatics in Mandera

"""
SCENARIO 1

Base line scenario
"""

P.β = rand(KenyaCoV.d_R₀) #Choose R₀ randomly from 2-3 95% PI range
P.c_t = t -> 1.
P.lockdown = false
P.schools_closed = false
P.M = M_Kenya

prob = KenyaCoV.create_KenyaCoV_non_neg_prob(u0,(0.,1*658.),P)

sims_baseline = KenyaCoV.run_consensus_simulations(P,prob,1000,CallbackSet())

@save joinpath(pwd(),"KenyaCoVOutputs/sims_consensus_baseline_vs2.jld2") sims_baseline

"""
SCENARIO 2

Full intervention scenario:
- Social mixing school 0%, home 120%, 100% other / work –
- Decline by 50% over two weeks
- Nbi / coast lockdown 90% applies 17th April
"""

P.β = rand(KenyaCoV.d_R₀) #Choose R₀ randomly from 2-3 95% PI range
P.c_t = ramp_down
P.lockdown = false
P.schools_closed = true
P.M = 1.2*M_Kenya_ho .+ M_Kenya_other .+ M_Kenya_work

prob = KenyaCoV.create_KenyaCoV_non_neg_prob(u0,(0.,1*658.),P)

sims_full_intervention = KenyaCoV.run_consensus_simulations(P,prob,1000,regional_lockdown_starts)

@save joinpath(pwd(),"KenyaCoVOutputs/sims_consensus_full_intervention_vs2.jld2") sims_full_intervention


#"""
#SCENARIO 3

#Regional lockdown ending. Schools stay shut
#"""

#P.β = rand(KenyaCoV.d_R₀) #Choose R₀ randomly from mean 2.5 (2-3) range
#P.c_t = ramp_down #This implements the social distancing over 14 days from time 0.

#P.lockdown = false
#P.schools_closed = true
#P.M = 1.2*M_Kenya_ho .+ M_Kenya_other .+ M_Kenya_work

#prob = KenyaCoV.create_KenyaCoV_non_neg_prob(u0,(0.,1*658.),P)

#sims_end_regional_lockdown = KenyaCoV.run_consensus_simulations(P::KenyaCoV.CoVParameters_AS,prob,1000,regional_lockdown_starts_and_finishes)

#@save joinpath(pwd(),"KenyaCoVOutputs/sims_consensus_end_lockdown.jld2") sims_end_regional_lockdown


"""
SCENARIO 3

Schools reopen in June
"""

P.β = rand(KenyaCoV.d_R₀) #Choose R₀ randomly from mean 2.5 (2-3) range
P.c_t = ramp_down #This implements the social distancing over 14 days from time 0.

P.lockdown = false
P.schools_closed = true
P.M = 1.2*M_Kenya_ho .+ M_Kenya_other .+ M_Kenya_work

prob = KenyaCoV.create_KenyaCoV_non_neg_prob(u0,(0.,1*658.),P)

sims_open_schools_june = KenyaCoV.run_consensus_simulations(P::KenyaCoV.CoVParameters_AS,prob,1000,measures_schools_open_june_2020)

@save joinpath(pwd(),"KenyaCoVOutputs/sims_consensus_open_schools_june.jld2") sims_open_schools_june


"""
SCENARIO 4
Schools reopen in August
"""

P.β = rand(KenyaCoV.d_R₀) #Choose R₀ randomly from mean 2.5 (2-3) range
P.c_t = ramp_down #This implements the social distancing over 14 days from time 0.

P.lockdown = false
P.schools_closed = true
P.M = 1.2*M_Kenya_ho .+ M_Kenya_other .+ M_Kenya_work

prob = KenyaCoV.create_KenyaCoV_non_neg_prob(u0,(0.,1*658.),P)

sims_open_schools_august = KenyaCoV.run_consensus_simulations(P::KenyaCoV.CoVParameters_AS,prob,1000,measures_schools_open_august_2020)

@save joinpath(pwd(),"KenyaCoVOutputs/sims_consensus_open_schools_august.jld2") sims_open_schools_august
