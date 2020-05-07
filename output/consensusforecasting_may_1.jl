push!(LOAD_PATH, joinpath(homedir(),"GitHub/KenyaCoV/src"))
using Plots,Parameters,Distributions,DifferentialEquations,JLD2,DataFrames,StatsPlots,FileIO,MAT,RecursiveArrayTools,CSV
using Revise
import KenyaCoV
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
for leaving_area in 1:20,arriving_area in 1:20
    if !(leaving_area in [Nairobi_index,Mombassa_index,Kwale_index,Kilifi_index]) && arriving_area in [Nairobi_index,Mombassa_index,Kwale_index,Kilifi_index]
        amount_reduced = 0.9*T_regional_lockdown[arriving_area,leaving_area]
        T_regional_lockdown[arriving_area,leaving_area] -= amount_reduced
        T_regional_lockdown[leaving_area,leaving_area] += amount_reduced #All avoided trips to locked down areas lead to staying at home
    end
end




@load "data/detection_rates_for_different_epsilons_model2.jld2" d_0 d_01 d_025 d_05 d_1
χ_zhang = vcat(0.34*ones(3),ones(10),1.47*ones(4))

#Initial infecteds
"""
Load age mixing matrices (these are all in to (row) from (col) format)
"""

@load "data/agemixingmatrix_Kenya_all_types.jld2" M_Kenya M_Kenya_ho M_Kenya_other M_Kenya_school M_Kenya_work


#Function for changing contacts so as to have -45% over 14 days
function ramp_down(t)
    if t < 0.
        return 1.
    elseif t>= 0. && t <= 14.
        return (1-t/14.) + 0.55*t/14.
    elseif t > 14.
        return 0.55
    end
end
using Dates
Date(2020,4,7) - Date(2020,3,13)
Date(2020,6,2) - Date(2020,3,13)
Date(2020,8,31) - Date(2020,3,13)
Date(2020,5,16) - Date(2020,3,13)

#Put in the regional lockdown on day 25 (April 7th)
function regional_lockdown_timing(u,t,integrator)
  !integrator.p.lockdown && t > 25.
end
function affect_regional_lockdown!(integrator)
  integrator.p.T = T_regional_lockdown
  integrator.p.lockdown = true
end
cb_regional_lockdown = DiscreteCallback(regional_lockdown_timing,affect_regional_lockdown!)

#Put in the regional lockdown on day 64 (May 16th)
function regional_lockdown_ending(u,t,integrator)
  integrator.p.lockdown && t > 64.
end
function affect_regional_lockdown_end!(integrator)
  integrator.p.T = T_normal
  integrator.p.lockdown = false
end
cb_regional_lockdown_end = DiscreteCallback(regional_lockdown_ending,affect_regional_lockdown_end!)

regional_lockdown_starts_and_finishes = CallbackSet(cb_regional_lockdown,cb_regional_lockdown_end)
#Open schools on June 2nd or August 31st
function open_schools_june(u,t,integrator)
  integrator.p.schools_closed && t > 81.
end
function open_schools_august(u,t,integrator)
  integrator.p.schools_closed && t > 171.
end
function affect_open_schools!(integrator)
  integrator.p.M = M_Kenya
  integrator.p.schools_closed = false
end
cb_open_schools_june = DiscreteCallback(open_schools_june,affect_open_schools!)
cb_open_schools_august = DiscreteCallback(open_schools_august,affect_open_schools!)


"""
Set up parameters
"""



u0,P,P_dest = KenyaCoV.model_ingredients_from_data("data/data_for_age_structuredmodel_with_counties.jld2",
                                            "data/flight_numbers.csv",
                                            "data/projected_global_prevelance.csv")
#Redistribute susceptibility JUST to rescale infectiousness so we get the correct R₀/r
P.χ = χ_zhang
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
P.χ .*= 1/max_eigval #This rescales everything so β is the same as R₀

u0[Nairobi_index,8,3] = 30 #10 initial pre-symptomatics in Nairobi
u0[Mombassa_index,8,3] = 10 #10 initial pre-symptomatics in Mombasa
u0[Mandera_index,8,3] = 5 #5 initial pre-symptomatics in Mandera

"""
Base line scenario
"""

P.β = rand(KenyaCoV.d_R₀) #Choose R₀ randomly from 2-3 range
P.lockdown = false
P.schools_closed = false
P.M = M_Kenya

prob = KenyaCoV.create_KenyaCoV_non_neg_prob(u0,(0.,1*365.),P)

sims_baseline = KenyaCoV.run_consensus_simulations(P,prob,1000,CallbackSet())

@save joinpath(homedir(),"Github/KenyaCoVOutputs/sims_consensus_baseline_vs2.jld2") sims_baseline


println("Finished baseline sims for consensus modelling ")



"""
SCENARIO 2 --- regional lockdown ending. Schools stay shut
"""



P.β = rand(KenyaCoV.d_R₀) #Choose R₀ randomly from mean 2.5 (2-3) range
P.c_t = ramp_down #This implements the social distancing over 14 days from time 0.

P.lockdown = false
P.schools_closed = true
P.M = M_Kenya_ho .+ M_Kenya_other .+ M_Kenya_work

prob = KenyaCoV.create_KenyaCoV_non_neg_prob(u0,(0.,1*365.),P)

sims_end_regional_lockdown = KenyaCoV.run_consensus_simulations(P::KenyaCoV.CoVParameters_AS,prob,1000,regional_lockdown_starts_and_finishes)

@save joinpath(homedir(),"Github/KenyaCoVOutputs/sims_consensus_end_lockdown.jld2") sims_end_regional_lockdown


"""
SCENARIO 3 --- Schools reopen in June
"""

P.β = rand(KenyaCoV.d_R₀) #Choose R₀ randomly from mean 2.5 (2-3) range
P.c_t = ramp_down #This implements the social distancing over 14 days from time 0.

P.lockdown = false
P.schools_closed = true
P.M = M_Kenya_ho .+ M_Kenya_other .+ M_Kenya_work

prob = KenyaCoV.create_KenyaCoV_non_neg_prob(u0,(0.,1*365.),P)

sims_open_schools_june = KenyaCoV.run_consensus_simulations(P::KenyaCoV.CoVParameters_AS,prob,1000,CallbackSet(cb_regional_lockdown,cb_open_schools_june))

@save joinpath(homedir(),"Github/KenyaCoVOutputs/sims_consensus_open_schools_june.jld2") sims_open_schools_june


"""
SCENARIO 4 --- Schools reopen in August
"""

P.β = rand(KenyaCoV.d_R₀) #Choose R₀ randomly from mean 2.5 (2-3) range
P.c_t = ramp_down #This implements the social distancing over 14 days from time 0.

P.lockdown = false
P.schools_closed = true
P.M = M_Kenya_ho .+ M_Kenya_other .+ M_Kenya_work

prob = KenyaCoV.create_KenyaCoV_non_neg_prob(u0,(0.,1*365.),P)

sims_open_schools_august = KenyaCoV.run_consensus_simulations(P::KenyaCoV.CoVParameters_AS,prob,1000,CallbackSet(cb_regional_lockdown,cb_open_schools_august))

@save joinpath(homedir(),"Github/KenyaCoVOutputs/sims_consensus_open_schools_august.jld2") sims_open_schools_august
