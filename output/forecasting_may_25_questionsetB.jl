push!(LOAD_PATH, joinpath(homedir(),"GitHub/KenyaCoV/src"))
using Plots,Parameters,Distributions,DifferentialEquations,JLD2,DataFrames,StatsPlots,FileIO,MAT,RecursiveArrayTools,CSV,Dates
import KenyaCoV
# include("/users/Ojal/Documents/Covid-19/jl_models/src/KenyaCoV.jl")
using LinearAlgebra:eigen
using Statistics: median, quantile


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
T_regional_lockdown[:,Nairobi_index] .*= 0.1;T_regional_lockdown[Nairobi_index,Nairobi_index] += 1 - sum(T_regional_lockdown[:,Nairobi_index])
#Mombasa
T_regional_lockdown[:,Mombassa_index] .*= 0.1;T_regional_lockdown[Mombassa_index,Mombassa_index] += 1 - sum(T_regional_lockdown[:,Mombassa_index])
#Kilifi
T_regional_lockdown[:,Kilifi_index] .*= 0.1;T_regional_lockdown[Kilifi_index,Kilifi_index] += 1 - sum(T_regional_lockdown[:,Kilifi_index])
#Kwale
T_regional_lockdown[:,Kwale_index] .*= 0.1;T_regional_lockdown[Kwale_index,Kwale_index] += 1 - sum(T_regional_lockdown[:,Kwale_index])


#Incoming travel
for leaving_area in 1:47,arriving_area in 1:47
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
@load "data/agemixingmatrix_china.jld2" M_China

candidate_M_school = zeros(size(M_Kenya_school))
primary_M_school = zeros(size(M_Kenya_school))
primary_M_school = zeros(size(M_Kenya_school))
secondary_M_school = zeros(size(M_Kenya_school)) #14-18
tertiary_M_school = zeros(size(M_Kenya_school)) #19-22


candidate_M_school[:,3] = 0.2*M_Kenya_school[:,3]#Standard 8
candidate_M_school[:,4] = 0.2*M_Kenya_school[:,4]#Form 4
primary_M_school[:,1:3] = M_Kenya_school[:,1:3]
secondary_M_school[:,3] = 0.2*M_Kenya_school[:,3]
secondary_M_school[:,4] = 0.8*M_Kenya_school[:,4]
tertiary_M_school[:,4] = 0.2*M_Kenya_school[:,4]
tertiary_M_school[:,5] = 0.6*M_Kenya_school[:,5]



include("intervention_combinations.jl");

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
#
# P.β = rand(KenyaCoV.d_R₀) #Choose R₀ randomly from 2-3 95% PI range
# P.c_t = t -> 1.
# P.lockdown = false
# P.schools_closed = false
# P.M = M_Kenya
#
# prob = KenyaCoV.create_KenyaCoV_non_neg_prob(u0,(0.,1*658.),P)
#
# sims_baseline = KenyaCoV.run_consensus_simulations(P,prob,1000,CallbackSet())
#
# # @save joinpath(pwd(),"KenyaCoVOutputs/sims_consensus_baseline_vs2.jld2") sims_baseline

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

@save joinpath(homedir(),"Github/KenyaCoVOutputs/sims_consensus_end_lockdown_30_perc_reduction.jld2") sims_end_regional_lockdown
println("Finished regional lockdown ending. Schools stay shut sims for consensus modelling ")


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

@save joinpath(homedir(),"Github/KenyaCoVOutputs/sims_consensus_open_schools_june_30_perc_reduction.jld2") sims_open_schools_june

println("Finished Schools reopen in June sims for consensus modelling ")

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

@save joinpath(homedir(),"Github/KenyaCoVOutputs/sims_consensus_open_schools_august_30_perc_reduction.jld2") sims_open_schools_august

println("Finished Schools reopen in August sims for consensus modelling ")

"""
What is the effect of reopening schools for Standard 8 and Form 4 candidates only?

Scenario IV: Candidates return to school June 2nd (Term 2) and August 31st (Term 3) and all other classes closed for rest of year. Candidates in school for only two terms (Term 2 and 3).

Scenario IVa:  For the candidates in school their school contacts are assumed to be at 90%.

Scenario IVb:  For the candidates in school their school contacts are assumed to be at 50%.

Scenario V: Candidates return to school August 31st and 4th January 2021 with all other classes closed. Candidates in school for two terms (Term 3 2020 and Term 1 2021).

Scenario Va:  For the candidates in school their school contacts are assumed to be at 90%.

Scenario Vb:  For the candidates in school their school contacts are assumed to be at 50%.

NB: After candidates sit exams and close schools the schools remained close for the remainder of the time the model is running.



What is the effect of reopening schools for primary school learners only (ages 0-14)?

Scenario VIa: Opening schools June 2nd and school contacts at 90%

Scenario VIb: Opening schools June 2nd and school contacts at 50%



What is the effect of reopening schools for secondary school learners only (ages 15-18)?

Scenario VIIa: Opening schools June 2nd and school contacts at 90%

Scenario VIIb: Opening schools June 2nd and school contacts at 50%



What is the effect of reopening schools for tertiary school learners only (ages 19-22)?

Scenario VIIIa: Opening schools June 2nd and school contacts at 90%

Scenario VIIIb: Opening schools June 2nd and school contacts at 50%

 """
