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
# χ_zhang = vcat(0.34*ones(3),ones(10),1.47*ones(4))

P.χ = ones(17)
P.rel_detection_rate = d_1
P.dt = 0.25;
P.ext_inf_rate = 0.;
P.ϵ = 1.
#Set the susceptibility vector --- just to specify the correct R₀
sus_matrix = repeat(ones(17),1,17)
R_A = P.ϵ*((1/P.σ₂) + (1/P.γ) ) #effective duration of asymptomatic
R_M = (P.ϵ/P.σ₂) + (P.ϵ_D/P.γ) #effective duration of mild
R_V = (P.ϵ/P.σ₂) + (P.ϵ_V/P.τ) #effective duration of severe
R_vector = [(1-P.rel_detection_rate[a])*R_A + P.rel_detection_rate[a]*(1-P.hₐ[a])*R_M + P.rel_detection_rate[a]*P.hₐ[a]*R_V for a = 1:17]
inf_matrix = repeat(R_vector',17,1)


eigs_kenya_schools_closed, = eigen(sus_matrix.*(1.2*M_Kenya_ho .+ 0.55*M_Kenya_other .+ 0.55*M_Kenya_work).*inf_matrix)
max_eigval_Kenya = Real(eigs_kenya_schools_closed[end])
P.χ .= ones(17) ./max_eigval_Kenya #This rescales everything so β is the same as R₀ for China

u0[Nairobi_index,8,3] = 500 #10 initial pre-symptomatics in Nairobi
u0[Mombassa_index,8,3] = 300 #10 initial pre-symptomatics in Mombasa
u0[Mandera_index,8,3] = 200 #5 initial pre-symptomatics in Mandera

P.dt = 0.25
"""
Scenario IV: Candidates return to school June 2nd (Term 2) and August 31st (Term 3) and all other classes closed for rest of year. Candidates in school for only two terms (Term 2 and 3).

Scenario IVa:  For the candidates in school their school contacts are assumed to be at 90%.

Scenario IVb:  For the candidates in school their school contacts are assumed to be at 50%.
"""
#
model_str_90 =
"""
Scenario IVa:  For the candidates in school their school contacts are assumed to be at 90%.
"""
model_str_50 =
"""
Scenario IVb:  For the candidates in school their school contacts are assumed to be at 50%.
"""

P.β = 1.22
P.c_t = t -> 1.
P.lockdown = false
P.schools_closed = false
P.before_week_two = true

P.M =0.8*M_Kenya_ho .+ M_Kenya_other .+ M_Kenya_work .+ 0.5*M_Kenya_school

prob = KenyaCoV.create_KenyaCoV_non_neg_prob(u0,(0.,1*658.),P)


scenariodata = KenyaCoV.run_scenario(P,prob,200,model_str_90,"_uniform_sus_candidates_june_90perc"," (candidates only 90%)",
                                    counties.county;
                                    interventions = measures_schools_open_june_2020_90pct_candidatesonly,
                                    make_new_directory=true);
#reset contact matrix and other variables
P.M =0.8*M_Kenya_ho .+ M_Kenya_other .+ M_Kenya_work .+ 0.5*M_Kenya_school
P.β = 1.22
P.c_t = t -> 1.
P.lockdown = false
P.schools_closed = false
P.before_week_two = true

scenariodata = KenyaCoV.run_scenario(P,prob,200,model_str_50,"_uniform_sus_candidates_june_50perc"," (candidates only 50%)",
                                    counties.county;
                                    interventions = measures_schools_open_june_2020_50pct_candidatesonly,
                                    make_new_directory=true);

println("Completed scenario (uniform susceptibility)  IV")

"""
Scenario V: Candidates return to school August 31st and 4th January 2021 with all other classes closed. Candidates in school for two terms (Term 3 2020 and Term 1 2021).

Scenario Va:  For the candidates in school their school contacts are assumed to be at 90%.

Scenario Vb:  For the candidates in school their school contacts are assumed to be at 50%.
"""

#
model_str_90 =
"""
Scenario Va:  For the candidates in school their school contacts are assumed to be at 90%.
"""
model_str_50 =
"""
Scenario Vb:  For the candidates in school their school contacts are assumed to be at 50%.
"""

P.β = 1.22
P.c_t = t -> 1.
P.lockdown = false
P.schools_closed = false
P.before_week_two = true

P.M =0.8*M_Kenya_ho .+ M_Kenya_other .+ M_Kenya_work .+ 0.5*M_Kenya_school

prob = KenyaCoV.create_KenyaCoV_non_neg_prob(u0,(0.,1*658.),P)


scenariodata = KenyaCoV.run_scenario(P,prob,200,model_str_90,"_uniform_sus_candidates_august_90perc"," (candidates only 90%)",
                                    counties.county;
                                    interventions = measures_schools_open_august_2020_90pct_candidatesonly,
                                    make_new_directory=true);
#reset contact matrix and other variables
P.M =0.8*M_Kenya_ho .+ M_Kenya_other .+ M_Kenya_work .+ 0.5*M_Kenya_school
P.β = 1.22
P.c_t = t -> 1.
P.lockdown = false
P.schools_closed = false
P.before_week_two = true

scenariodata = KenyaCoV.run_scenario(P,prob,200,model_str_50,"_uniform_sus_candidates_august_50perc"," (candidates only 50%)",
                                    counties.county;
                                    interventions = measures_schools_open_august_2020_50pct_candidatesonly,
                                    make_new_directory=true);

println("Completed scenario (uniform susceptibility)  V")



"""
What is the effect of reopening schools for primary school learners only (ages 0-14)?

Scenario VIa: Opening schools June 2nd and school contacts at 90%

Scenario VIb: Opening schools June 2nd and school contacts at 50%

"""


#
model_str_90 =
"""
Scenario VIa: Opening schools June 2nd and school contacts at 90%
"""
model_str_50 =
"""
Scenario VIb: Opening schools June 2nd and school contacts at 50%
"""

P.β = 1.22
P.c_t = t -> 1.
P.lockdown = false
P.schools_closed = false
P.before_week_two = true

P.M =0.8*M_Kenya_ho .+ M_Kenya_other .+ M_Kenya_work .+ 0.5*M_Kenya_school

prob = KenyaCoV.create_KenyaCoV_non_neg_prob(u0,(0.,1*658.),P)


scenariodata = KenyaCoV.run_scenario(P,prob,200,model_str_90,"_uniform_sus_primary_june_90perc"," (Primary only 90%)",
                                    counties.county;
                                    interventions = measures_schools_open_june_2020_90pct_primaryonly,
                                    make_new_directory=true);
#reset contact matrix and other variables
P.M =0.8*M_Kenya_ho .+ M_Kenya_other .+ M_Kenya_work .+ 0.5*M_Kenya_school
P.β = 1.22
P.c_t = t -> 1.
P.lockdown = false
P.schools_closed = false
P.before_week_two = true

scenariodata = KenyaCoV.run_scenario(P,prob,200,model_str_50,"_uniform_sus_primary_june_50perc"," (Primary only 50%)",
                                    counties.county;
                                    interventions = measures_schools_open_june_2020_50pct_primaryonly,
                                    make_new_directory=true);

println("Completed scenario (uniform susceptibility)  VI")

"""
What is the effect of reopening schools for secondary school learners only (ages 15-18)?

Scenario VIIa: Opening schools June 2nd and school contacts at 90%

Scenario VIIb: Opening schools June 2nd and school contacts at 50%

"""

model_str_90 =
"""
Scenario VIIa: Opening schools June 2nd and school contacts at 90%
"""
model_str_50 =
"""
Scenario VIIb: Opening schools June 2nd and school contacts at 50%
"""

P.β = 1.22
P.c_t = t -> 1.
P.lockdown = false
P.schools_closed = false
P.before_week_two = true

P.M =0.8*M_Kenya_ho .+ M_Kenya_other .+ M_Kenya_work .+ 0.5*M_Kenya_school

prob = KenyaCoV.create_KenyaCoV_non_neg_prob(u0,(0.,1*658.),P)


scenariodata = KenyaCoV.run_scenario(P,prob,200,model_str_90,"_uniform_sus_secondary_june_90perc"," (Secondary only 90%)",
                                    counties.county;
                                    interventions = measures_schools_open_june_2020_90pct_secondaryonly,
                                    make_new_directory=true);
#reset contact matrix and other variables
P.M =0.8*M_Kenya_ho .+ M_Kenya_other .+ M_Kenya_work .+ 0.5*M_Kenya_school
P.β = 1.22
P.c_t = t -> 1.
P.lockdown = false
P.schools_closed = false
P.before_week_two = true

scenariodata = KenyaCoV.run_scenario(P,prob,200,model_str_50,"_uniform_sus_secondary_june_50perc"," (Secondary only 50%)",
                                    counties.county;
                                    interventions = measures_schools_open_june_2020_50pct_secondaryonly,
                                    make_new_directory=true);

println("Completed scenario (uniform susceptibility)  VII")


"""
What is the effect of reopening schools for tertiary school learners only (ages 19-22)?

Scenario VIIIa: Opening schools June 2nd and school contacts at 90%

Scenario VIIIb: Opening schools June 2nd and school contacts at 50%

"""

model_str_90 =
"""
Scenario VIIIa: Opening schools June 2nd and school contacts at 90%
"""
model_str_50 =
"""
Scenario VIIIb: Opening schools June 2nd and school contacts at 50%
"""

P.β = 1.22
P.c_t = t -> 1.
P.lockdown = false
P.schools_closed = false
P.before_week_two = true

P.M =0.8*M_Kenya_ho .+ M_Kenya_other .+ M_Kenya_work .+ 0.5*M_Kenya_school

prob = KenyaCoV.create_KenyaCoV_non_neg_prob(u0,(0.,1*658.),P)


scenariodata = KenyaCoV.run_scenario(P,prob,200,model_str_90,"_uniform_sus_tertiary_june_90perc"," (Tertiary only 90%)",
                                    counties.county;
                                    interventions = measures_schools_open_june_2020_90pct_tertiaryonly,
                                    make_new_directory=true);
#reset contact matrix and other variables
P.M =0.8*M_Kenya_ho .+ M_Kenya_other .+ M_Kenya_work .+ 0.5*M_Kenya_school
P.β = 1.22
P.c_t = t -> 1.
P.lockdown = false
P.schools_closed = false
P.before_week_two = true

scenariodata = KenyaCoV.run_scenario(P,prob,200,model_str_50,"_uniform_sus_tertiary_june_50perc"," (Tertiary only 50%)",
                                    counties.county;
                                    interventions = measures_schools_open_june_2020_50pct_tertiaryonly,
                                    make_new_directory=true);

println("Completed scenario (uniform susceptibility)  VIII")


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
