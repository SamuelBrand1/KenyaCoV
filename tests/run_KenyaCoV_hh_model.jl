push!(LOAD_PATH, joinpath(homedir(),"GitHub/KenyaCoV/src"))

using Plots,Parameters,Distributions,DifferentialEquations,JLD2,DataFrames,CSV,RecursiveArrayTools,DelimitedFiles,SparseArrays
using Revise
import KenyaCoV
using LinearAlgebra:eigen,normalize

function transmissionrateinHH(ϵ,γ,σ₂,SAR)
        return (-(ϵ*γ + σ₂) + sqrt((ϵ*γ + σ₂)^2 + 4*ϵ*σ₂*γ*(SAR/(1-SAR))))/(2*ϵ)
end

function eff_popsize!(P::KenyaCoV.CoVParameters_HH,immobile_age_indices)
    P.N̂ = P.T*P.N
    P.N̂[:,immobile_age_indices] .= P.N_region_age[:,immobile_age_indices]



end

n_s = 12
n_ta = 8
"""
function model_ingredients_from_data(filename)

model_ingredients_from_data reads in a JLD2 file which includes all the relevant information required for simulation.
"""
function model_ingredients_from_data(datafile)
    JLD2.@load(datafile,
            N_region_age,
            M_Kenya,
            movements_per_person,
            P_dest,
            ρ,T,σ,rel_detection_rates,
            M_Kenya_ho,
            hosp_rate_by_age,
            ICU_rate_by_age_cond_hosp)

    n,n_a = size(N_region_age)
    index_as = CartesianIndices((1:n, 1:n_a,1:n_s))
    index_as_events = CartesianIndices((1:n, 1:n_a,1:n_ta))
    linear_as= LinearIndices((1:n, 1:n_a,1:n_s))
    linear_as_events = LinearIndices((1:n, 1:n_a,1:n_ta))
    #Population state array
    suspop_kenya = zeros(Int64,n,n_a,n_s) #Array by area, age group and disease state
    for i = 1:n,j=1:n_a
        suspop_kenya[i,j,1] = N_region_age[i,j]
    end

    #Define the change matrix
    dc = sparse(zeros(Int64,n*n_a*n_s,n*n_a*n_ta))
    KenyaCoV.change_matrix!(dc,index_as_events,linear_as)
    # change_matrix(dc)
    #
    # #Parameter definition
    P = KenyaCoV.CoVParameters_HH(n_age = n_a,n_area=n,n_state = n_s,
                        T = T,ρ = ρ,χ = σ,symptomatic_rate = rel_detection_rates[:,2], #This assumes that the detection rates match ϵ = 0.1
                        M = M_Kenya,
                        M_ho = [M_Kenya_ho for i = 1:47],
                        hₐ = hosp_rate_by_age,
                        ICUₐ = ICU_rate_by_age_cond_hosp,
                        dc=dc,
                        N̂=N_region_age,
                        Î = zeros(n,n_a),
                        λ = zeros(n,n_a),
                        λ_loc = zeros(n,n_a),
                        dN= zeros(Int64,n*n_a*n_ta),
                        poi_rates = zeros(n*n_a*n_ta),
                        du_linear = zeros(Int64,n*n_a*n_s),
                        index_as = index_as,
                        index_as_events = index_as_events,
                        linear_as= linear_as,
                        linear_as_events = linear_as_events)

    return suspop_kenya,P,P_dest
end

suspop_kenya,P,P_dest = model_ingredients_from_data("data/data_for_age_structuredmodel_with_counties.jld2")


P.βᵢ = transmissionrateinHH(P.ϵ,P.γ,P.σ₂,0.2)


"""
Load age structured data
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

@load "data/agemixingmatrix_Kenya_all_types.jld2" M_Kenya M_Kenya_ho M_Kenya_other M_Kenya_school M_Kenya_work

counties = CSV.read("data/2019_census_age_pyramids_counties.csv")
names = counties.county
Nairobi_index = findfirst(counties.county .== "Nairobi")
Mombassa_index = findfirst(counties.county .== "Mombasa")
Kwale_index = findfirst(counties.county .== "Kwale")
Kilifi_index = findfirst(counties.county .== "Kilifi")
Mandera_index  = findfirst(counties.county .== "Mandera")

χ_zhang = vcat(0.34*ones(3),ones(10),1.47*ones(4))

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



eigs_kenya_schools_closed, = eigen(sus_matrix.*(1.2*M_Kenya_ho .+ 0.55*M_Kenya_other .+ 0.55*M_Kenya_work).*inf_matrix)
max_eigval_Kenya = Real(eigs_kenya_schools_closed[end])
P.χ .= χ_zhang ./max_eigval_Kenya #This rescales everything so β is the same as R₀ for China

u0[Nairobi_index,8,3] = 500 #10 initial pre-symptomatics in Nairobi
u0[Mombassa_index,8,3] = 300 #10 initial pre-symptomatics in Mombasa
u0[Mandera_index,8,3] = 200 #5 initial pre-symptomatics in Mandera

P.dt = 0.25

@load "data/posterior_distribution_R0.jld2" posterior_R₀



"""
Run model

"""
P.β = rand(posterior_R₀)

P.M = 1.2*M_Kenya_ho .+ 0.55*M_Kenya_other .+ 0.55*M_Kenya_work

prob = KenyaCoV.create_KenyaCoV_non_neg_prob(u0,(0.,60.),P)
# @time sol = solve(prob,FunctionMap(),dt = P.dt)
#
# sims_test = KenyaCoV.run_simulations(P,prob,10)
# output = extract_information_from_simulations(sims_test);
model_str =
"""
This is a test String
for model description.
"""

sims = KenyaCoV.run_simulations(P,prob,10;interventions=CallbackSet())
output = KenyaCoV.extract_information_from_simulations(sims);
scenariodata = generate_report(output,model_str,"_test"," (test)",counties.county;make_new_directory=false);
scenariodata.prevalence_ICU_ts.med
scenariodata = KenyaCoV.run_scenario(P,prob,10,model_str,"_test"," (test)",counties.county;interventions = CallbackSet(),make_new_directory=true)


data = KenyaCoV.run_scenario(P,prob,10,model_str,"_test"," (test)",counties.county;interventions = CallbackSet(),make_new_directory = false)
