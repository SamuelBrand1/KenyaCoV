#Script for building the population pyramid into data
using DataFrames,Plots,JLD2,MAT,LinearAlgebra,Distances,StatsPlots,CSV
using LinearAlgebra:normalize,normalize!
using TransformVariables, DynamicHMC, DynamicHMC.Diagnostics, Distributions
age_cats = ["0-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39","40-44","45-49","50-54","55-59","60-64","65-69","70-74","75-79","80+"]
# age_cats_in
# Age specific rates of hospitalisation if becoming a case
#Estimate from CDC COVID-19 Response Team. Severe Outcomes Among Patients with Coronavirus Disease 2019 (COVID-19) - United States, February 12-March 16, 2020. MMWR Morb. Mortal. Wkly. Rep. 69, 343–346 (2020).
hosp_rate_CDC = [mean([1.6,2.5]),#0-19 year olds
                 mean([14.3,20.8]),#20-44 yos
                 mean([21.2,28.3]),#45-54
                 mean([20.5,30.1]),#55-64
                 mean([28.6,43.5]),#65-74
                 mean([30.5,58.7]),#75-84
                 mean([31.3,70.3])]./100 #85+
hosp_rate_by_age = [hosp_rate_CDC[1],hosp_rate_CDC[1],hosp_rate_CDC[1],hosp_rate_CDC[1],
                    hosp_rate_CDC[2],hosp_rate_CDC[2],hosp_rate_CDC[2],hosp_rate_CDC[2],hosp_rate_CDC[2],
                    hosp_rate_CDC[3],hosp_rate_CDC[3],
                    hosp_rate_CDC[4],hosp_rate_CDC[4],
                    hosp_rate_CDC[5],hosp_rate_CDC[5],
                    hosp_rate_CDC[6],hosp_rate_CDC[7]]

"""
Sub county population data
"""
# subcounty_population = readtable("data/distribution-of-population-by-age-sex-county-and-sub-county-kenya-2019-census-volume-iii.csv")
# subcounties = unique(subcounty_population.sub_county)
# unique(subcounty_population.Age)
# subC_to_C = Dict{String,String}()
# # for
"""
Load the MCMC outputs for the age-specific relative symptomatic rate
"""
@load "data/MCMC_results_fit_d.jld2" results_d_tau_0 results_d_tau_01 results_d_tau_025 results_d_tau_05 results_d_tau_1
trans = as((θ = as(Array, asℝ₊, 17),))#Transforms the relative values into definitely positive numbers 16th value is fixed as 1
gr()
function MCMCchain2estimates(results)
    posterior = first.(transform.(trans, results.chain))
    posterior_array = zeros(length(posterior),17)
    for i = 1:length(posterior),j=1:17
        posterior_array[i,j] = posterior[i][j]
    end
    for i = 1:length(posterior),j=1:17
        posterior_array[i,j] = posterior_array[i,j]/posterior_array[i,17]
    end
    θ̂ = [mean(posterior_array[:,i]) for i = 1:17]
    lb = [quantile(posterior_array[:,i],0.025) for i = 1:17]
    ub = [quantile(posterior_array[:,i],0.975) for i = 1:17]
    return θ̂,θ̂ .- lb,ub.-θ̂
end

d_0,lb_0,ub_0 = MCMCchain2estimates(results_d_tau_0)
d_01,lb_01,ub_01 = MCMCchain2estimates(results_d_tau_01)
d_025,lb_025,ub_025 = MCMCchain2estimates(results_d_tau_025)
d_05,lb_05,ub_05 = MCMCchain2estimates(results_d_tau_05)
d_1,lb_1,ub_1 = MCMCchain2estimates(results_d_tau_1)

# fig_d = groupedbar(hcat(θ̂_0,θ̂_01,θ̂_025,θ̂_05,θ̂_1),fillalpha = 0.7,yscale = :log10)
# fig_d = scatter(θ̂_0,fillalpha = 0.7,yscale = :log10)


fig_d = scatter(d_0,yerr = (lb_0,ub_0),fillalpha = 0.3,lab="rel. infectiousness = 0%",legend = :bottomright,
        xticks = (1:17,age_cats),
        ylabel = "Relative symptomatic rate", title = "Symptomatic rates for subclinical rel. infectiousness scenarios",
        yscale = :log10,
        ms = 10)
scatter!(fig_d,d_01,yerr = (lb_01,ub_01),fillalpha = 0.3,lab="rel. infectiousness = 10%",ms=10)
scatter!(fig_d,d_025,yerr = (lb_025,ub_025),fillalpha = 0.3,lab="rel. infectiousness = 25%",ms=10)
scatter!(fig_d,d_05,yerr = (lb_05,ub_05),fillalpha = 0.3,lab="rel. infectiousness = 50%",ms=10)
scatter!(fig_d,d_1,yerr = (lb_1,ub_1),fillalpha = 0.3,lab="rel. infectiousness = 100%",ms=10)
scatter!(size=(700,400))
xlabel!("Age group (years)")
savefig("plotting/symptomatic_rates.pdf")

"""
Load age-specific susceptibilities (for use in that scenario) and age-specific relative detection rates
for different relative infectiousness values τ = 0,0.1,0.25,0.5,1
"""
@load "data/susceptibility_rates.jld2" σ
@load "data/detection_rates_for_different_taus.jld2" d_0 d_01 d_025 d_05 d_1
rel_detection_rates = hcat(d_0,d_01,d_025,d_05,d_1)


"""
Load population sizes for each risk region
"""

popsize_risk_region = readtable("data/2019_census_age_pyramids_riskregions.csv")
N_region_age = zeros(Int64,20,17)
for i = 1:20,j = 1:17
    N_region_age[i,j] = round(Int64,popsize_risk_region[i,j+2])
end
heatmap(N_region_age)


"""
Load age mixing matrix for Kenya and convert into JLD2 array, do the same for China for baseline comparisons
"""
prop_75_79_amongst_75_plus_china = 0.7927
prop_75_79_amongst_75_plus_Kenya = sum(N_region_age[:,16])/(sum(N_region_age[:,16:17]))

function extend_and_convert_Prem_matrix(M,proportion)
    M_from_to = zeros(17,17)
    for i = 1:16,j=1:15
        M_from_to[i,j] = M[i,j]
    end
    for i = 1:16
        M_from_to[i,16] = M[i,16]*proportion
        M_from_to[i,17] = M[i,16]*(1-proportion)
    end
    for j = 1:17
        M_from_to[17,j] = M_from_to[16,j]
    end
    return Matrix(M_from_to')
end

# agemixingmatrix_table_china = readtable("data/china_baseline_age_mixing.csv")
# agemixingmatrix_china = zeros(16,16)
# for i = 1:16,j=1:16
#     agemixingmatrix_china[i,j] = agemixingmatrix_table_china[i,j]
# end
# heatmap(agemixingmatrix_china,clims = (0.,2.))
# M_China = extend_and_convert_Prem_matrix(agemixingmatrix_china,prop_75_79_amongst_75_plus_china)
# heatmap(M_China,clims = (0.,2.5))
@load "data/agemixingmatrix_china.jld2" M_China

agemixingmatrix_table = readtable("data/Kenya_age_matrix.csv")
agemixingmatrix = zeros(16,16)
for i = 1:16,j=1:16
    agemixingmatrix[i,j] = agemixingmatrix_table[i,j]
end
M_Kenya = extend_and_convert_Prem_matrix(agemixingmatrix,prop_75_79_amongst_75_plus_Kenya)
heatmap(M_Kenya,clims = (0.,2.))

@save "data/agemixingmatrix_Kenya_norestrictions.jld2" M_Kenya


agemixingmatrix_table = readtable("data/Kenya_age_matrix_home_only.csv")
agemixingmatrix = zeros(16,16)
for i = 1:16,j=1:16
    agemixingmatrix[i,j] = agemixingmatrix_table[i,j]
end
M_Kenya_ho = extend_and_convert_Prem_matrix(agemixingmatrix,prop_75_79_amongst_75_plus_Kenya)
heatmap(M_Kenya_ho,clims = (0.,2.))
@save "data/agemixingmatrix_Kenya_homeonly.jld2" M_Kenya_ho

agemixingmatrix_table = readtable("data/Kenya_age_matrix_other.csv")
agemixingmatrix = zeros(16,16)
for i = 1:16,j=1:16
    agemixingmatrix[i,j] = agemixingmatrix_table[i,j]
end
M_Kenya_other = extend_and_convert_Prem_matrix(agemixingmatrix,prop_75_79_amongst_75_plus_Kenya)

agemixingmatrix_table = readtable("data/Kenya_age_matrix_school_only.csv")
agemixingmatrix = zeros(16,16)
for i = 1:16,j=1:16
    agemixingmatrix[i,j] = agemixingmatrix_table[i,j]
end
M_Kenya_school = extend_and_convert_Prem_matrix(agemixingmatrix,prop_75_79_amongst_75_plus_Kenya)

agemixingmatrix_table = readtable("data/Kenya_age_matrix_work_only.csv")
agemixingmatrix = zeros(16,16)
for i = 1:16,j=1:16
    agemixingmatrix[i,j] = agemixingmatrix_table[i,j]
end
M_Kenya_work = extend_and_convert_Prem_matrix(agemixingmatrix,prop_75_79_amongst_75_plus_Kenya)

@save "data/agemixingmatrix_Kenya_all_types.jld2" M_Kenya M_Kenya_ho M_Kenya_other M_Kenya_school M_Kenya_work
"""
Some heatmap plots
"""
@load "data/agemixingmatrix_Kenya_all_types.jld2" M_Kenya M_Kenya_ho M_Kenya_other M_Kenya_school M_Kenya_work
heatmap(M_Kenya,clims = (0.,1))
plt_w = bar(M_Kenya_work[:,17],xticks = (1:17,age_cats),size = (800,300),title = "work contacts",ylims = (0,0.6))
plt_o = bar(M_Kenya_other[:,17],xticks = (1:17,age_cats),size = (800,300),title = "other contacts",ylims = (0,0.6))
plt_h = bar(M_Kenya_ho[:,17],xticks = (1:17,age_cats),size = (800,300),title = "home contacts",ylims = (0,0.6))
plt_s = bar(M_Kenya_school[:,17],xticks = (1:17,age_cats),size = (800,300),title = "school contacts",ylims = (0,0.6))

plt_before = heatmap(M_Kenya,xticks = (1:17,age_cats),yticks = (1:17,age_cats),size = (800,600),title = "Estimated contacts before social distancing",clims = (0,3),
                    xlabel = "Contact from individual in age group",ylabel = "Contact received by someone in age group")
plt_after = heatmap(1.5*M_Kenya_ho .+ 0.25*M_Kenya_other .+ 0.75*M_Kenya_work,xticks = (1:17,age_cats),yticks = (1:17,age_cats),size = (800,600),title = "Estimated contacts after social distancing",clims = (0,3),
            xlabel = "Contact from individual in age group",ylabel = "Contact received by someone in age group")

plt_before_and_after = groupedbar(hcat(M_Kenya[:,17], 1.5*M_Kenya_ho[:,17] .+ 0.25*M_Kenya_other[:,17] .+ 0.75*M_Kenya_work[:,17]),
                                    xticks = (1:17,age_cats),
                                    size = (800,500),
                                    title = "Estimated contact rate for a 75+ year old before and after social distancing",
                                    lab = ["Before social distancing" "After social distancing"],
                                    xlabel = "Age group 75+ year old is contacting",
                                    ylabel = "Daily number of contacts")
savefig(plt_before,"plotting/social_contacts_before_SD.pdf")
savefig(plt_after,"plotting/social_contacts_after_SD.pdf")
savefig(plt_before_and_after,"plotting/contacts_made_by_oldpeople_before_and_after_SD.pdf")
layout = @layout [a b]
plot(plt_before,plt_after,layout=layout)
plot!(size = (1200,800))


"""
Load movement matrix - then derive the P,ρ and T values
"""
@load "data/mv_matrix.jld2" mv_matrix
movements_per_person = zeros(20,20)
for i = 1:20,j=1:20
    movements_per_person[i,j] = mv_matrix[j,i]/1000#change to from col to row format for ease of multiplication with row vectors
end
heatmap(movements_per_person)
ρ = [sum(movements_per_person,dims = 1)[j]*5/30 for j = 1:20]
P_dest = zeros(20,20)#probability distribution of destination
for j = 1:20
    P_dest[:,j] = normalize(movements_per_person[:,j],1)
end
T = zeros(20,20) # combined location density matrix
for i = 1:20,j=1:20
    if i == j
        T[i,j] = 1-ρ[i]
    else
        T[i,j] = ρ[j]*P_dest[i,j]
    end
end
heatmap(T,clims = (0.,0.1))


"""
Save all the data
"""
@save "data/data_for_age_structuredmodel.jld2" N_region_age M_Kenya movements_per_person P_dest ρ T σ rel_detection_rates M_Kenya_ho hosp_rate_by_age
