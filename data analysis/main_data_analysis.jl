#Script for building the population pyramid into data
using DataFrames,Plots,JLD2,MAT,LinearAlgebra,Distances,StatsPlots,CSV,DelimitedFiles,RData
using LinearAlgebra:normalize,normalize!
# using TransformVariables, DynamicHMC, DynamicHMC.Diagnostics, Distributions

JLD2.@load("data/agemixingmatrix_Kenya_homeonly.jld2")

heatmap(M_Kenya_ho)
RData.load("data/hh_mixing_matrices_80.rda")





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


agemixingmatrix_table = DataFrame!(CSV.File("data/Kenya_age_matrix_home_only.csv"))
agemixingmatrix = zeros(16,16)
for i = 1:16,j=1:16
    agemixingmatrix[i,j] = agemixingmatrix_table[i,j]
end
heatmap(agemixingmatrix)
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


# Some heatmap plots of age mixing matrices

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



## Load movement matrix - then derive the P,ρ and T values

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
@save "data/data_for_age_structuredmodel.jld2" N_region_age M_Kenya movements_per_person P_dest ρ T σ rel_detection_rates M_Kenya_ho hosp_rate_by_age ICU_rate_by_age_cond_hosp


"""
After using analysis of movements --- data for county model
"""
N_region_age = N_pop
T = T_opt
ρ = ρ_county
P_dest = P_opt
@load "data/agemixingmatrix_Kenya_norestrictions.jld2" M_Kenya
@load "data/agemixingmatrix_Kenya_homeonly.jld2" M_Kenya_ho

σ =  vcat(0.34*ones(3),ones(10),1.47*ones(4))
@load "data/detection_rates_for_different_epsilons_model2.jld2" d_0 d_01 d_025 d_05 d_1
rel_detection_rates = hcat(d_0,d_01,d_025,d_05,d_1)

@save "data/data_for_age_structuredmodel_with_counties.jld2" N_region_age M_Kenya movements_per_person P_dest ρ T σ rel_detection_rates M_Kenya_ho hosp_rate_by_age ICU_rate_by_age_cond_hosp
