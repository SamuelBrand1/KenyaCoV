#Script for building the population pyramid into data
using DataFrames,Plots,JLD2,MAT,LinearAlgebra,Distances,StatsPlots
using LinearAlgebra:normalize,normalize!

"""
Age dependent reporting for MERS-like scenario

Age    Susceptibility   Reporting
0-4      0.00282027   0.0142415
5-9      0.00364015   0.0196252
10-14  0.00288739   0.0114435
15-19  0.0120954   0.0273375
20-24  0.0377834   0.122321
25-29  0.0579848   0.303001
30-34  0.0685311   0.42322
35-39  0.067254   0.390164
40-44  0.0615844   0.317926
45-49  0.0965786   0.465243
50-54  0.118037   0.480286
55-59  0.15693   0.459481
60-64  0.284027   0.727757
65-69  0.393072   0.874986
70-74  0.411986   0.555446
75-79  0.703186   0.320811
80-84  1.32085   0.185529
85-89  1.08871   0.0708618
90-94  0.577218   0.0115956
95-99  0.830651   0.00386519
"""
reporting_rate = [0.0142415, 0.0196252,0.0114435,0.0273375,0.122321,
                    0.303001,0.42322,0.390164,0.317926,0.465243,
                    0.480286,0.459481,0.727757,0.874986,0.555446,0.320811]

@save "data/reporting_rate_for_MERS_like_scenario.jld2" reporting_rate


"""
World Pop data: First 47 rows are counties, then next 20 rows are risk regions
"""
# poppyramids_tbl = readtable("data/agePyramid_with_riskregions.csv")
# # poppyramids_tbl = poppyramids_tbl[1:47,6:22]
# poppyramids_counties = zeros(47,16);
# for i = 1:47,a = 1:16
#     if a < 16
#         poppyramids_counties[i,a] = poppyramids_tbl[i,a+1]
#     else
#         poppyramids_counties[i,16] = poppyramids_tbl[i,16+1] + poppyramids_tbl[i,17+1]
#     end
# end
# for i = 1:47
#     poppyramids_counties[i,:] .= normalize(poppyramids_counties[i,:],1)
# end
# poppyramids_riskregions = zeros(20,16);
# for i = 1:20,a = 1:16
#     if a < 16
#         poppyramids_riskregions[i,a] = poppyramids_tbl[i+47,a+1]
#     else
#         poppyramids_riskregions[i,16] = poppyramids_tbl[i+47,16+1] + poppyramids_tbl[i+47,17+1]
#     end
# end
# for i = 1:20
#     poppyramids_riskregions[i,:] .= normalize(poppyramids_riskregions[i,:],1)
# end
# @save "data/populationpyramids_by_riskregion.jld2" poppyramids_riskregions
"""
Worldpop pyramid for the whole of Kenya --- this is now
"""

# poppyramid = readtable("data/Kenya_population_pyramid.csv")
#
# N = 0.
# d1, = size(poppyramid)
# for i = 1:d1
#     global N
#     N += poppyramid[i,2] + poppyramid[i,3]
# end
# prop_by_agegroup = [(poppyramid[i,2] + poppyramid[i,3])/N for i = 1:15]
# push!(prop_by_agegroup,1 - sum(prop_by_agegroup)) #Older people groupped in to last category
# bar(prop_by_agegroup)
# @save "data/populationpyramid.jld2" prop_by_agegroup

"""
Load age mixing matrix for Kenya and convert into JLD2 array, do the same for China for baseline comparisons
"""
agemixingmatrix_table_china = readtable("data/china_baseline_age_mixing.csv")
agemixingmatrix_china = zeros(16,16)
for i = 1:16,j=1:16
    agemixingmatrix_china[i,j] = agemixingmatrix_table_china[i,j]
end
heatmap(agemixingmatrix_china,clims = (0.,5.))
M_China = agemixingmatrix_china
@save "data/agemixingmatrix_china.jld2" M_China

agemixingmatrix_table = readtable("data/Kenya_age_matrix.csv")
agemixingmatrix = zeros(16,16)
for i = 1:16,j=1:16
    agemixingmatrix[i,j] = agemixingmatrix_table[i,j]
end


@save "data/agemixingmatrix.jld2" agemixingmatrix
"""
Load population sizes for each risk region
"""

popsize_risk_region = readtable("data/2019_census_age_pyramids_riskregions.csv")

N_region_age = zeros(Int64,20,16)
for i = 1:20,j = 1:16
    if j < 16
        N_region_age[i,j] = round(Int64,popsize_risk_region[i,j+2])
    else
        N_region_age[i,16] = round(Int64,popsize_risk_region[i,18] + popsize_risk_region[i,19])
    end
end
heatmap(N_region_age)
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
    P_dest[:,j] = LinearAlgebra.normalize(movements_per_person[:,j],1)
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
Load age-specific susceptibility
"""
age_specific_tbl = readtable("data/sus_profile_yang.csv")
age_specific_sus = [age_specific_tbl[i,1] for i = 1:16]
"""
Save all the data
"""
@save "data/data_for_age_structuredmodel.jld2" N_region_age agemixingmatrix movements_per_person P_dest ρ T age_specific_sus
