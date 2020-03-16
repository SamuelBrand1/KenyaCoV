#Script for building the population pyramid into data
using DataFrames,Plots,JLD2,MAT,LinearAlgebra,Distances,StatsPlots
using LinearAlgebra:normalize,normalize!


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

agemixingmatrix_table_china = readtable("data/china_baseline_age_mixing.csv")
agemixingmatrix_china = zeros(16,16)
for i = 1:16,j=1:16
    agemixingmatrix_china[i,j] = agemixingmatrix_table_china[i,j]
end
heatmap(agemixingmatrix_china,clims = (0.,2.))
M_China = extend_and_convert_Prem_matrix(agemixingmatrix_china,prop_75_79_amongst_75_plus_china)
heatmap(M_China,clims = (0.,2.5))
@save "data/agemixingmatrix_china.jld2" M_China

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

reduced_treatment_rates = [(0.,1),(1/3.5,0.5)]

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
Save all the data
"""
@save "data/data_for_age_structuredmodel.jld2" N_region_age M_Kenya movements_per_person P_dest ρ T σ rel_detection_rates M_Kenya_ho
