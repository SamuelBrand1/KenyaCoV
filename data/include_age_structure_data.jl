#Script for building the population pyramid into data
using DataFrames,Plots,JLD2,MAT,LinearAlgebra

"""
Load population pyramid and convert into percentages for each of 16 age groups used in
the age structured models
"""
poppyramid = readtable("data/Kenya_population_pyramid.csv")

N = 0.
d1, = size(poppyramid)
for i = 1:d1
    global N
    N += poppyramid[i,2] + poppyramid[i,3]
end
prop_by_agegroup = [(poppyramid[i,2] + poppyramid[i,3])/N for i = 1:15]
push!(prop_by_agegroup,1 - sum(prop_by_agegroup)) #Older people groupped in to last category
bar(prop_by_agegroup)
@save "data/populationpyramid.jld2" prop_by_agegroup

"""
Load age mixing matrix and convert into JLD2 array
"""
agemixingmatrix_table = readtable("data/Kenya_age_matrix.csv")
agemixingmatrix = zeros(16,16)
for i = 1:16,j=1:16
    agemixingmatrix[i,j] = agemixingmatrix_table[i,j]
end

heatmap(agemixingmatrix,clims = (0.,5.))
@save "data/agemixingmatrix.jld2" agemixingmatrix
"""
Load population sizes for each risk region
"""
popsize_risk_region = readtable("data/population_risk_regions_2019.csv")
N_region = [round(Int64,popsize_risk_region[i,2]) for i = 1:20 ]
N_region_age = zeros(Int64,20,16)
for i = 1:20,j = 1:16
    N_region_age[i,j] = round(Int64,popsize_risk_region[i,2]*prop_by_agegroup[j])
end
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
P = zeros(20,20)#probability distribution of destination
for j = 1:20
    P[:,j] = LinearAlgebra.normalize(movements_per_person[:,j],1)
end
T = zeros(20,20) # combined location density matrix
for i = 1:20,j=1:20
    if i == j
        T[i,j] = 1-ρ[i]
    else
        T[i,j] = ρ[j]*P[i,j]
    end
end
heatmap(T,clims = (0.,0.1))
"""
Save all the data
"""
@save "data/data_for_age_structuredmodel.jld2" N_region_age agemixingmatrix movements_per_person P ρ T
