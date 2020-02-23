#Script for building the population pyramid into data
using DataFrames,Plots,JLD2

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
@save "data/populationpyramid.jld" prop_by_agegroup

"""
Load age mixing matrix and convert into JLD2 array
"""
agemixingmatrix_table = readtable("data/Kenya_age_matrix.csv")
agemixingmatrix = zeros(16,16)
for i = 1:16,j=1:16
    agemixingmatrix[i,j] = agemixingmatrix_table[i,j]
end
heatmap(agemixingmatrix)
@save "data/agemixingmatrix.jld2" agemixingmatrix
