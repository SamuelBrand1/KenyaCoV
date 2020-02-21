#Script for building the population pyramid into data
x = 1
using DataFrames,Plots
poppyramid = readtable("data/Kenya_population_pyramid.csv")
"""
Convert into percentages for each of 16 age groups used in
the age structured models
"""
N = 0.
d1, = size(poppyramid)
for i = 1:d1
    global N
    N += poppyramid[i,2] + poppyramid[i,3]
end
prop_by_agegroup = [(poppyramid[i,2] + poppyramid[i,3])/N for i = 1:15]
push!(prop_by_agegroup,1 - sum(prop_by_agegroup)) #Older people groupped in to last category
bar(prop_by_agegroup)
