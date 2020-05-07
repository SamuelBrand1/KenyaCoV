# File used for model fitting to the Kenyan data

using TransformVariables, LogDensityProblems, DynamicHMC, DynamicHMC.Diagnostics, Distributions,MAT,Plots,StatsPlots
using MCMCDiagnostics, LinearAlgebra
using Parameters, Statistics, Random
import ForwardDiff
using CSV,DataFrames,JLD2

#Age-structured contact matrix from Prem et al
@load "data/agemixingmatrix_Kenya_all_types.jld2" M_Kenya # general contacts in kenya (home, school, work, other, This is already extended into 17 age groups

#Kenyan case data for confirmed cases broken down by age and symptomatic vs. not symptomatic
# To continue  from here
file = matopen("CCDC_Case_Dist_Feb_11.mat")
CCDC_Case_Dist_Feb_11 = read(file, "CCDC_Case_Dist_11_feb")
close(file)
#Chinese population pyramid
file = matopen("Pop_Pyramids.mat")
China_PP = read(file, "China_PP")
close(file)
#Group the population pyramid according to the CCDC report age groups
_China_PP = zeros(9)
for i = 1:8
    _China_PP[i] = China_PP[2*(i-1)+ 1] + China_PP[2*(i-1)+ 2]
end
_China_PP[9] = sum(China_PP[19:end])
LinearAlgebra.normalize!(_China_PP,1);

age_cats = ["0-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39",
            "40-44","45-49","50-54","55-59","60-64","65-69","70-74","75-79","80+"];
