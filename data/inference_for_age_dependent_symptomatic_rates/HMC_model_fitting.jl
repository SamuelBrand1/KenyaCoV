# File used for model fitting to the Kenyan data

using TransformVariables, LogDensityProblems, DynamicHMC, DynamicHMC.Diagnostics, Distributions,MAT,Plots,StatsPlots
using MCMCDiagnostics, LinearAlgebra
using Parameters, Statistics, Random
import ForwardDiff
using CSV,DataFrames,JLD2

#Age-structured contact matrix from Prem et al
@load "data/agemixingmatrix_Kenya_all_types.jld2" M_Kenya # general contacts in kenya (home, school, work, other, This is already extended into 17 age groups

#Kenyan case data for confirmed cases broken down by age and symptomatic vs. not symptomatic
Kenya_Case_Dis = CSV.read("/Users/Ojal/Documents/Covid-19/R_models/Kenya_Case_Data.csv")

#Kenyan population pyramid
Kenya_PP_County =  CSV.read("data/2019_census_age_pyramids_counties.csv")
Kenya_PP = sum(Matrix(Kenya_PP_County[2:18]), dims=1)
Kenya_PP = convert(Vector{Float64},vec(Kenya_PP))
LinearAlgebra.normalize!(Kenya_PP,1)

age_cats = ["0-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39",
            "40-44","45-49","50-54","55-59","60-64","65-69","70-74","75-79","80+"];

plt_heatmap = heatmap(M_Kenya,xticks=(1:2:17,age_cats[1:2:17]),yticks=(1:2:17,age_cats[1:2:17]),
    title = "Age-structured contact rates: Kenya (all contacts)",
    xlabel = "Contact from individual in age group",
    ylabel = "Contact received by someone in age group ")
savefig(plt_heatmap,"plotting/Kenya_all_contacts.pdf")

#Calculate eigenvalues and eigenvectors of the contact matrix. The (L_1) normalised leading eigenvector
#is the case distribution **if**
evals,evects = eigen(M_Kenya);

R₀= Real(evals[end])
v = Real.(evects[:,end])
Case_Dist = LinearAlgebra.normalize(convert(Vector{Float64},vec(sum(Matrix(Kenya_Case_Dis[:]),dims=2))),1) # Sums symtomatic and asymptomatic case by age then normalises

plt = groupedbar(hcat(Kenya_PP,v,Case_Dist),lab = ["Kenya pop." "Age indep. pred." "Conf. cases"],
        xticks=(1:17,age_cats),
        title = "Age profile of confirmed Kenyan COVID-19 cases (26th April)",
        ylabel = "Proportion of cases",
        xlabel = "Age of case")
savefig(plt,"plotting/age_profile_cases.pdf") #Uniform attack rate does not seem to expalin the case distribution. Attack rate based solely on social contacts overstimates the contribution of young age groups to cases
