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

# Model fitting

Kenya_Case_Dis_MAT = Matrix(Kenya_Case_Dis[:])

#This is a struct for the case data in numbers rather than proportions
struct CaseDistribution
    C::Array{Int64,2} #Total number of cases in each age group by symptomatic category
    n::Int64 #Total number of cases
    prediction::Function #This is the function that makes a prediction about case distribution
end

#This is the log-likelihood function (a callable function on CaseDistribution structs).
# Use a binomial likelihood at the age group level

function (cases::CaseDistribution)(θ) #Susceptibility form of likelihood
    @unpack θ = θ               # extract the parameters θ is a generic label for any set of params, the difference is in the prediction function
    @unpack n, C, prediction = cases   # extract the data and prediction function
    T = eltype(θ)
    p = prediction(θ) #predicted number of cases, asymptomatic and symptomatic
    obs_tot_cases_age = vec(sum(C,dims=2))  # Observed total cases by age
    pred_tot_cases_age = vec(sum(p,dims=2))  # predicted total cases by age
    pred_agedist = vec(pred_tot_cases_age./sum(p)) # predicted case (both symp and asymp) distribution by age
    pred_symp = p[:,2]./pred_tot_cases_age # predicted symptomatic rate by age

    try
        logL = 0.
        for i = 1:17
            logL += C[i,1]*log(1-pred_symp[i]) +  C[i,2]*log(pred_symp[i]) + obs_tot_cases_age[i]*log(pred_agedist[i])  # C[,1] has number asymptomatic and C[,2] has number symptomatic. The third element of the log lik is a multinomial distribution of cases by age.
        end
        return logL
    catch errtype
        return T(-Inf)
    end

end

#Function for calculating the case distribution
#Prediction  function assumes mean 2 days of pre-symptomatic transmission, mean 7 days infectious period for
#mild infecteds and that the eventually severe cases don't transmit more than mild cases

function pred_case_distribution_using_iter_K_model2_splitAsymp(χ::Vector,d::Vector,ϵ::Vector,C)
    d1, = size(C) # number of age groups
    K = zeros(eltype(d),d1,d1)
    for a = 1:d1,b=1:d1
        K[a,b] = χ[a]*C[a,b]*(2*ϵ[b] + 7*d[b] + 7*(1-d[b])*ϵ[b])
    end
    v = (K^10)*ones(d1)
    asymp = d.*v  # predicted asymptomatics cases
    syp = (1-d).*v  # predicted symptomatics cases
    symp = v - asymp # predicted symptomatics cases

    return hcat(asymp,symp)
end

# NEED to add functions here to call model run and get expected case distibution given interventions current fit would assume interventions have had full effect and we are in a new steady state

#This function returns the HMC chains with tree statistics for the run
function HMC_for_detection_rate(χ,ϵ,n_draws)
    cases_to_fit = CaseDistribution(Array{Int64,2}(Kenya_Case_Dis_MAT),
                        Int64(sum(Kenya_Case_Dis_MAT)),
                        d -> pred_case_distribution_using_iter_K_model2_splitAsymp(χ,d,ϵ*ones(17),M_Kenya))

    trans = as((θ = as(Array, asℝ₊, 17),)) # All parameters are transformed to be positive.
    P = TransformedLogDensity(trans, cases_to_fit) #This creates a transformed log-likelihood
    ∇P = ADgradient(:ForwardDiff, P) #This automatically generates a log-likelihood gradient at the same time as the likelihood is called
    return results = mcmc_with_warmup(Random.GLOBAL_RNG, ∇P, n_draws) # reporter = NoProgressReport()
end

χ_zhang = vcat(0.34*ones(3),ones(10),1.47*ones(4))  # age specific susceptibility

#d_chain_model2_epsilon_0 = HMC_for_detection_rate(χ_zhang,0.,10000)
#d_chain_model2_epsilon_01 = HMC_for_detection_rate(χ_zhang,0.1,10000)
#d_chain_model2_epsilon_025 = HMC_for_detection_rate(χ_zhang,0.25,10000)
#d_chain_model2_epsilon_05 = HMC_for_detection_rate(χ_zhang,0.5,10000)
d_chain_model2_epsilon_1 = HMC_for_detection_rate(χ_zhang,1.,10)

@save "HMC_chains_for_model2.jld2" d_chain_model2_epsilon_0 d_chain_model2_epsilon_01 d_chain_model2_epsilon_025 d_chain_model2_epsilon_05 d_chain_model2_epsilon_1


zhang_sus_dist = pred_case_distribution_using_iter_K_model2(χ_zhang,ones(17),ones(17),China_to_from)
_zhang_sus_dist = convert_to_data_agegroups(zhang_sus_dist)
plt = groupedbar(hcat(_zhang_sus_dist,_v,Case_Dist),lab = ["Zhang sus. profile only" "Age indep. pred." "Conf. cases"],
        xticks=(1:9,["0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+"]),
        title = "Age profile of confirmed Chinese COVID-19 cases (11/2)",
        ylabel = "Proportion of cases",
        xlabel = "Age of case")
