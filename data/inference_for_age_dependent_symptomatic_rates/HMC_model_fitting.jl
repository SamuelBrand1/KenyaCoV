# File used for model fitting to the Kenyan data

using TransformVariables, LogDensityProblems, DynamicHMC, DynamicHMC.Diagnostics, Distributions,MAT,Plots,StatsPlots
using MCMCDiagnostics, LinearAlgebra
using Parameters, Statistics, Random
import ForwardDiff
using CSV,DataFrames,JLD2
using BenchmarkTools

#Age-structured contact matrix from Prem et al
@load "data/agemixingmatrix_Kenya_all_types.jld2" M_Kenya # general contacts in kenya (home, school, work, other, This is already extended into 17 age groups

#Kenyan case data for confirmed cases broken down by age and symptomatic vs. not symptomatic
Kenya_Case_Dis = CSV.read("data/inference_for_age_dependent_symptomatic_rates/Kenya_Case_Data.csv")

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
# savefig(plt_heatmap,"plotting/Kenya_all_contacts.pdf")

#Calculate eigenvalues and eigenvectors of the contact matrix. The (L_1) normalised leading eigenvector
#is the case distribution **if**
evals,evects = eigen(M_Kenya);

R₀= Real(evals[end])
v = Real.(evects[:,end])
Case_Dist = LinearAlgebra.normalize(convert(Vector{Float64},vec(sum(Matrix(Kenya_Case_Dis[:]),dims=2))),1) # Sums symtomatic and asymptomatic case by age then normalises

plt = groupedbar(hcat(Kenya_PP,v,Case_Dist),lab = ["Kenya pop." "Age indep. pred." "Conf. cases"],
        xticks=(1:17,age_cats),
        title = "Age profile of confirmed Kenyan COVID-19 cases (21 May 2020)",
        ylabel = "Proportion of cases",
        xlabel = "Age of case")
# savefig(plt,"plotting/age_profile_cases.png") #Uniform attack rate does not seem to expalin the case distribution. Attack rate based solely on social contacts overstimates the contribution of young age groups to cases

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
    @unpack χ,d,ϵ = θ               # extract the parameters θ is a generic label for any set of params, the difference is in the prediction function
    @unpack C, n, prediction = cases   # extract the data and prediction function
    T = eltype(χ)
    p = prediction(χ,d,ϵ) #predicted number of cases, asymptomatic and symptomatic
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
        K[a,b] = χ[a]*C[a,b]*(2*ϵ[1] + 7*d[b] + 7*(1-d[b])*ϵ[1])
    end
    v = (K^10)*ones(d1)
    symp = d.*v  # predicted symptomatics cases
    asymp = v - symp # predicted asymptomatics cases
    return hcat(asymp,symp)
end

# NEED to add functions here to call model run and get expected case distibution given interventions current fit would assume interventions have had full effect and we are in a new steady state

#This function returns the HMC chains with tree statistics for the run
function HMC_for_detection_rate(n_draws)
    cases_to_fit = CaseDistribution(Array{Int64,2}(Kenya_Case_Dis_MAT),
        Int64(sum(Kenya_Case_Dis_MAT)),
        (χ,d,ϵ) -> pred_case_distribution_using_iter_K_model2_splitAsymp(χ,d,ϵ,M_Kenya))
    trans = as((χ = as(Array, asℝ₊, 17),d= as(Array, asℝ₊, 17), ϵ=as(Array, asℝ₊, 1))) # All parameters are transformed to be positive.
    P = TransformedLogDensity(trans, cases_to_fit) #This creates a transformed log-likelihood
    ∇P = ADgradient(:ForwardDiff, P) #This automatically generates a log-likelihood gradient at the same time as the likelihood is called
    #q₀ = fill(0.5,17)  # initial values
    results = mcmc_with_warmup(Random.GLOBAL_RNG, ∇P, n_draws,
                                        initialization=(q = vcat(zeros(17),log.(0.1*ones(17)),log(1.)),) ) # reporter = NoProgressReport()
    # return TransformVariables.transform.(trans,results.chain)
    return results
end


@time results = HMC_for_detection_rate(10)
@save "MCMC_results_.jld2" results

# Diagnostics
trans = as((χ = as(Array, asℝ₊, 17),d= as(Array, asℝ₊, 17), ϵ=as(Array, asℝ₊, 1)))
  # Tree statistics
summarize_tree_statistics(results.tree_statistics)

  # Effective sample sizes
posterior = first.(transform.(trans, results.chain)) #Apply transformation to chain, use the first paramter χ
eff_sample_sizes = Float64[]
for a = 1:17
    posterior_samples_for_group_a = [posterior[i][a] for i = 1:length(posterior)]
    push!(eff_sample_sizes,effective_sample_size(posterior_samples_for_group_a)   )
end

plt=histogram(eff_sample_sizes,bins = 100,lab = "",xlabel = "Eff. sample size",
            title = "Frequency of eff. sample size by age group",
            xticks = (1:1:10),yticks = 0:1:15)

savefig(plt,"plotting/effective_sample_sizes.png")

# Posterior distribution
trans = as((χ = as(Array, asℝ₊, 17),d= as(Array, asℝ₊, 17), ϵ=as(Array, asℝ₊, 1)))
posterior = TransformVariables.transform.(trans,results.chain)

posterior_χ =[c.χ for c in posterior]
posterior_d =[c.d for c in posterior]
posterior_ϵ =[c.ϵ for c in posterior]

posterior_array_ϵ = zeros(length(posterior),1)
for i = 1:length(posterior_ϵ)
    posterior_array_ϵ[i] = posterior_ϵ[i][1]
end
median(posterior_array_ϵ)
quantile(posterior_array_ϵ[:],0.025)
quantile(posterior_array_ϵ[:],0.975)

posterior_array_χ = zeros(length(posterior),17)
for i = 1:length(posterior),j=1:17
     posterior_array_χ[i,j] = posterior_χ[i][j]
end
med= [mean(posterior_array_χ[:,i]) for i = 1:17]
lb = [quantile(posterior_array_χ[:,i],0.025) for i = 1:17]
ub = [quantile(posterior_array_χ[:,i],0.975) for i = 1:17]


fig = scatter(med,yerr = (lb,ub),fillalpha = 0.3,legend = nothing,
        xticks = (1:17,age_cats), size=(800,400),
        ylabel = "Susceptibility",xlabel="Age groups",
        yscale = :log10)
savefig(fig,"plotting/fitted_age_susceptibility.png")



posterior_array_d = zeros(length(posterior),17)
for i = 1:length(posterior),j=1:17
     posterior_array_d[i,j] = posterior_d[i][j]
end
med= [mean(posterior_array_d[:,i]) for i = 1:17]
lb = [quantile(posterior_array_d[:,i],0.025) for i = 1:17]
ub = [quantile(posterior_array_d[:,i],0.975) for i = 1:17]

fig = scatter(med,yerr = (lb,ub),fillalpha = 0.3,legend = nothing,
        xticks = (1:17,age_cats), size=(800,400), ylims=(0,1),
        ylabel = "Symptomatic rate",xlabel="Age groups")
savefig(fig,"plotting/fitted_symptomatic_rates.png")


# Plot fit to data
function posterior_prediction(result)
    posterior = transform.(trans, result.chain)
    pred = map((χ,d,ϵ) -> pred_case_distribution_using_iter_K_model2(χ,d,ϵ,M_Kenya),posterior);
    return groupedbar(hcat(Case_Dist,mean(pred)),lab = ["True distrib." "Post. prediction"])
end
plt1 = posterior_prediction(results);
