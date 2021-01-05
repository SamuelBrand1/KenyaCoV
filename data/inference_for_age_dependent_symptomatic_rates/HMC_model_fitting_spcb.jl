# File used for model fitting to the Kenyan data

using TransformVariables, LogDensityProblems, DynamicHMC, DynamicHMC.Diagnostics, Distributions,MAT,Plots,StatsPlots
using MCMCDiagnostics, LinearAlgebra
using Parameters, Statistics, Random
import ForwardDiff
using CSV,DataFrames,JLD2
using BenchmarkTools

#Age-structured contact matrix from Prem et al
# @load "data/agemixingmatrix_Kenya_all_types.jld2" M_Kenya # general contacts in kenya (home, school, work, other, This is already extended into 17 age groups
@load "data/agemixingmatrix_Kenya_all_types.jld2"
@load "data/data_hh_model.jld2" hh_mat_vector #hh_mat_vector is a 47-vector of the household joint occupancy

"""
Method that calculates the density dependent transmission rate in HHs so that the
secondary attack rate from a symptomatic case within the household is SAR... this seems to be estimated at about 20% SAR
"""
function transmissionrateinHH(ϵ,γ,σ₂,SAR)
    @unpack ϵ,γ,σ₂ = p
    return (-(ϵ*γ + σ₂) + sqrt((ϵ*γ + σ₂)^2 + 4*ϵ*σ₂*γ*(SAR/(1-SAR))))/(2*ϵ)
end


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
        size = (800,400),
        title = "Age profile of confirmed Kenyan COVID-19 cases (19th June 2020)",
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
    logprior::Function #Using log-prior direct because Distributions package has efficient implementations of the logpdf for each of its distributions
end


#This is the log-likelihood function (a callable function on CaseDistribution structs).
# Use a binomial likelihood at the age group level

function (cases::CaseDistribution)(θ) #Susceptibility form of likelihood
    @unpack χ,d,ϵ = θ               # extract the parameters θ is a generic label for any set of params, the difference is in the prediction function
    @unpack C, n, prediction,logprior = cases   # extract the data and prediction function
    T = eltype(χ)
    p = prediction(χ,d,ϵ) #predicted number of cases, asymptomatic and symptomatic
    obs_tot_cases_age = vec(sum(C,dims=2))  # Observed total cases by age
    pred_tot_cases_age = vec(sum(p,dims=2))  # predicted total cases by age
    pred_agedist = vec(pred_tot_cases_age./sum(p)) # predicted case (both symp and asymp) distribution by age
    pred_symp = p[:,2]./pred_tot_cases_age # predicted symptomatic rate by age

    try
        logL = 0.
        for i = 1:17
            logL += C[i,1]*log(1-pred_symp[i]) +  C[i,2]*log(pred_symp[i]) + obs_tot_cases_age[i]*log(pred_agedist[i]) + logprior(χ,d,ϵ)  # C[,1] has number asymptomatic and C[,2] has number symptomatic. The third element of the log lik is a multinomial distribution of cases by age.
        end
        return logL
    catch errtype
        return T(-Inf)
    end

end

#Function for calculating the case distribution
#Prediction  function assumes mean 2 days of pre-symptomatic transmission, mean 7 days infectious period for
#mild infecteds and that the eventually severe cases don't transmit more than mild cases

function logprior_zhang(χ,d,ϵ)
    ll = 0.
    for a = 1:3
        ll += logpdf(Gamma(1,1),χ[a])
        #ll += logpdf(Gamma(9,0.0378),χ[a])
    end
    for a = 4:13
        #ll += logpdf(Gamma(9,0.1111),χ[a])
        ll += logpdf(Gamma(1,1),χ[a])
    end
    for a = 14:17
        #ll += logpdf(Gamma(9,0.163333),χ[a])
        ll += logpdf(Gamma(1,1),χ[a])
    end
    return ll
end

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

#Proposing a slightly different model which constrains the attack rate for symptomatic cases in the household to fit an estimate (usually ~20%)
function pred_case_distribution_HH(χ::Vector,d::Vector,ϵ::Vector,C_out,C_in,SAR)
    d1, = size(C) # number of age groups
    K = zeros(eltype(d),d1,d1)
    for a = 1:d1,b=1:d1
        K[a,b] = χ[a]*C_out[a,b]*(2*ϵ[1] + 7*d[b] + 7*(1-d[b])*ϵ[1]) + χ[a]*SAR*C_in[a,b]
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
        (χ,d,ϵ) -> pred_case_distribution_using_iter_K_model2_splitAsymp(χ,d,ϵ,0.5*M_Kenya_work .+ 0.5*M_Kenya_other .+ M_Kenya_ho),
        logprior_zhang)
    trans = as((χ = as(Array, asℝ₊, 17),d= as(Array, asℝ₊, 17), ϵ=as(Array, asℝ₊, 1))) # All parameters are transformed to be positive.
    P = TransformedLogDensity(trans, cases_to_fit) #This creates a transformed log-likelihood
    ∇P = ADgradient(:ForwardDiff, P) #This automatically generates a log-likelihood gradient at the same time as the likelihood is called
    #q₀ = fill(0.5,17)  # initial values
    results = mcmc_with_warmup(Random.GLOBAL_RNG, ∇P, n_draws,
                                        initialization=(q = vcat(zeros(17),log.(0.1*ones(17)),log(1.)),))#,reporter = NoProgressReport()) # reporter = NoProgressReport()
    # return TransformVariables.transform.(trans,results.chain) # Transform back to normal numbers
    return results
end


@time results = HMC_for_detection_rate(10000)
trans = as((χ = as(Array, asℝ₊, 17),d= as(Array, asℝ₊, 17), ϵ=as(Array, asℝ₊, 1)))
results_trans = TransformVariables.transform.(trans,results.chain)
res = results_trans[1]
pred_distribs = [sum(pred_case_distribution_using_iter_K_model2_splitAsymp(res.χ,res.d,res.ϵ,0.5*M_Kenya_work .+ 0.5*M_Kenya_other .+ M_Kenya_ho),dims = 2)/sum(pred_case_distribution_using_iter_K_model2_splitAsymp(res.χ,res.d,res.ϵ,0.5*M_Kenya_work .+ 0.5*M_Kenya_other .+ M_Kenya_ho)) for res in results_trans]

data_distribution = sum(Kenya_Case_Dis_MAT,dims = 2)/sum(Kenya_Case_Dis_MAT)
groupedbar(hcat(data_distribution,mean(pred_distribs)))
bar(mean([res.d for res in results_trans]))


posterior_ϵ = [res.ϵ[1] for res in results_trans]
plot(posterior_ϵ,lab = "epsilon")
histogram(posterior_ϵ)


@save "plotting/MCMC_results_transformed_with_priors.jld2" results
#@load "plotting/MCMC_results_transformed_with_priors.jld2" results

# Diagnostics
  # Tree statistics
summarize_tree_statistics(results.tree_statistics)

  # Effective sample sizes
posterior = first.(transform.(trans, results.chain)) #Apply transformation to chain, use the first paramter χ
eff_sample_sizes = Float64[]
for a = 1:17
    posterior_samples_for_group_a = [posterior[i][a] for i = 1:length(posterior)]
    push!(eff_sample_sizes,effective_sample_size(posterior_samples_for_group_a)   )
end

@save  "data/inference_for_age_dependent_symptomatic_rates/effective_sample_sizes.jld2" eff_sample_sizes
plt=histogram(eff_sample_sizes,bins = 100,lab = "",xlabel = "Eff. sample size",
            title = "Frequency of eff. sample size by age group",yticks = 0:1:15)

savefig(plt,"plotting/effective_sample_sizes.png")


# Posterior distribution

posterior = TransformVariables.transform.(trans,results.chain)

posterior_χ =[c.χ for c in posterior]
posterior_d =[c.d for c in posterior]
posterior_ϵ =[c.ϵ for c in posterior]

posterior_array_ϵ = zeros(length(posterior),1)
for i = 1:length(posterior_ϵ)
    posterior_array_ϵ[i] = posterior_ϵ[i][1]
end

plt=plot(posterior_array_ϵ, xlabel="Iteration",ylabel="Epsilon", legend=nothing)
histogram(posterior_array_ϵ)
#savefig(plt,"plotting/chain_epsilon.png")

median(posterior_array_ϵ)
quantile(posterior_array_ϵ[:],0.025)
quantile(posterior_array_ϵ[:],0.975)

posterior_array_χ = zeros(length(posterior),17)
for i = 1:length(posterior),j=1:17
     posterior_array_χ[i,j] = posterior_χ[i][j]
end

plt=plot(posterior_array_χ, xlabel="Iteration",ylabel="Chi", size=(800,400),legend=nothing)
#savefig(plt,"plotting/chain_chi.png")

med= [mean(posterior_array_χ[:,i]) for i = 1:17]
lb = [quantile(posterior_array_χ[:,i],0.025) for i = 1:17]
ub = [quantile(posterior_array_χ[:,i],0.975) for i = 1:17]


fig = scatter(med,yerr = (lb,ub),fillalpha = 0.3,legend = nothing,
        xticks = (1:17,age_cats), size=(800,400),
        ylabel = "Susceptibility",xlabel="Age groups")
#savefig(fig,"plotting/fitted_age_susceptibility.png")



posterior_array_d = zeros(length(posterior),17)
for i = 1:length(posterior),j=1:17
     posterior_array_d[i,j] = posterior_d[i][j]
end

plt=plot(posterior_array_d, xlabel="Iteration",ylabel="Symptomatic rate", size=(800,400),legend=nothing)
#savefig(plt,"plotting/chain_Symp_rate.png")

med= [mean(posterior_array_d[:,i]) for i = 1:17]
lb = [quantile(posterior_array_d[:,i],0.025) for i = 1:17]
ub = [quantile(posterior_array_d[:,i],0.975) for i = 1:17]

fig = scatter(med,yerr = (lb,ub),fillalpha = 0.3,legend = nothing,
        xticks = (1:17,age_cats), size=(800,400), ylims=(0,1),
        ylabel = "Symptomatic rate",xlabel="Age groups")
#savefig(fig,"plotting/fitted_symptomatic_rates.png")


# Plot fit to data
obs_symp = vec( Kenya_Case_Dis_MAT[:,2]./sum(Kenya_Case_Dis_MAT,dims=2))
obs_agedist = Case_Dist

posterior_symp_rates = zeros(17,length(posterior))
posterior_case_dist  = zeros(17,length(posterior))
for i in 1:length(posterior)
    res = pred_case_distribution_using_iter_K_model2_splitAsymp(posterior_array_χ[i,:],posterior_array_d[i,:],posterior_array_ϵ[i,:],0.5*M_Kenya_work .+ 0.5*M_Kenya_other .+ M_Kenya_ho)
    pred_tot_cases_age = vec(sum(res,dims=2))  # predicted total cases by age
    pred_symp = res[:,2]./pred_tot_cases_age # predicted symptomatic rate by age
    pred_agedist = vec(pred_tot_cases_age./sum(res)) # predicted case (both symp and asymp) distribution by age
    posterior_symp_rates[:,i] = pred_symp
    posterior_case_dist[:,i] = pred_agedist

end

posterior_symp_rates = median(posterior_symp_rates,dims=2)
posterior_case_dist = median(posterior_case_dist,dims=2)

plt1= groupedbar(hcat(Case_Dist,posterior_case_dist),lab = ["Obs. case dist." "Pred. case dist."],xticks=(1:17,age_cats),size=(800,400))
#savefig(plt1,"plotting/obs_and_fitted_case_distribution.png")
plt2= groupedbar(hcat(obs_symp,posterior_symp_rates),lab = ["True symp rate." "Pred. symp rate"],xticks=(1:17,age_cats), size=(800,400))
#savefig(plt2,"plotting/obs_and_fitted_symptomatic_rates.png")
