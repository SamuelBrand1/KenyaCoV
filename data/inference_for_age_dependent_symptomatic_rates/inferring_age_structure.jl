
using TransformVariables, LogDensityProblems, DynamicHMC, DynamicHMC.Diagnostics, Distributions,MAT,Plots,StatsPlots
using MCMCDiagnostics, LinearAlgebra
using Parameters, Statistics, Random
import ForwardDiff
using CSV,DataFrames,JLD2

#Home-setting age-structured contact matrix from Prem et al
file = matopen("data/inference_for_age_dependent_symptomatic_rates/China_Mixing_ho.mat")
China_Mixing_ho = read(file, "China_Mixing_ho")
close(file)
#CCDC case data for confirmed cases broken down by age
file = matopen("data/inference_for_age_dependent_symptomatic_rates/CCDC_Case_Dist_Feb_11.mat")
CCDC_Case_Dist_Feb_11 = read(file, "CCDC_Case_Dist_11_feb")
close(file)
#Chinese population pyramid
file = matopen("data/inference_for_age_dependent_symptomatic_rates/Pop_Pyramids.mat")
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

#Extend the Prem contact matrix as described above
China_from_to = zeros(17,17);
ratio = China_PP[19]/sum(China_PP[19:end])
for i = 1:16,j=1:15
    China_from_to[i,j] = China_Mixing_ho[i,j]
end
for i = 1:16
    China_from_to[i,16] = China_Mixing_ho[i,16]*ratio
    China_from_to[i,17] = China_Mixing_ho[i,16]*(1-ratio)
end
for j = 1:17
    China_from_to[17,j] = China_from_to[16,j]
end
China_to_from = Matrix(China_from_to')

#Method for converting from the 17 5-yearly age groups to the 9 10-yearly age groups in the data
function convert_to_data_agegroups(v)
    _v = zeros(9)
    for i = 1:8
        _v[i] = v[2*(i-1)+1] + v[2*(i-1)+2]
    end
    _v[9] = v[17]
    return LinearAlgebra.normalize!(_v,1)
end
plt_heatmap = heatmap(China_from_to',xticks=(1:2:17,age_cats[1:2:17]),yticks=(1:2:17,age_cats[1:2:17]),
                    title = "Age-structured contact rates: China (contacts at home)",
                    xlabel = "Contact from individual in age group",
                    ylabel = "Contact received by someone in age group ")
savefig(plt_heatmap,"plotting/China_home_only.pdf")


#Calculate eigenvalues and eigenvectors of the contact matrix. The (L_1) normalised leading eigenvector
#is the case distribution **if**
evals,evects = eigen(China_to_from);

R₀ = Real(evals[end])
v = Real.(evects[:,end])
_v = convert_to_data_agegroups(v)

Case_Dist = LinearAlgebra.normalize(CCDC_Case_Dist_Feb_11[:],1)


plt = groupedbar(hcat(_China_PP,_v,Case_Dist),lab = ["China pop." "Age indep. pred." "Conf. cases"],
        xticks=(1:9,["0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+"]),
        title = "Age profile of confirmed Chinese COVID-19 cases (11/2)",
        ylabel = "Proportion of cases",
        xlabel = "Age of case")
savefig(plt,"plotting/age_profile_cases.pdf")


#This is a struct for the case data in numbers rather than proportions
struct CaseDistribution
    C::Vector{Int64} #Total number of cases in each age group
    n::Int64 #Total number of cases
    prediction::Function #This is the function that makes a prediction about case distribution
end


#This is the log-likelihood function (a callable function on CaseDistribution structs). Its an unnormalised
#independent multinomial likelihood
function (cases::CaseDistribution)(θ) #Susceptibility form of likelihood
    @unpack θ = θ               # extract the parameters θ is a generic label for any set of params, the difference is in the prediction function
    @unpack n, C, prediction = cases   # extract the data and prediction function
    T = eltype(θ)
    p = prediction(θ) #predicted distribution of cases
    try
        logL = 0.
        for i = 1:8
            logL += C[i]*log(p[2*(i-1)+1] + p[2*(i-1)+2])
        end
        logL += C[9]*log(p[17])
        return logL  - 200*(θ[17] - 1)^2 #Prior the last sigma to tether the posterior distribution to a choice
    catch errtype
        return T(-Inf)
    end

end

#Modified prediction function --- this assumes mean 2 days of pre-symptomatic transmission, mean 7 days infectious period for
#mild infecteds and that the eventually severe cases don't transmit more than mild cases

function pred_case_distribution_using_iter_K_model2(χ::Vector,d::Vector,ϵ::Vector,C)
    d1, = size(C)
    K = zeros(eltype(d),d1,d1)
    for a = 1:d1,b=1:d1
        K[a,b] = χ[a]*C[a,b]*(2*ϵ[b] + 7*d[b] + 7*(1-d[b])*ϵ[b])
    end
    v = (K^10)*ones(d1)

    return d.*v/sum(d.*v)
end

#This function returns the HMC chains with tree statistics for the run
function HMC_for_detection_rate(χ,ϵ,n_draws)
    cases_to_fit = CaseDistribution(Int64.(CCDC_Case_Dist_Feb_11[:]),
                        Int64(sum(CCDC_Case_Dist_Feb_11)),
                        d -> pred_case_distribution_using_iter_K_model2(χ,d,ϵ*ones(17),China_to_from))

    trans = as((θ = as(Array, asℝ₊, 17),)) # All parameters are transformed to be positive.
    P = TransformedLogDensity(trans, cases_to_fit) #This creates a transformed log-likelihood
    ∇P = ADgradient(:ForwardDiff, P) #This automatically generates a log-likelihood gradient at the same time as the likelihood is called
    return results = mcmc_with_warmup(Random.GLOBAL_RNG, ∇P, n_draws,reporter = NoProgressReport())
end

χ_zhang = vcat(0.34*ones(3),ones(10),1.47*ones(4))

d_chain_model2_epsilon_0 = HMC_for_detection_rate(χ_zhang,0.,10000)
d_chain_model2_epsilon_01 = HMC_for_detection_rate(χ_zhang,0.1,10000)
d_chain_model2_epsilon_025 = HMC_for_detection_rate(χ_zhang,0.25,10000)
d_chain_model2_epsilon_05 = HMC_for_detection_rate(χ_zhang,0.5,10000)
d_chain_model2_epsilon_1 = HMC_for_detection_rate(χ_zhang,1.,10000)

@save "data/inference_for_age_dependent_symptomatic_rates/HMC_chains_for_model2.jld2" d_chain_model2_epsilon_0 d_chain_model2_epsilon_01 d_chain_model2_epsilon_025 d_chain_model2_epsilon_05 d_chain_model2_epsilon_1


zhang_sus_dist = pred_case_distribution_using_iter_K_model2(χ_zhang,ones(17),ones(17),China_to_from)
_zhang_sus_dist = convert_to_data_agegroups(zhang_sus_dist)
plt = groupedbar(hcat(_zhang_sus_dist,_v,Case_Dist),lab = ["Zhang sus. profile only" "Age indep. pred." "Conf. cases"],
        xticks=(1:9,["0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80+"]),
        title = "Age profile of confirmed Chinese COVID-19 cases (11/2)",
        ylabel = "Proportion of cases",
        xlabel = "Age of case")

@load "data/inference_for_age_dependent_symptomatic_rates/HMC_chains_for_model2.jld2" d_chain_model2_epsilon_0 d_chain_model2_epsilon_01 d_chain_model2_epsilon_025 d_chain_model2_epsilon_05 d_chain_model2_epsilon_1
    trans = as((θ = as(Array, asℝ₊, 17),)); #Posterior draws are saved in transformed state, trans is used to convert back to untransformed draws.

function MCMCchain2estimates(results)
            posterior = first.(transform.(trans, results.chain))
            posterior_array = zeros(length(posterior),17)
            for i = 1:length(posterior),j=1:17
                posterior_array[i,j] = posterior[i][j]
            end
            for i = 1:length(posterior),j=1:17
                posterior_array[i,j] = posterior_array[i,j]/posterior_array[i,17] #Renormalise d_80+ to 1 and other values to relative to this
            end
            θ̂ = [mean(posterior_array[:,i]) for i = 1:17]
            lb = [quantile(posterior_array[:,i],0.025) for i = 1:17]
            ub = [quantile(posterior_array[:,i],0.975) for i = 1:17]
            return θ̂,θ̂ .- lb,ub.-θ̂
        end

d_0,lb_0,ub_0 = MCMCchain2estimates(d_chain_model2_epsilon_0)
d_01,lb_01,ub_01 = MCMCchain2estimates(d_chain_model2_epsilon_01)
d_025,lb_025,ub_025 = MCMCchain2estimates(d_chain_model2_epsilon_025)
d_05,lb_05,ub_05 = MCMCchain2estimates(d_chain_model2_epsilon_05)
d_1,lb_1,ub_1 = MCMCchain2estimates(d_chain_model2_epsilon_1)


x = 5:30:(16*30+5)

fig_d = scatter(x.-8,d_0,yerr = (lb_0,ub_0),fillalpha = 0.3,lab="epsilon = 0",legend = :topleft,
            xticks = (x,["0-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39","40-44","45-49","50-54","55-59","60-64","65-70","70-74","75-79","80+"]),
                ylabel = "Relative detection rate", title = "Detection rate and rel. transmissibility of undetected cases",
                yscale = :log10)
scatter!(fig_d,x.-4,d_01,yerr = (lb_01,ub_01),fillalpha = 0.3,lab="epsilon = 0.1")
scatter!(fig_d,x,d_025,yerr = (lb_025,ub_025),fillalpha = 0.3,lab="epsilon = 0.25")
scatter!(fig_d,x.+4,d_05,yerr = (lb_05,ub_05),fillalpha = 0.3,lab="epsilon = 0.5")
scatter!(fig_d,x.+8,d_1,yerr = (lb_1,ub_1),fillalpha = 0.3,lab="epsilon = 1")
scatter!(size=(700,400),xlabel = "Age group")

savefig(fig_d,"plotting/model_2_detection_rates_old_next_gen.png")
@save "data/inference_for_age_dependent_symptomatic_rates/detection_rates_for_different_epsilons_model2.jld2" d_0 d_01 d_025 d_05 d_1



# Model disgnostics
function posterior_prediction(result,tau)
    posterior = first.(transform.(trans, result.chain))
    pred = map(d -> pred_case_distribution_using_iter_K_model2(χ_zhang,d,tau*ones(17),China_to_from),posterior);
    _pred = map(convert_to_data_agegroups,pred)
    return groupedbar(hcat(Case_Dist,mean(_pred)),lab = ["True distrib." "Post. prediction"],title = "Check prediction for tau =$tau")
end

# plt_pred_05 = groupedbar(hcat(Case_Dist,mean(_pred)),lab = ["True distrib." "Post. prediction"],title = "Check prediction for tau =0.25")
plt1 = posterior_prediction(d_chain_model2_epsilon_0,0);
plt2 = posterior_prediction(d_chain_model2_epsilon_01,0.1);
plt3 = posterior_prediction(d_chain_model2_epsilon_025,0.25);
plt4 = posterior_prediction(d_chain_model2_epsilon_05,0.5);
plt5 = posterior_prediction(d_chain_model2_epsilon_1,1);

layout = @layout [a b;c d;e]
plt = plot(plt1,plt2,plt3,plt4,plt5,layout = layout)
plot!(plt,size = (1000,1000))

summarize_tree_statistics(d_chain_model2_epsilon_1.tree_statistics)

#Group results
results_group = [d_chain_model2_epsilon_0,d_chain_model2_epsilon_01,d_chain_model2_epsilon_025,d_chain_model2_epsilon_05,d_chain_model2_epsilon_1];
eff_sample_sizes = Float64[]
for result in results_group
    posterior = first.(transform.(trans, result.chain)) #Apply transformation to chain
    for a = 1:17
        posterior_samples_for_group_a = [posterior[i][a] for i = 1:length(posterior)]
    push!(eff_sample_sizes,effective_sample_size(posterior_samples_for_group_a)   )
    end
end

histogram(eff_sample_sizes,bins = 100,lab = "",xlabel = "Eff. sample size",
            title = "Frequency of eff. sample size by age group and scenario",
            xticks = (3000:1000:10000),yticks = 0:1:15)
