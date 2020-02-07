
push!(LOAD_PATH, "/Users/Sam/GitHub/KenyaCoV/src")
using Plots,Parameters,Distributions,KenyaCoV,DifferentialEquations,StatsPlots,JLD2,FileIO,MAT,RecursiveArrayTools

u0,P,transport_matrix = KenyaCoV.model_ingredients_from_data("data/combined_population_estimates.csv",
                                                             "data/optimal_transition_matrix.jld2",
                                                            "data/optimal_movement_matrix.jld2",
                                                            "data/flight_numbers.csv",
                                                            "data/projected_global_prevelance.csv")

"""
Current estimates from Read and Jewell --- assuming that ascertainment is a decent estimate for symptomatic rate

SCENARIO 1: scenario where asymptomatics are just as infectious as symptomatics, this broadly matches Jewell/Read where
undetected infecteds are just as infectious. No need for real-time growth rate matching.
After the first infected assume no more introductions
"""

P.τ = 0.;#No treatment
P.σ = 1/4; #Latent period mean 4 days
P.γ = 1/1.61; #Fast infectious duration
P.β = 1.94;
P.ϵ = 1.;
P.ext_inf_rate = 0.;
u0[30,3,1] += 1#One asymptomatic in Nairobi
jump_prob_tl = create_KenyaCoV_prob(u0,(0.,365.),P)

function output_infecteds_by_county(sol,i)
    [sum(reshape(u,n,n_s,2)[:,3:4,1:2,],dims=[2,3])[:,1,1] for u in sol.u],false
end

CoV_ens_prob = EnsembleProblem(jump_prob_tl,output_func = output_infecteds_by_county)
sim = solve(CoV_ens_prob,SimpleTauLeaping(),dt = 0.25,trajectories = 500)
forecasts = [VectorOfArray(sim.u[i]) for i = 1:length(sim.u)]
# @save "output/forecasts.jld2" forecasts
times = collect(0:0.25:365.)
total_peaktimes = [0. for i = 1:500]
for i = 1:500
    y = sum(forecasts[i][:,:],dims =1)[:]
    total_peaktimes[i] = times[argmax(y)]
end
@save "output/total_peaktimes.jld2" total_peaktimes

peaktimes_by_county = zeros(500,47)
for i = 1:500,j=1:47
    y = forecasts[i][j,:]
    peaktimes_by_county[i,j] = times[argmax(y)]
end



# β_range = range(1.5*(P.γ),4.5*(P.γ),length = 6)
# data = []
# for (i,β) in enumerate(β_range)
#     println(i)
#     P.β = β
#     jump_prob_tl = create_KenyaCoV_prob(u0,(0.,60.),P)
#     CoV_ens_prob = EnsembleProblem(jump_prob_tl,output_func =rerun_early_extinctions)
#     sim = solve(CoV_ens_prob,SimpleTauLeaping(),dt = 0.25,trajectories = 500)
#     push!(data,(β,sim))
# end


# @save "vary_beta_rho_$(P.ρ).jld2" data
# D_01Nai = load("vary_beta_rho_0.01.jld2")["data"]
# D_001Nai = load("vary_beta_rho_0.001.jld2")["data"]
# # quantile(data[6][2].u,[0.05, 0.25, 0.5, 0.75, 0.95])
# gr()
# p1 = plot()
# for (i,β) in enumerate(β_range)
#     if i == 1
#         boxplot!(p1,data_Mom_001[i][2].u,lab = "Mombassa seed",fc = :red,mc=:red)
#         boxplot!(p1,D_001Nai[i][2].u,lab = "Nairobi seed",fc = :blue,mc=:blue)
#     else
#         boxplot!(p1,data_Mom_001[i][2].u,lab = "",fc = :red,mc=:red)
#         boxplot!(p1,D_001Nai[i][2].u,lab = "",fc = :blue,mc=:blue)
#     end
#
# end
# plot!(p1,xticks = (1.5:2:11.5,round.(collect(β_range)/P.γ,digits = 3 )),
#         yscale = :log10,
#         title = "Cumulative cases by 60 days (rho = $(P.ρ))",legend = :bottomright)
# display(p1)
# p2 = plot()
# for (i,β) in enumerate(β_range)
#     boxplot!(p2,data_Mom_01[i][2].u,lab = "",fc = :red,mc=:red)
#     boxplot!(p2,D_01Nai[i][2].u,lab = "",fc = :blue,mc=:blue)
# end
# plot!(p2,xticks = (1.5:2:11.5,round.(collect(β_range)/P.γ,digits = 3 )),
#         yscale = :log10,
#         title = "Cumulative cases by 60 days (rho = 0.01)",xlabel = "Reproductive ratio")
#
# l = @layout [a;b]
#
# plt = plot(p1,p2,layout = l)
# savefig(plt,"Cum_cases.png")
#
# #Now look at Mombassa seed
# u0[30,3,1] = 0
# u0[28,3,1] +=1
# data_Mom_01 = []
# for (i,β) in enumerate(β_range)
#     println(i)
#     P.β = β
#     jump_prob_tl = create_KenyaCoV_prob(u0,(0.,60.),P)
#     CoV_ens_prob = EnsembleProblem(jump_prob_tl,output_func =rerun_early_extinctions)
#     sim = solve(CoV_ens_prob,SimpleTauLeaping(),dt = 0.25,trajectories = 500)
#     push!(data_Mom_01,(β,sim))
# end
#
# KenyaCoV.transportstructure_params!(P,0.001,transport_matrix)
#
# data_Mom_001 = []
# for (i,β) in enumerate(β_range)
#     println(i)
#     P.β = β
#     jump_prob_tl = create_KenyaCoV_prob(u0,(0.,60.),P)
#     CoV_ens_prob = EnsembleProblem(jump_prob_tl,output_func =rerun_early_extinctions)
#     sim = solve(CoV_ens_prob,SimpleTauLeaping(),dt = 0.25,trajectories = 500)
#     push!(data_Mom_001,(β,sim))
# end


# function rerun_early_extinctions(sol,i)
#     num_cases = sum(reshape(sol[end],n,n_s,2)[:,7:8,:])
#     if num_cases >= 6
#         return (num_cases,false)
#     else
#         return (num_cases,true)
#     end
# end
