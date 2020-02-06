
push!(LOAD_PATH, "/Users/Sam/GitHub/KenyaCoV/src")
using Plots,Parameters,Distributions,KenyaCoV,DifferentialEquations,StatsPlots,JLD2,FileIO

KenyaCoV.transportstructure_params!(P,0.01,transport_matrix)
P.τ = 0.
P.β = 2.5*(P.γ)
u0[30,3,1] += 1#One asymptomatic in Nairobi
jump_prob_tl = create_KenyaCoV_prob(u0,(0.,60.),P)
@time sol = solve(jump_prob_tl,SimpleTauLeaping(),dt = 0.25)


function rerun_early_extinctions(sol,i)
    num_cases = sum(reshape(sol[end],n,n_s,2)[:,7:8,:])
    if num_cases >= 6
        return (num_cases,false)
    else
        return (num_cases,true)
    end
end



# CoV_ens_prob = EnsembleProblem(jump_prob_tl,output_func =rerun_early_extinctions)
# sim = solve(CoV_ens_prob,SimpleTauLeaping(),dt = 0.25,trajectories = 500)
β_range = range(1.5*(P.γ),4.5*(P.γ),length = 6)
data = []
for (i,β) in enumerate(β_range)
    println(i)
    P.β = β
    jump_prob_tl = create_KenyaCoV_prob(u0,(0.,60.),P)
    CoV_ens_prob = EnsembleProblem(jump_prob_tl,output_func =rerun_early_extinctions)
    sim = solve(CoV_ens_prob,SimpleTauLeaping(),dt = 0.25,trajectories = 500)
    push!(data,(β,sim))
end


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
