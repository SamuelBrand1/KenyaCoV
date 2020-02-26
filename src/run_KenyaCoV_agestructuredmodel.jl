push!(LOAD_PATH, joinpath(homedir(),"GitHub/KenyaCoV/src"))
using Plots,Parameters,Distributions,DifferentialEquations,JLD2,DataFrames
using Revise
import KenyaCoV
using LinearAlgebra:eigen


"""
Load age structured data
"""

u0,P,P_dest = KenyaCoV.model_ingredients_from_data("data/data_for_age_structuredmodel.jld2",
                                            "data/flight_numbers.csv",
                                            "data/projected_global_prevelance.csv")
"""
Can adjust β to match a desired R₀ value by evaluating the leading eigenvalue
The idea is to match to the chinese epidemic -- it might be different in Kenya
"""
M_china_tbl = readtable("data/china_baseline_age_mixing.csv")
M_china = zeros(16,16)
for i = 1:16,j=1:16
    M_china[i,j] = M_china_tbl[i,j]
end
P.ϵ = 1
sus_matrix = repeat(P.χ,1,16)
eigs_china, = eigen((P.δ + (1-P.δ)*P.ϵ)*sus_matrix.*M_china)
max_eigval_china = Real(eigs_china[end])
P.β = 2.2*P.γ
eigs, = eigen((P.β/P.γ)*sus_matrix.*P.M)
# @time KenyaCoV.rates(P.poi_rates,u0,P,1.)
#dc = P.dc
# @time KenyaCoV.nonneg_tauleap(deepcopy(u0),u0,P,0.)

u0[KenyaCoV.ind_nairobi_as,15,3] += 100#10 diseased

prob = KenyaCoV.create_KenyaCoV_non_neg_prob(u0,(0.,365.),P)
@time sol = solve(prob,FunctionMap(),dt = 0.25)
sum(sol[end][:,:,7:8])

prob_ode = KenyaCoV.create_KenyaCoV_ode_prob(u0,(0.,365.),P)
@time sol_ode = solve(prob_ode,Tsit5())

times = 0:1:365
I_area = zeros(Int64,20,length(sol.t))
for i = 1:20,(ind,t) in enumerate(sol.t)
    I_area[i,ind] = sum(sol(t)[i,:,3:4])
end
# I = [sum(sol(t)[:,:,3:4]) for t in sol.t]
plt = plot(sol.t,I_area[1,:], lab = 1,xlims = (0.,180));
for i = 2:20
    plot!(plt,sol.t,I_area[i,:],lab = i);
end
display(plt)
