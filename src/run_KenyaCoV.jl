push!(LOAD_PATH, "/Users/Sam/GitHub/KenyaCoV/src")
using DifferentialEquations,Plots,DataFrames,Parameters,LinearAlgebra,Distributions,Distributed
include("kenya_data.jl");
include("types.jl");
include("gravity_model.jl");
ρ = 0.001
location_matrix = similar(transport_matrix)

for i = 1:n,j = 1:n
    if i != j
        location_matrix[i,j] = ρ*transport_matrix[i,j]
    else
        location_matrix[i,j] = 1-ρ
    end
end


P = CoVParameters(T = location_matrix,ρ = ρ,β = 2.5/3.6,τ =0.
                ,Î=zeros(n),N̂=zeros(n),λ_urb=zeros(n),λ_rur = zeros(n) )
"""
States:
1 -> S
2 -> E
3 -> I_subclinical
4 -> I_diseased
5 -> H(ospitalised)
6 -> Recovered
7 -> Cumulative I_sub
8 -> Cumulative I_dis
9 -> Cumulative Dead
"""
u0 = zeros(Int64,n,9,2) #Array by area, state and urban vs rural
for i = 1:n
    u0[i,1,1] = KenyaTbl[:Urban][i]
    if !ismissing(KenyaTbl[:Rural][i])
        u0[i,1,2] = KenyaTbl[:Rural][i]
    end
end
N = sum(u0)
N_urb = [sum(u0[i,:,1]) for i = 1:n]
N_rural = [sum(u0[i,:,2]) for i = 1:n]
N̂ = location_matrix*N_urb + N_rural


u0[30,3,1] += 1#One asymptomatic in Nairobi

# include("events.jl");
# prob = DiscreteProblem(u0,(0.0,60.0),P)
# jump_prob = JumpProblem(prob,DirectFW(),jump_urb_trans,
#                                     jump_rural_trans,
#                                     jump_incubation,
#                                     jump_recovery,
#                                     jump_hosp,
#                                     jump_death,
#                                     save_positions=(false,false))
# @time sol = solve(jump_prob,FunctionMap(),saveat = 1.)
# # integ = init(jump_prob,FunctionMap(),saveat = 7.)
# sol[:,end][30,:,1]
#
# CoVensemble_prob = EnsembleProblem(jump_prob)
# # addprocs(3)
# CoVensemble = solve(CoVensemble_prob,FunctionMap(),EnsembleThreads(),trajectories = 100)
include("regularjumps.jl");
u0_vec = u0[:]
# reshape(u0_vec,n,n_s,2)  #Reverse operation vector -> array form

prob_tl = DiscreteProblem(u0_vec,(0.,365.),P)
jump_prob_tl = JumpProblem(prob_tl,Direct(),reg_jump)
@time sol_tl = solve(jump_prob_tl,SimpleTauLeaping();dt = 1.)
sol_tl.u
ũ = [reshape(u,n,n_s,2) for u in sol_tl.u]
susceptibles = [sum(u[:,1,:])/sum(u) for u in ũ ]
infecteds_A = [sum(u[:,3,:]) for u in ũ ]
infecteds_D = [sum(u[:,4,:]) for u in ũ ]
recovereds = [sum(u[:,5,:])/sum(u) for u in ũ ]
cum_infecteds =  [sum(u[:,7:8,:])/sum(u) for u in ũ ]

plot(sol_tl.t,susceptibles,lab="S")
plot(sol_tl.t,infecteds_A,lab ="I_A")
plot!(sol_tl.t,infecteds_D,lab ="I_D")
plot!(sol_tl.t,recovereds,lab="R")
plot!(sol_tl.t,infecteds_D,lab ="I_D")
plot!(sol_tl.t,cum_infecteds,lab ="cum. I")
