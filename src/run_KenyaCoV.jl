
using DifferentialEquations,Plots,DataFrames,Parameters,LinearAlgebra,Distributions
include("kenya_data.jl");
include("types.jl");
include("gravity_model.j")
ρ = 0.01
transmission_matrix = similar(transport_matrix)
for i = 1:n,j = 1:n
    if i != j
        transmission_matrix[i,j] = ρ*transport_matrix[i,j]
    else
        transmission_matrix[i,j] = 0. #within county transmission dealt with seperately
    end
end
KenyaTbl[:County][30]
heatmap(transport_matrix)
transport_matrix[30,28]
transport_matrix[28,30]
KenyaTbl[:Urban][30]/KenyaTbl[:Urban][28]

P = CoVParameters2(T = transmission_matrix,ρ = ρ,δ = 1.,β = 2.)
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
N_urb = sum(u0[:,:,1],dims = 2)
N_rural = sum(u0[:,:,2],dims = 2)

u0[1,4,1] += 100
prob = DiscreteProblem(u0,(0.0,21.0),P)
jump_prob = JumpProblem(prob,Direct(),jump_uu_trans,
                                    jump_wct_trans,
                                    jump_incubation,
                                    jump_recovery,
                                    jump_hosp,
                                    jump_leave_hosp,
                                    save_positions=(false,false))
@time sol = solve(jump_prob,FunctionMap(),saveat = 7.)
# integ = init(jump_prob,FunctionMap(),saveat = 7.)
