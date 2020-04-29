push!(LOAD_PATH, joinpath(homedir(),"GitHub/KenyaCoV/src"))
using Plots,Parameters,Distributions,DifferentialEquations,JLD2,DataFrames
using Revise
import KenyaCoV
using LinearAlgebra:eigen,normalize


"""
Load age structured data
"""

u0,P,P_dest = KenyaCoV.model_ingredients_from_data("data/data_for_age_structuredmodel.jld2",
                                            "data/flight_numbers.csv",
                                            "data/projected_global_prevelance.csv")
N = sum(u0[:,:,1],dims = 1)


@load "data/agemixingmatrix_china.jld2" M_China
@load "data/agemixingmatrix_Kenya_norestrictions.jld2" M_Kenya
@load "data/agemixingmatrix_Kenya_homeonly.jld2" M_Kenya_ho


"""
Can adjust β to match a desired R₀ value by evaluating the leading eigenvalue
The idea is to match to the chinese epidemic R₀ -- it will be different in Kenya
"""

P.ϵ = 0.1
P.χ = ones(KenyaCoV.n_a)
R₀_scale = KenyaCoV.calculate_R₀_scale(P)
P.χ = ones(KenyaCoV.n_a)/R₀_scale
P.β = 2*P.γ

KenyaCoV.calculate_R₀(P)
KenyaCoV.calculate_R₀_homeonly(P)

P.dt = 0.25

"""
Can vary the spatial contact assumptions as well
"""
# P.ρ = zeros(20)
# KenyaCoV.transportstructure_params!(P,P.ρ,P_dest)


"""
Run model

"""

u0[KenyaCoV.ind_nairobi_as,5,4] = 5#10 diseased
P.ϵ_D = 1
function ramp_down(t)
    if t < 60.
        return (1-t/60) + 0.5*t/60
    else
        return 0.5
    end
end
P.c_t = ramp_down
P.c_t = t -> 1.
prob = KenyaCoV.create_KenyaCoV_non_neg_prob(u0,(0.,365.),P)
@time sol = solve(prob,FunctionMap(),dt = P.dt)

cum_inc = [sum(sol(t)[:,:,8]) for t = 0.:1.:365]
plot(diff(cum_inc))
plot!(diff(cum_inc).+1,yscale = :log10)


sum(sol[end][:,:,8])/sum(u0)
sum(sol[end][:,:,8])

prob_ode = KenyaCoV.create_KenyaCoV_ode_prob(u0,(0.,365.),P)
@time sol_ode = solve(prob_ode,Tsit5())

times = 0:1:365
I_area = zeros(Int64,20,length(sol.t))
for i = 1:20,(ind,t) in enumerate(sol.t)
    I_area[i,ind] = sum(sol(t)[i,:,3:4])
end
# I = [sum(sol(t)[:,:,3:4]) for t in sol.t]
plt = plot(sol.t,I_area[1,:], lab = 1,xlims = (0.,30),ylims = (0.,100));
for i = 2:20
    plot!(plt,sol.t,I_area[i,:],lab = i);
end
display(plt)

"""
Example of adding an event
"""
function isolation_limit(u,t,integrator) # Event when event_f(u,t) == 0
  integrator.p.isolating_detecteds && sum(u[:,:,8]) < 10000
end
function affect_isolation_limit!(integrator)
  integrator.p.τ = 0
  integrator.p.isolating_detecteds = false
end
cb = DiscreteCallback(isolation_limit,affect_isolation_limit!)

P.isolating_detecteds = true
P.τ = 1/3
@time sol_cb = solve(prob,FunctionMap(),dt = P.dt,callback = cb)
sum(sol_cb[end][:,:,7:8])

"""
Example of adding a constant rate jump to the method
"""

#1. Define the jump rate and affect on integrator
jumprate(u,p,t) = 1.
function jumpaffect!(integrator)
    integrator.u[4,5,3] += 1
end
#2. Declare the ConstantRateJump
example_jump = ConstantRateJump(jumprate,jumpaffect!)
#3. Make a jump problem that BOTH inherits the DiscreteProblem (prob) AND is given an "aggregator" method
#In this case the aggregator method "Direct()" is the Gillespie alg

Jmp_prob = JumpProblem(prob,Direct(),example_jump)

#Solve the whole problem using FunctionMap() this time steps forward discretely BUT if a ConstantRateJump
# event occurs in the interval it includes that as well
sol = solve(Jmp_prob,FunctionMap(),dt = P.dt)
sol.t
