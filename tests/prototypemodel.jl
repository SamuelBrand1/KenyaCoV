push!(LOAD_PATH,joinpath(homedir(),"Github/KenyaCoV/src"))
# using Revise
# import KenyaCoV
using DifferentialEquations,ModelingToolkit,Latexify,SparseArrays,StaticArrays,LinearAlgebra,Plots
using DelimitedFiles,JLD2,BenchmarkTools
JLD2.@load("data/agemixingmatrix_Kenya_norestrictions.jld2")
# heatmap(M_Kenya)


const σ = 1/3.1
const γ = 1/2.4
S₀ = [1e6 for a = 1:17]
I₀ = zeros(17)
I₀[5] = 100
N_age = [S₀[a] + I₀[a] for a = 1:17]
N = sum(N_age)
M = repeat(1. ./N_age,1,17).*M_Kenya

evals,evects = eigen((1/γ)*M_Kenya)
const β = 1/Real(evals[end])

# @variables S[1:17] E[1:17] I[1:17] R[1:17]
# λ = β.*γ.*M_Kenya*I

# f_sym = vcat(-R₀*S*λ,R₀*S - σ.*E,σ.*E - γ.*I,γ.*I)
λ₀ = zeros(17)
function AgemixingSEIR(du,u,p,t)
        S = @view u[:,1]
        E = @view u[:,2]
        I = @view u[:,3]
        R = @view u[:,4]

        R₀,N = p
        λ = β.*γ.*(M*(I .+ I₁))
        du[:,1] .= -R₀.*S.*λ
        du[:,2] .= (R₀.*S.*λ) .- σ.*E
        du[:,3] .= σ.*E .- γ.*I
        du[:,4] .= γ.*I
        return nothing
end

u0 = hcat(S₀,zeros(17),I₀,zeros(17))
tspan = (0.,365*1)
p = [2.5,N,zeros(17)]
AgemixingSEIR(similar(u0),u0,p,0.)


prob = ODEProblem(AgemixingSEIR,u0,tspan,p)
@benchmark sol = solve(prob,Rosenbrock23(),p = [3,N])
@benchmark sol = solve(prob,Tsit5(),p = [3,N,λ₀])

plot(sol,vars = collect((2*17 +1 ):(2*17 + 17)))

sol.u[end][(3*17 +1 ):(3*17 + 17)]./N

sys = modelingtoolkitize(prob)
jac = eval(ModelingToolkit.generate_jacobian(sys)[2])
tgrad = eval(ModelingToolkit.generate_tgrad(sys)[2])
sp_pattern =ModelingToolkit.jacobian_sparsity(sys)
J =zeros(68,68)
jac(J,u0,p,0.)
J
AgemixingSEIR_sym = ODEFunction(AgemixingSEIR,jac=jac,tgrad=tgrad,sparsity = sp_pattern)
prob_sym = ODEProblem(AgemixingSEIR_sym,u0,tspan,p)
@benchmark sol2 = solve(prob_sym,Rosenbrock23(),p = [2.5,N])
plot(sol2,vars = collect((3*17 +1 ):(3*17 + 17)))


DiffEqBase.has_jac(AgemixingSEIR_sym)
f.jac(J, uprev, p, t)
uprev = copy(u0)
J = zeros(68,68)
AgemixingSEIR_sym.jac(J,u0,p,0.)
jac(J,u0,[0.,1],0.)
