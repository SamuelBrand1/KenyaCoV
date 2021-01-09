push!(LOAD_PATH,joinpath(homedir(),"Github/KenyaCoV/src"))
# using Revise
# import KenyaCoV
using DifferentialEquations,ModelingToolkit,Latexify,SparseArrays,StaticArrays,LinearAlgebra,Plots
using DelimitedFiles,JLD2,Interpolations
JLD2.@load("data/agemixingmatrix_Kenya_norestrictions.jld2")
# heatmap(M_Kenya)
M_Kenya[2,1]
JLD2.@load("data/data_for_age_structuredmodel_with_counties.jld2")
T
const σ_inc = 1/3.1
const γ_rec = 1/2.4
S₀ = [1e6/17 for i = 1:47,a = 1:17]
I₀ = zeros(47,17)
I₀[30,5] = 1000.
N_mat = S₀ .+ I₀
N̂ = T*N_mat
movements = T'*T


evals,evects = eigen((1/γ_rec)*M_Kenya)
const β = 1/Real(evals[end])

# @variables S[1:17] E[1:17] I[1:17] R[1:17]
# λ = β.*γ.*M_Kenya*I

# f_sym = vcat(-R₀*S*λ,R₀*S - σ.*E,σ.*E - γ.*I,γ.*I)
function calculate_forceofinfection(β,γ,T,M,I)
        Î = (T*I)./N̂
        β*γ*T'*Î*M
end
function AgemixingSEIR(du,u,p,t)
        S = @view u[:,:,1]
        E = @view u[:,:,2]
        I = @view u[:,:,3]
        R = @view u[:,:,4]
        R₀,N = p
        # λ  = β*γ*(M*I)
        λ = calculate_forceofinfection(β,γ_rec,T,M_Kenya,I)
        @. du[:,:,1] = -R₀*S*λ
        @. du[:,:,2] = (R₀*S*λ) - σ_inc*E
        @. du[:,:,3] = σ_inc*E - γ_rec*I
        @. du[:,:,4] = γ_rec*I
        return nothing
end
AgemixingSEIR(similar(u0),u0,p,0.)

u0 = cat(S₀,zeros(47,17),I₀,zeros(47,17),dims = 3)
tspan = (0.,365*1)
p = [2.5,N]
prob = ODEProblem(AgemixingSEIR,u0,tspan,p)
@time sol = solve(prob,AutoTsit5(Rosenbrock23()),p = [1.2,N])

plot(sol.t,[sum(u[:,:,3]) for u in sol.u])
plot(sol.t,[sum(u[:,:,4]) for u in sol.u]./N)


sys = modelingtoolkitize(prob)
jac = eval(ModelingToolkit.generate_jacobian(sys)[2])
tgrad = eval(ModelingToolkit.generate_tgrad(sys)[2])
sp_pattern =ModelingToolkit.jacobian_sparsity(sys)
J =zeros(68,68)
jac(J,u0,p,0.)
J
AgemixingSEIR_sym = ODEFunction(AgemixingSEIR,jac=jac,tgrad=tgrad,jac_prototype = similar(J),sparsity = sp_pattern)
prob_sym = ODEProblem(AgemixingSEIR_sym,u0,tspan,p)
@time sol2 = solve(prob_sym,AutoTsit5(Rosenbrock23()),p = [2.5,N])
plot(sol2,vars = collect((2*17 +1 ):(2*17 + 17)))


DiffEqBase.has_jac(AgemixingSEIR_sym)
f.jac(J, uprev, p, t)
uprev = copy(u0)
J = zeros(68,68)
AgemixingSEIR_sym.jac(J,u0,p,0.)
jac(J,u0,[0.,1],0.)
