push!(LOAD_PATH,joinpath(homedir(),"Github/KenyaCoV/src"))
# using Revise
# import KenyaCoV
using DifferentialEquations,ModelingToolkit,Latexify,SparseArrays,StaticArrays,LinearAlgebra,Plots
using DelimitedFiles,JLD2
JLD2.@load("data/agemixingmatrix_Kenya_norestrictions.jld2")
heatmap(M_Kenya)

# const σ = 1/3.1
# const γ = 1/2.4
# evals,evects = eigen((1/γ)*M_Kenya)
# const β = 1/Real(evals[end])

@parameters t R₀
@variables S[1:17](t) E[1:17](t) I[1:17](t) R[1:17](t)
@derivatives D'~t
λ = β.*γ.*M_Kenya*I

@parameters t σ ρ β
@variables x(t) y(t) z(t)
@derivatives D'~t
eqs = [D(x) ~ σ*(y-x),
       D(y) ~ x*(ρ-z)-y,
       D(z) ~ x*y - β*z]
de = ODESystem(eq)

D([S,E]) ~ [-R₀*S.*λ, R₀*S.*λ - σ.*E]
eqs = vcat(D(S) ~ -R₀*S,
        D(E) ~ R₀*S - σ.*E,
        D(I) ~ σ.*E - γ.*I,
        D(R) ~ γ.*I)
eqs = [D(S) ~ -R₀*S  ]

f_expr = build_function(eqs)

de = ODESystem(eqs)
ModelingToolkit.calculate_jacobian(de)
u0 = [S => 1e6*ones(17),
      E => zeros(17),
      I => 100*ones(17),
      R => zeros(17)]

p  = [R₀ => 1.50]
tspan = (0.0,200.0)

prob = ODEProblem(de,u0,tspan,p)
