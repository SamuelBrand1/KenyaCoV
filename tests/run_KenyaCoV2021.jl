

push!(LOAD_PATH, "./src")
    using DifferentialEquations,Parameters,DataFrames, JLD2, FileIO, LinearAlgebra, CSV, NamedArrays, SparseArrays, Revise, Plots
    #using Distributions,RecursiveArrayTools
    import KenyaCoV



@time sims = KenyaCoV.run_sims()
sims.u[1].incidence_I
plot(sum(sims.u[1].incidence_I,dims=1)[1,:])
    plot!(sum(sims.u[1].incidence_A,dims=1)[1,:])
    plot!(sum(sims.u[1].incidence_M,dims=1)[1,:])
    plot!(sum(sims.u[1].incidence_V,dims=1)[1,:])
sum(sims.u[1].incidence_I)


β₀ = 1.
@load("data/baseline_fitted_parameters.jld2")
p_baseline.

M_h,M_o,M_s,M_w = KenyaCoV.get_rescaled_contact_matrices()
N_kenya = KenyaCoV.get_population_size_matrix()

I_A = rand(47,17)
I_P = rand(47,17)
I_M = rand(47,17)
I_V = rand(47,17)
# λ::Matrix{Float64} = zeros(n,n_a)
(p_baseline.β₀   .* (p_baseline.β_home.*M_h .+ p_baseline.β_other.*M_o .+ p_baseline.β_work.*M_w .+ p_baseline.β_school.*M_s) )  .*((p_baseline.ϵ.*I_A .+ I_P .+ I_M .+ I_V)./permutedims(N_kenya))

A=rand(17,17)
B=rand(47,17)
A*B
dc = sparse(zeros(Int64,KenyaCoV.n*KenyaCoV.n_a*KenyaCoV.n_s,KenyaCoV.n*KenyaCoV.n_a*KenyaCoV.n_ta))
KenyaCoV.change_matrix(dc)

I_A./permutedims(N_kenya)

#
