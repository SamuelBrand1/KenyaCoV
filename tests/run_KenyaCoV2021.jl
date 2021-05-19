

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







#
