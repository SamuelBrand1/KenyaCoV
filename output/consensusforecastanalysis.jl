push!(LOAD_PATH, joinpath(homedir(),"GitHub/KenyaCoV/src"))
using Plots,Parameters,Distributions,DifferentialEquations,JLD2,DataFrames,StatsPlots,FileIO,MAT,RecursiveArrayTools
using Statistics: median, quantile
using LinearAlgebra: eigen

#Load data
@load joinpath(homedir(),"Github/KenyaCoVOutputs/sims_consensus_baseline.jld2") sims_baseline
@load joinpath(homedir(),"Github/KenyaCoVOutputs/sims_consensus_baseline_scaled.jld2") sims_baseline_scaled
@load joinpath(homedir(),"Github/KenyaCoVOutputs/sims_consensus_control.jld2") sims_controls
@load joinpath(homedir(),"Github/KenyaCoVOutputs/sims_consensus_control.jld2") sims_controls_scaled

#
sims_baseline.u[2][50][:,:,1:3]

#First grab incidence
function total_incidence_for_each_sim(sims)
    cum_incidence_for_each_sim = zeros(1000,366)
    for k = 1:1000,t = 1:366
        cum_incidence_for_each_sim[k,t] = sum(sims.u[k][t][:,:,1:3])
    end
    return diff(cum_incidence_for_each_sim,dims = 2)
end
incidence = total_incidence_for_each_sim(sims_baseline)

incidence_for_each_sim = diff(cum_incidence_for_each_sim,dims = 2)
y = median(incidence,dims = 1)
plot(y[:].+1,yscale = :log10)
