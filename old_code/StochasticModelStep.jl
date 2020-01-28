# Initialise and define the stochastic model


function DistributeEvents(rates::Vector{Float64},MaxNumEvents::Vector{Int64})
    return [min(rand(Poisson(rates[i])),MaxNumEvents[i]) for i = 1:length(rates)]
end


#Initialise matrix --- useful for using the DifferentialEquations plotting add-ons etc
#Rows are S1,I1,R1,S2,I2,R2, cols are different counties


function EpidemicEvents!(t,u,p::VariableRSVParameters)
    du = zeros(Int64,n,6)
    #Rates
    N = sum(u,2)
    N_at = (1-p.ρ)*N + p.ρ*F*N
    χ = [exp(p.ξ*cospi(2*(t-KenyaTbl[:Phase][i])/365.25)) for i = 1:n]
    I_at = inf_2*((1-p.ρ)*u[:,5] + p.ρ*F*u[:,5]) + ((1-p.ρ)*u[:,2] + p.ρ*F*u[:,2])
    I_at = χ.*I_at
    λ = p.β*(1-p.ρ)*(I_at./N_at) + p.β*p.ρ*((I_at'./N_at')*F)'
    #Implement events
    #Births
    BirthRates = [ Δt*BR(t)*KenyaTbl[:BirthRate_2014][i]*N[i] for i = 1:n]
    Births = DistributeEvents(BirthRates,1000*ones(Int64,length(BirthRates)))
    du[:,1] += Births
    #Deaths
    DeathRates = [ Δt*μ0*μ(t)*u[i,j] for j = 1:6 for i = 1:n]
    Deaths = DistributeEvents(DeathRates,u[:])
    Deaths = reshape(Deaths,n,6)
    du -= Deaths
    #Recoveries
    Recs1 = DistributeEvents(Δt*γ1*u[:,2],u[:,2])
    Recs2 = DistributeEvents(Δt*γ2*u[:,5],u[:,5])
    du[:,2] -= Recs1
    du[:,3] += Recs1
    du[:,5] -= Recs2
    du[:,6] += Recs2
    #Reversions
    Revs1 = DistributeEvents(Δt*ν*u[:,3],u[:,3])
    Revs2 = DistributeEvents(Δt*ν*u[:,6],u[:,6])
    du[:,3] -= Revs1
    du[:,4] += Revs1
    du[:,6] -= Revs2
    du[:,4] += Revs2
    #incidence
    InfRate1 = [Δt*u[i,1]*λ[i] for i = 1:n]
    InfRate2 = [Δt*sus_2*u[i,4]*λ[i] for i = 1:n]
    Infs1 = DistributeEvents(InfRate1,u[:,1])
    Infs2 = DistributeEvents(InfRate2,u[:,4])
    du[:,1] -= Infs1
    du[:,2] += Infs1
    du[:,4] -= Infs2
    du[:,5] += Infs2
    u = u + du
end



function EquilibriumEpidemicEvents!(t,u,p::VariableRSVParameters)
    du = zeros(Int64,n,6)
    #Rates
    N = sum(u,2)
    N_at = (1-p.ρ)*N + p.ρ*F*N
    χ = [exp(p.ξ*cospi(2*(t-KenyaTbl[:Phase][i])/365.25)) for i = 1:n]
    I_at = inf_2*((1-p.ρ)*u[:,5] + p.ρ*F*u[:,5]) + ((1-p.ρ)*u[:,2] + p.ρ*F*u[:,2])
    I_at = χ.*I_at
    λ = p.β*(1-p.ρ)*(I_at./N_at) + p.β*p.ρ*((I_at'./N_at')*F)'
    #Implement events
    #Births
    BirthRates = [ Δt*μ0*KenyaTbl[:Total][i] for i = 1:n]
    Births = DistributeEvents(BirthRates,1000*ones(Int64,length(BirthRates)))
    du[:,1] += Births
    #Deaths
    DeathRates = [ Δt*μ0*u[i,j] for j = 1:6 for i = 1:n]
    Deaths = DistributeEvents(DeathRates,u[:])
    Deaths = reshape(Deaths,n,6)
    du -= Deaths
    #Recoveries
    Recs1 = DistributeEvents(Δt*γ1*u[:,2],u[:,2])
    Recs2 = DistributeEvents(Δt*γ2*u[:,5],u[:,5])
    du[:,2] -= Recs1
    du[:,3] += Recs1
    du[:,5] -= Recs2
    du[:,6] += Recs2
    #Reversions
    Revs1 = DistributeEvents(Δt*ν*u[:,3],u[:,3])
    Revs2 = DistributeEvents(Δt*ν*u[:,6],u[:,6])
    du[:,3] -= Revs1
    du[:,4] += Revs1
    du[:,6] -= Revs2
    du[:,4] += Revs2
    #incidence
    InfRate1 = [Δt*u[i,1]*λ[i] for i = 1:n]
    InfRate2 = [Δt*sus_2*u[i,4]*λ[i] for i = 1:n]
    Infs1 = DistributeEvents(InfRate1,u[:,1])
    Infs2 = DistributeEvents(InfRate2,u[:,4])
    du[:,1] -= Infs1
    du[:,2] += Infs1
    du[:,4] -= Infs2
    du[:,5] += Infs2
    u = u + du
end
