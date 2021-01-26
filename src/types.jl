"""
    mutable struct CoVParameters

Struct for containing relevant epidemilogical parameters for the age-structured version of KenyaCoV
"""
Base.@kwdef mutable struct CoVParameters_AS
    #Epidemiological parameters and social contact/mixing rates
    β::Float64 = 1. #Basic transmission probability per contact OR infectious contact rate (can be greater than one depending on scaling)
    c_t::Function = t -> 1. #Time varying basic contact rate
    γ::Float64 = 1/7. #recovery rate for mild and asymptomatic cases (Consensus estimate mean 9 days infectious, default is 1/σ₂ + 1/γ = 9 days)
    σ₁::Float64 = 1/3. #end of latency rate (Consensus estimate mean 3 days)
    σ₂::Float64 = 1/2. #end of pre-symptomatic rate (Consensus estimate mean 4 days infectious, default is 1/σ₂ + 1/γ = 4 days)
    δ::Float64 = 0.9#Proportion of over 80s who get identified
    τ::Float64 = 1/5. #Rate of hospitalisation treatment, conditional on eventually needing it (V category)
    τ_initial::Float64 = 0. # isolation rate for symptomatics at beginning of epidemic
    clear_quarantine = 0. # Average time to end isolation NOT USED IN THIS VERSION
    μ₁::Float64 = 0.0#Excess mortality due to disease NOT USED IN THIS VERSION
    ϵ::Float64 = 0.1 #Relative infectiousness of undetectable/undetected infecteds both pre-symptomatic and asymptomatic
    ϵ_D::Float64 = 1.#Relative infectiousness of mild infecteds due to social avoidance
    ϵ_V::Float64 = 7. /5.#Relative infectiousness of mild, then severe, infecteds due to social avoidance --- default is to set so that M and V have same number of infectious contacts per age group
    rel_detection_rate::Vector{Float64} = ones(n_a) #relative symptomatic rate
    hₐ::Vector{Float64} = zeros(n_a)   #proportion of severe cases if symptomatic
    ICUₐ::Vector{Float64} = zeros(n_a) #proportion of hospitalised cases that become critical
    χ::Vector{Float64} = ones(n_a) #relative susceptibility
    ρ::Vector{Float64} = [0.01 for i in 1:n] #Time spent outside of area
    T::Matrix{Float64} = zeros(n,n)#Probability distributions of location for mobile individuals
    M::Matrix{Float64} = zeros(n_a,n_a) #Age mixing matrix
    M_ho::Matrix{Float64} = zeros(n_a,n_a) #Age mixing matrix --- if only home contacts are made
    ext_inf_rate::Float64 = 0. #Scales the contact rate with the infecteds arriving via air
    into_mom::Vector{Int}#Number of people flying into Mombassa each day
    into_nai::Vector{Int}#Number of people flying into Nairobi each day
    global_prev::Vector{Float64}#Chance that each person arriving was infected on that day
    #Control variables
    isolating_detecteds::Bool = false #This determines if people are still being isolated
    lockdown::Bool = false  #This determines if social distancing and travel restrictions are still in force
    schools_closed::Bool = true
    before_week_two::Bool = true
    #Calculation variables
    dt::Float64 = 1. #Useful for the non-negative method
    Î::Matrix{Float64} = zeros(n,n_a) #For inplace calculations
    N̂::Matrix{Float64} = zeros(n,n_a)#For inplace calculations
    λ::Matrix{Float64} = zeros(n,n_a)#For inplace calculations
    λ_loc::Matrix{Float64} = zeros(n,n_a)#For inplace calculations
    dN::Vector{Int64} = zeros(Int64,n*n_a*n_ta)#For inplace calculations
    poi_rates::Vector{Float64} = zeros(n*n_a*n_ta)#For inplace calculations
    dc::SparseMatrixCSC{Int64,Int64} = sparse(zeros(Int64,n*n_a*n_s,n*n_a*n_ta))#For inplace calculations
    du_linear::Vector{Int64} = zeros(Int64,n*n_a*n_s) #for inplace calculations
end

"""
    mutable struct CoVParameters_ODE

Struct for containing relevant epidemilogical parameters for the age-structured ODE version of KenyaCoV
"""
Base.@kwdef mutable struct CoVParameters_ODE
    #Epidemiological parameters and social contact/mixing rates
    β₀::Float64 = 3. #Basic transmission probability per contact OR infectious contact rate (can be greater than one depending on scaling)
    c_t::Function = t -> 1. #Time varying basic contact rate for other activities
    c_t_home::Function = t -> 1. #Time varying basic contact rate at home
    c_t_schools::Function = t -> 1. #Time varying basic contact rate at schools
    γA::Float64 = 1/2.4 #recovery rate for asymptomatic cases
    γM::Float64 = 1/2.4 #recovery rate for mild cases
    γV::Float64 = 1/2.4 #hospitalisation rate for severe cases
    α::Float64 = 1/3.1 #Incubation rate to become infectious
    αP::Float64 = 1/2. #Incubation rate of pre-symptomatic transmissible phase
    ϵ::Float64 = 0.4 #Relative infectiousness of undetectable/undetected infecteds both pre-symptomatic and asymptomatic
    δ::Vector{Float64} = ones(16) #Symptomatic rate by age
    σ::Vector{Float64} = ones(16) #Relative susceptibility rate by age

    M_ho::Matrix{Float64} = zeros(16,16) #Age mixing matrix --- if only home contacts are made
    M_other::Matrix{Float64} = zeros(16,16) #Age mixing matrix --- if only home contacts are made
    M_schools::Matrix{Float64} = zeros(16,16) #Age mixing matrix --- if only home contacts are made

    # ext_inf_rate::Float64 = 0. #Scales the contact rate with the infecteds arriving via air
    # into_mom::Vector{Int}#Number of people flying into Mombassa each day
    # into_nai::Vector{Int}#Number of people flying into Nairobi each day
    # global_prev::Vector{Float64}#Chance that each person arriving was infected on that day
end
