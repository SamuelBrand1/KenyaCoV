abstract type AbstractCoVParameters end

"""
    mutable struct CoVParameters_AS

Struct for containing relevant epidemilogical parameters for the age-structured version of KenyaCoV
"""
@with_kw mutable struct CoVParameters_AS <: AbstractCoVParameters
    #Epidemiological parameters and social contact/mixing rates
    β::Float64 = 1. #Basic transmission probability per contact OR infectious contact rate (can be greater than one depending on scaling)
    c_t::Function = t -> 1. #Time varying basic contact rate
    γ::Float64 = 1/7. #recovery rate for mild and asymptomatic cases (Consensus estimate mean 9 days infectious, default is 1/σ₂ + 1/γ = 9 days)
    σ₁::Float64 = 1/3.1 #end of latency rate (Consensus estimate mean 3 days)
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
    mutable struct CoVParameters_HH

Struct for containing relevant epidemilogical parameters for the age-structured version of KenyaCoV with explict
    differences between within and without household transmission
"""
@with_kw mutable struct CoVParameters_HH <: AbstractCoVParameters
    #Number of age, area and state types
    "Number of age types"
    n_age::Int64
    "Number of area types"
    n_area::Int64
    "Number of state types"
    n_state::Int64
    #Epidemiological parameters and social contact/mixing rates
    "Basic transmission probability per contact OR infectious contact rate (can be greater than one depending on scaling) OUTSIDE HH"
    βₒ::Float64 = 1.
    "Basic transmission probability per contact OR infectious contact rate (can be greater than one depending on scaling) INSIDE HH"
    βᵢ::Float64 = 1.
    "Time varying basic contact rate"
    c_t::Function = t -> 1.
    "recovery rate for mild and asymptomatic cases (Consensus estimate mean 9 days infectious, default is 1/σ₂ + 1/γ = 9 days)"
    γ::Float64 = 1/7.
    "end of latency rate (Consensus estimate mean 3 days)"
    σ₁::Float64 = 1/3.1
    "end of pre-symptomatic rate (Consensus estimate mean 4 days infectious, default is 1/σ₂ + 1/γ = 4 days)"
    σ₂::Float64 = 1/2.
    "Proportion of over 80s who are symptomatic"
    δ::Float64 = 0.9
    "Rate of hospitalisation treatment, conditional on eventually needing it (V category)"
    τ::Float64 = 1/5.
    " isolation rate for symptomatics at beginning of epidemic"
    τ_initial::Float64 = 0.
    " Average time to end isolation NOT USED IN THIS VERSION"
    clear_quarantine = 0.
    "Excess mortality due to disease NOT USED IN THIS VERSION"
    μ₁::Float64 = 0.0
    "Relative infectiousness of undetectable/undetected infecteds both pre-symptomatic and asymptomatic"
    ϵ::Float64 = 0.1
    "Relative infectiousness of mild infecteds due to social avoidance"
    ϵ_M::Float64 = 1.
    "Relative infectiousness of mild, then severe, infecteds due to social avoidance --- default is to set so that M and V have same number of infectious contacts per age group"
    ϵ_V::Float64 = 7. /5.
    "relative symptomatic rate"
    symptomatic_rate::Vector{Float64} = ones(n_a)
    "proportion of severe cases if symptomatic"
    hₐ::Vector{Float64} = zeros(n_a)
    "proportion of hospitalised cases that become critical NOT USED IN THIS VERSION"
    ICUₐ::Vector{Float64} = zeros(n_a)
    "relative susceptibility"
    χ::Vector{Float64} = ones(n_a)
    "Time spent outside of area"
    ρ::Vector{Float64} = [0.01 for i in 1:n]
    "Probability distributions of location for mobile individuals"
    T::Matrix{Float64} = zeros(n,n)
    "Age mixing matrix"
    M::Matrix{Float64} = zeros(n_a,n_a)
    "Age mixing matrix --- if only home contacts are made this changes area by area"
    M_ho::Vector{Matrix{Float64}} = [zeros(n_a,n_a) for i = 1:47]
    "Arrival rates of new infectious individuals by area"
    ext_inf_arrivals::Vector{Function} = [t -> 0. for i = 1:47]
    #Control variables
    "This determines if people are still being isolated"
    isolating_detecteds::Bool = false
    "This determines if social distancing and travel restrictions are still in force"
    lockdown::Bool = false
    schools_closed::Bool = true
    before_week_two::Bool = true
    #Calculation variables
    dt::Float64 = 1. #Useful for the non-negative method
    Î::Matrix{Float64} = zeros(n,n_a) #For inplace calculations
    N::Matrix{Int64} = zeros(Int64,n,n_a)
    N̂::Matrix{Float64} = zeros(n,n_a)#For inplace calculations
    λ::Matrix{Float64} = zeros(n,n_a)#For inplace calculations
    λ_loc::Matrix{Float64} = zeros(n,n_a)#For inplace calculations
    dN::Vector{Int64} = zeros(Int64,n*n_a*n_ta)#For inplace calculations
    poi_rates::Vector{Float64} = zeros(n*n_a*n_ta)#For inplace calculations
    dc::SparseMatrixCSC{Int64,Int64} = sparse(zeros(Int64,n*n_a*n_s,n*n_a*n_ta))#For inplace calculations
    du_linear::Vector{Int64} = zeros(Int64,n*n_a*n_s) #for inplace calculations
    index_as::CartesianIndices{3,Tuple{UnitRange{Int64},UnitRange{Int64},UnitRange{Int64}}} = CartesianIndices((1:n, 1:n_a,1:n_s))
    index_as_events::CartesianIndices{3,Tuple{UnitRange{Int64},UnitRange{Int64},UnitRange{Int64}}} = CartesianIndices((1:n, 1:n_a,1:n_ta))
    linear_as::LinearIndices{3,Tuple{UnitRange{Int64},UnitRange{Int64},UnitRange{Int64}}} = LinearIndices((1:n, 1:n_a,1:n_s))
    linear_as_events::LinearIndices{3,Tuple{UnitRange{Int64},UnitRange{Int64},UnitRange{Int64}}} = LinearIndices((1:n, 1:n_a,1:n_ta))
end
