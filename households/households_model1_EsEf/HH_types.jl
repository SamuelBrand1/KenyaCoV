

@with_kw mutable struct CoVParameters
    β::Float64 = 2.5/3.6
    γ::Float64 = 1/3.6
    σ::Float64 = 1/2.
    δ::Float64 = 0.05#Proportion of symptomatic/diseased vs non-symptomatic cases
    τ::Float64 = 1/15. #treatment/isolation rate for symptomatics
    μ₁::Float64 = 0.01#Excess mortality due to disease
    ϵ::Float64 = 0.1 #Relative infectiousness of undetectable infecteds
    ρ::Vector{Float64} = [0.01 for i in 1:n] #spatial coupling
    T::Matrix{Float64} = zeros(n,n)#transmission matrix
    ext_inf_rate::Float64 = 0. #Scales the contact rate with the infecteds arriving via air
    into_mom::Vector{Int}#Number of people flying into Mombassa each day
    into_nai::Vector{Int}#Number of people flying into Nairobi each day
    global_prev::Vector{Float64}
    dt::Float64 = 1. #Useful for the non-negative method
    Î::Vector{Float64} = zeros(n) #For inplace calculations
    N̂::Vector{Float64} = zeros(n)#For inplace calculations
    λ_urb::Vector{Float64} = zeros(n)#For inplace calculations
    λ_rur::Vector{Float64} = zeros(n)#For inplace calculations
    dN::Vector{Int64} = zeros(Int64,n*n_t)#For inplace calculations
    poi_rates::Vector{Float64} = zeros(n*n_t)#For inplace calculations
    dc::SparseMatrixCSC{Int64,Int64} = sparse(zeros(Int64,n*n_s*2,n*n_t))#For inplace calculations
end

@with_kw mutable struct CoVParameters_AS
    #β::Float64 = 1.            #HH
    βˢ::Float64 = 1.            #HH β subsequent (or β household)
    βᶠ::Float64 = 1.            #HH β first (or β other)
    c_t::Function = t -> 1. #Time varying basic contact rate
    γ::Float64 = 1/2.5
    σ::Float64 = 1/5.
    δ::Float64 = 0.9#Proportion of over 80s who get identified
    τ::Float64 = 0. #current isolation rate for symptomatics
    τ_initial::Float64 = 0. # isolation rate for symptomatics at beginning of epidemic
    clear_quarantine = 1/14. # Two weeks on average to end isolation
    μ₁::Float64 = 0.0#Excess mortality due to disease
    ϵ::Float64 = 0.1 #Relative infectiousness of undetectable infecteds
    ϵ_D::Float64 = 1.#Relative infectiousness of detectable infecteds after intervention
    rel_detection_rate::Vector{Float64} = ones(n_a) #relative detection rate
    χ::Vector{Float64} = ones(n_a) #relative susceptibility
    ρ::Vector{Float64} = [0.01 for i in 1:n] #Time spent outside of area
    T::Matrix{Float64} = zeros(n,n)#Probability distributions of location for mobile individuals
    M::Matrix{Float64} = zeros(n_a,n_a) #Age mixing matrix
    M_ho::Matrix{Float64} = zeros(n_a,n_a) #Age mixing matrix --- if only home contacts are made
    ext_inf_rate::Float64 = 0. #Scales the contact rate with the infecteds arriving via air
    into_mom::Vector{Int}#Number of people flying into Mombassa each day
    into_nai::Vector{Int}#Number of people flying into Nairobi each day
    global_prev::Vector{Float64}
    isolating_detecteds::Bool = true #This determines if people are still being isolated
    lockdown::Bool = true  #This determines if social distancing and travel restrictions are still in force
    dt::Float64 = 1. #Useful for the non-negative method
    Î::Matrix{Float64} = zeros(n_wa,n_a) #For inplace calculations
    N̂::Matrix{Float64} = zeros(n_wa,n_a)#For inplace calculations

    # λ::Matrix{Float64} = zeros(n_wa,n_a)#For inplace calculations     #HH
    λᶠ::Matrix{Float64} = zeros(n_wa,n_a)                               #HH
    λˢ::Matrix{Float64} = zeros(n_wa,n_a)                               #HH

    # λ_loc::Matrix{Float64} = zeros(n_wa,n_a)#For inplace calculations #HH
    λᶠ_loc::Matrix{Float64} = zeros(n_wa,n_a)                           #HH
    λˢ_loc::Matrix{Float64} = zeros(n_wa,n_a)                           #HH

    dN::Vector{Int64} = zeros(Int64,n_wa*n_a*n_ta)#For inplace calculations
    poi_rates::Vector{Float64} = zeros(n_wa*n_a*n_ta)#For inplace calculations
    dc::SparseMatrixCSC{Int64,Int64} = sparse(zeros(Int64,n_wa*n_a*n_s,n_wa*n_a*n_ta))#For inplace calculations
    du_linear::Vector{Int64} = zeros(Int64,n_wa*n_a*n_s) #for inplace calculations
end
