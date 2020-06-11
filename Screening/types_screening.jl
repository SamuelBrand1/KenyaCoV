"""
    mutable struct CoVParameters_Screening

Struct for containing relevant epidemilogical parameters for the age-structured version of KenyaCoV + screening and contact tracing interventions
"""
@with_kw mutable struct CoVParameters_Screening
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



    #=Screening intervention performed:
        -1 -> No intervention
        1->Contact tracing of hospitalized only
        2->Symptomatic screening
        3->Symptomatic screening with contact tracing
        4->Mass screening
        5->Mass screening with contact tracing=#
    dN_S::Vector{Int64} = zeros(Int64,n*n_a*n_ta)#For inplace calculations same as dN but for the quaratining events

    #Contact tracing parameters
    CT_dur::Int64=14    #how  far back do we look for contacts
    C₁_set::Bool=false
    C₁::Array{Float64,2}=zeros(n_a,n_a)#M .* CT_dur #for inplace calculations of daily number of contacteds of age aj to row ai during all CT_dur
    C₂::Array{Float64,2}=zeros(n_a,n_a) #for inplace calculations of transposed(newVH[r,:])*C₁
    N_a::Vector{Float64} =zeros(n_a) #for inplace calculations of number of total contacteds per age by all detected V for all CT_dur
    C₄::Array{Int64,2} =zeros(n_a,CT_dur) #for inplace calculations of number of total contacteds per age by all detected V for each day of CT_dur #CT_dur[1]=t-1=yesterday  and CT_dur[14]=t-14=14 days ago
    CT_n::Int64=100    #mean number of traced contacts per detected
    detected::Array{Int64,2} =zeros(n,n_a) #Number of individuals detected per region and age

    #Mass screening parameters:
    S_strategy::Array{Int64,2}=zeros(Int64,n,370)    #number of tests performed per day per region
    screening_delay::Int64=2  #number of days the tests results are delayed
    selection_Pa::Array{Float64,3}=zeros(370,n,n_a)  # for inplace calculations of probability of an individual of age a to be selected (for testing) each day, for a given location (∑selection_Pa[d,r,:] = 1)
    test_sensitivity::Float64=.0    #proba(I|I)     #drawn from distribution?
    test_specificity::Float64=.0    #proba(!I|!I)   #drawn from distribution?
    Sympt_scr_strategies::Vector{Float64}=zeros(370)

    #For interventions based on history
    selection_P::Array{Float64,4}=zeros(370,n,n_a,8)  # for inplace calculations of probability of an individual of state s to be selected (for testing or as a contact) each day, for a given location and age (∑selection_P[d,r,a,:] = 1)
    selection_N::Array{Int64,2} = zeros(Int64,CT_dur,8#=n_s=#)  #for inplace calculations of number of contacted per age, date and state in a specific region #selection_N[a,day,s]
    state_mvt_P::Array{Float64,4}=zeros(370,n,n_a,15)    # for inplace calculations of numbers/probabilities of movements between states
    toQ::Array{Int64,3}=zeros(Int64,n,n_a,8#=n_s=#)  #number of individuals to quarantine
    selection_P_symptomatics::Array{Float64,4}=zeros(370,n,n_a,8)  # for inplace calculations of probability of an individual of state M or V to be selected (for testing) each day, for a given location and age (∑selection_P[d,r,a,1:2] = 1)
    selection_Pa_symptomatics::Array{Float64,3}=zeros(370,n,n_a)  # for inplace calculations of probability of a symptomatic individual of age a to be selected (for testing) each day, for a given location (∑selection_Pa[d,r,:] = 1)
    multinomial_vector2::Vector{Int64}=zeros(Int64,2) # for inplace calculations of multinomial distribution in update_states!
    multinomial_vector4::Vector{Int64}=zeros(Int64,4) # for inplace calculations of multinomial distribution in update_states!
end
