
#####Contact tracing structs
@with_kw mutable struct Contact             #***
    wa::Int=0
    a::Int=0
    s::Int=1                    #state
    contact_t::Float64=0        #contact time
    #contact_dur::Float64=0      #contact duration
end
@with_kw mutable struct IQ_Person      #***
    wa::Int=0
    a::Int=0
    detection_dur::Float64=0.
    contacts::Vector{Contact}=[]
end
@with_kw mutable struct SessionParams
    sc_nb::Int=0
    n_traj::Int=100
    τₚ::Float64=.0;  stop_Q::Bool=true;  Κ_max_capacity12::Int=0;  κ_per_event4::Int=50;  IDs_cfirst::Bool=false;
    dt=.5;  ext_inf_rate::Float64=0.;    ϵ::Float64=.0;  δ::Float64=.0;  γ::Float64=.0; σ::Float64=.0; β::Float64=.0 #β = r_R₀*γ/(δ + ϵ*(1-δ))
    τ::Float64=1/2.;  κ::Int=12;  κₘ::Int=10;   Δₜ::Int=10;
    Κ_max_capacity::Int=0; Κ_max_capacity4::Int=0;
end
######
@with_kw mutable struct CoVParameters_AS
    β::Float64 = 1.
    c_t::Function = t -> 1. #Time varying basic contact rate
    γ::Float64 = 1/2.5
    σ::Float64 = 1/5.
    δ::Float64 = 0.9#Proportion of over 80s who get identified
    τ::Float64 = 0. #current isolation rate for symptomatics
    τ_initial::Float64 = 0. # isolation rate for symptomatics at beginning of epidemic

    clear_quarantine = 1/14. # Two weeks on average to end isolation #**** same as 1/Q_dur

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
    λ::Matrix{Float64} = zeros(n_wa,n_a)#For inplace calculations
    λ_loc::Matrix{Float64} = zeros(n_wa,n_a)#For inplace calculations
    dN::Vector{Int64} = zeros(Int64,n_wa*n_a*n_ta)#For inplace calculations
    poi_rates::Vector{Float64} = zeros(n_wa*n_a*n_ta)#For inplace calculations
    dc::SparseMatrixCSC{Int64,Int64} = sparse(zeros(Int64,n_wa*n_a*n_s,n_wa*n_a*n_ta))#For inplace calculations
    du_linear::Vector{Int64} = zeros(Int64,n_wa*n_a*n_s) #for inplace calculations

    #####Contact tracing parameters
    τₚ::Float64 = τ/(τ+γ)                   #**** probability of detection given that you're diseased (Iᴰ)
    κ::Float64 = 5                          #**** Mean number of contacts per day
    κₘ::Float64 = 7                        #**** nb days IQ remembers his contacts
    Mₚ::Matrix{Float64} = zeros(n_a,n_a)    #**** Age mixing probabilities matrix
    Δₜ::Float64=7                           #describes the tracing period per days (how far back is the detected I asked about his contacts)
    κ_per_event4::Int=30                    #**** Number of contacts traced per event IQ->H
    Κ_max_capacity::Array{Float64,1}=zeros(n_wa)     #**** Tracing capacity, maximum number of traceds                          #!!!! vector for all wa
    Κ_current::Array{Float64,1}=zeros(n_wa)          #**** Total current number of traced contacts (between timestep 0 and t)   #!!!! vector for all wa
    l_IQ::Vector{IQ_Person}=[]  #****  Each IQ is added with a generated period of 1/τ [[wa,a,exponential(1/τ),[contacts]],...], with their contact s@ each timestep
    uₚ::Array{Float64,3}=zeros(n_wa,n_a,n_s) #**** For inplace calculations. Matrix of probabilities: when contacting someone with a specific wa and a, what is the chance of him being S,E,IQ,Iᴰ,Iᴬ,.. WE DO NOT MEET H!
    t_max_capacity::Float64=-1              #**** When the tracing stopped
    stop_Q::Bool=false      #****5 do we stop detecting after we stopped tracing
    IDs_cfirst::Bool=false
end

########
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
