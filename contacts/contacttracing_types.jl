
@with_kw mutable struct Contact             #***
    wa::Int=0
    a::Int=0
    s::Int=1                    #state
    contact_t::Float64=0        #contact time
    contact_dur::Float64=0      #contact duration
end
@with_kw mutable struct IQ_Person      #***
    wa::Int=0
    a::Int=0
    #detection_dur::Float64=0.
    contacts::Vector{Contact}=[]
end
#=@with_kw mutable struct Q_Person      #****5
    wa::Int=0
    a::Int=0
    s::Int=0
    q_dur::Vector{Float64}=[]
end=#
@with_kw mutable struct SessionParams
    sc_nb::Int=0
    n_traj::Int=100
    R₀::Float64=.0;τₚ::Vector{Float64}=zeros(n_a);  #=stop_Q::Bool=true;  Κ_max_capacity12::Int=0;  =#κ_per_event4::Int=50;  IDs_cfirst::Bool=false;
    dt=.5;  ext_inf_rate::Float64=0.;    ϵ::Float64=.0;  δ::Float64=.0;  γ::Float64=.0; σ::Float64=.0; β::Float64=.0 #β = r_R₀*γ/(δ + ϵ*(1-δ))
    τ::Float64=1/2.;  κ::Int=12;  κₘ::Int=10;   Δₜ::Int=10;
    #Κ_max_capacity::Int=0; Κ_max_capacity4::Int=0;
    cumI::Float64=-1;cumI_diff::Float64=-1;#for plotting
    CT_Imin::Array{Float64,1}=zeros(n_wa);CT_dur::Array{Float64,1}=zeros(n_wa)
end

@with_kw mutable struct CoVParameters_AS
    β::Float64 = 2.5/3.6
    γ::Float64 = 1/3.6
    σ::Float64 = 1/2.
    δ::Float64 = 0.05#Proportion of symptomatic/diseased vs non-symptomatic cases
    τ::Float64 = 1/15. #treatment/isolation rate for symptomatics
    μ₁::Float64 = 0.0#Excess mortality due to disease
    ϵ::Float64 = 0.1 #Relative infectiousness of undetectable infecteds
    χ::Vector{Float64} = ones(n_a)
    ρ::Vector{Float64} = [0.01 for i in 1:n] #Time spent outside of area
    T::Matrix{Float64} = zeros(n,n)#Probability distributions of location for mobile individuals
    M::Matrix{Float64} = zeros(n_a,n_a) #Age mixing matrix
    ext_inf_rate::Float64 = 0. #Scales the contact rate with the infecteds arriving via air
    into_mom::Vector{Int}#Number of people flying into Mombassa each day
    into_nai::Vector{Int}#Number of people flying into Nairobi each day
    global_prev::Vector{Float64}
    dt::Float64 = 1. #Useful for the non-negative method
    Î::Matrix{Float64} = zeros(n_wa,n_a) #For inplace calculations
    N̂::Matrix{Float64} = zeros(n_wa,n_a)#For inplace calculations
    λ::Matrix{Float64} = zeros(n_wa,n_a)#For inplace calculations
    λ_loc::Matrix{Float64} = zeros(n_wa,n_a)#For inplace calculations
    dN::Vector{Int64} = zeros(Int64,n_wa*n_a*n_ta)#For inplace calculations
    poi_rates::Vector{Float64} = zeros(n_wa*n_a*n_ta)#For inplace calculations
    dc::SparseMatrixCSC{Int64,Int64} = sparse(zeros(Int64,n_wa*n_a*n_s,n_wa*n_a*n_ta))#For inplace calculations
    du_linear::Vector{Int64} = zeros(Int64,n_wa*n_a*n_s) #for inplace calculations

    τₚ::Vector{Float64} = zeros(n_a)        #****7 detection probability per area           #**** probability of detection given that you're diseased (Iᴰ)
    κ::Float64 = 5                          #**** Mean number of contacts per day
    κₘ::Float64 = 7                        #**** nb days IQ remembers his contacts
    Mₚ::Matrix{Float64} = zeros(n_a,n_a)    #**** Age mixing probabilities matrix
    #contacts::Array{Tuple{Int64,Int64,Float64},1}=[] #Array of tuples (wa_infector,wa_infected,time)
    Δₜ::Float64=7                           #describes the tracing period per days (how far back is the detected I asked about his contacts)
    κ_per_event4::Int=30                    #**** Number of contacts traced per event IQ->H
    #****7 Κ_max_capacity::Array{Float64,1}=zeros(n_wa)     #**** Tracing capacity, maximum number of traceds                          #!!!! vector for all wa
    #****7 Κ_current::Array{Float64,1}=zeros(n_wa)          #**** Total current number of traced contacts (between timestep 0 and t)   #!!!! vector for all wa

    l_IQ::Vector{IQ_Person}=[]  #****  Each IQ is added with a generated period of 1/τ [[wa,a,exponential(1/τ),[contacts]],...], with their contact s@ each timestep
    uₚ::Array{Float64,3}=zeros(n_wa,n_a,n_s) #**** For inplace calculations. Matrix of probabilities: when contacting someone with a specific wa and a, what is the chance of him being S,E,IQ,Iᴰ,Iᴬ,.. WE DO NOT MEET H!
    #****7 t_max_capacity::Float64=-1              #**** When the tracing stopped

    Q_dur::Float64=14       #****4 Quarantine duration in time
    #l_Q=[Q_Person(wa,a,s,[]) for wa=1:n_wa,a=1:n_a,s=1:n_s]   #****5
    #stop_Q::Bool=false      #****5 do we stop detecting after we stopped tracing
    IDs_cfirst::Bool=false

    CT_Imin::Array{Float64,1}=zeros(n_wa) #****7 minimum number of detecteds (ID+IQ) PER REGION after which we start contact tracing
    CT_dur::Array{Float64,1}=zeros(n_wa)     #****7 duration in days of contact tracing after it started when we got n_min_startCT
    t_startedCT::Array{Float64,1}=[-1  for i=1:n_wa] #****7 the time CT started per region
end
#=
States:
1 -> S
2 -> E
3 -> I_subclinical
4 -> I_diseased                                         #****  Iᴰ  Not to be hospitalised
5 -> Q(uarantined)                                      #***! We renamed H to Q (in ALL CURRENT FILE)
6 -> Recovered
7 -> Cumulative I_sub
8 -> Cumulative I_dis
9 -> Cumulative I_Q                                    #****2  Cumulative dead BECOMES Cumulative I_Q
10 -> I_Q diseased to be quarantined                     #****  IQ
11 -> Q_S Susceptibles in quarantine                    #****3      (WAS: 11 -> C(daily contacteds)     #****2)
12 -> Cumulative infected contacteds                     #****2 #****6 for infecteds only E IA ID IQ
13 -> Cumulative contacteds (all)
14 -> Cumulative deaths

Events for each wider area and age group:

1-> Transmission
2-> Incubation into asymptomatic E->A
3-> Incubation into diseased E->D                       #****
4-> Diseased become hospitalised/treated                #****  Iᴰ becomes hospitalised : was Iᴰ->Q BECOMES IQ->Q
5-> Quarantined recover Q->R
6-> Diseased recover D->R
7-> Asymptomatics recover A->R
8-> Quarantined->death
9-> Incubation into diseased to be detected  E->IQ    #****
10->Diseased to be Quarantined who recover DQ->R       #**** IQ->R
11->S to Q Contacts
12->E to Q
13->I_A to Q
14->I_D to Q
15->I_Q to Q
16->Qs->S                     #****4
=#

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
