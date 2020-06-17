#cd("/home/rabiaaziza/julia_projects/2020-06-10_Screening_cleaned/")
    push!(LOAD_PATH, "./src")
    push!(LOAD_PATH, "./Screening")
    using Plots,Parameters,Distributions,DifferentialEquations,JLD2,DataFrames,StatsPlots,FileIO,MAT,RecursiveArrayTools,CSV
    using Revise
    import KenyaCoV_screening
    using LinearAlgebra:eigen
    using Statistics: median, quantile
    include("../callbacks_screening.jl")

########
n_traj=500
    counties = CSV.read("data/2019_census_age_pyramids_counties.csv")
    u0_0,P_0=KenyaCoV_screening.run_scenarios()
    P_0.β=rand(KenyaCoV_screening.d_R₀)
    P_0.screening_delay=1
    P_0.clear_quarantine=14
    P_0.test_sensitivity=.9
    P_0.test_specificity=.98
    P_0.CT_dur=10
    P_0.CT_n=30
    nb_months=1
    S_strategies=[zeros(370)  for i=1:10]
        for i=3:7      for day=(i-1)*30+1:(i-1)*30+30*nb_months     S_strategies[i-2][day]=1e3;      end end
        for i=3:7      for day=(i-1)*30+1:(i-1)*30+30*nb_months     S_strategies[i+5-2][day]=5e3;    end end
    println();println();

######## Intervention 5 : Mass screening with contact tracing
session=02  ;    scenarios=[i  for i=1:10]
    folder="./Screening/output/5_MSCT/session"*string(session)*"_"*string(n_traj)*"sims/";if !isdir(folder)   mkdir(folder)   end
    @time for i=1:size(scenarios,1)
        println("Intervention=MSCT\tsession=",session,"\tsc=",scenarios[i],"\tn_traj=",n_traj)
        u0,P=deepcopy(u0_0),deepcopy(P_0)
        for r=1:KenyaCoV_screening.n_wa P.S_strategy[r,:]=S_strategies[i];  end # same in all kenya
        #P.S_strategy[30,:]=S_strategies[i];P.S_strategy[28,:]=S_strategies[i];# in Nairobi and Monbassa only
        cb=CallbackSet(callback_selection_and_mvt_probabilities,callback_MS,callback_CT)
        prob = KenyaCoV_screening.create_KenyaCoV_non_neg_prob(u0,(0.,1*365.),P)
        @time results = KenyaCoV_screening.run_βrand(P,prob,n_traj;interventions = cb)
        @save folder*"MSCT_sc"*string(scenarios[i])*"_"*string(n_traj)*"sims.jld2" results
    end
