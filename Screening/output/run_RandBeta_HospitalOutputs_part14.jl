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
    #counties = CSV.read("data/2019_census_age_pyramids_counties.csv")
    u0_0,P_0=KenyaCoV_screening.run_scenarios()
    P_0.β=rand(KenyaCoV_screening.d_R₀)
    P_0.screening_delay=1
    P_0.clear_quarantine=14
    P_0.test_sensitivity=.9
    P_0.test_specificity=.98
    P_0.CT_dur=10
    P_0.CT_n=30
    nb_months=1
    S_strategies=[zeros(370)  for i=1:15]
        for i=4:8      for day=(i-1)*30+1:(i-1)*30+30*nb_months     S_strategies[i-3][day]=1e3;      end end
        for i=4:8      for day=(i-1)*30+1:(i-1)*30+30*nb_months     S_strategies[i+5-3][day]=5e3;      end end
        for i=4:8      for day=(i-1)*30+1:(i-1)*30+30*nb_months     S_strategies[i+10-3][day]=1e4;      end end
    println();println();

######## Intervention 6 Social distancing
session=01  ;    scenarios=[1]   ;
    folder="./Screening/output/6_SD/session"*string(session)*"_"*string(n_traj)*"sims/";if !isdir(folder)   mkdir(folder)   end
    @time for i=1:size(scenarios,1)
        println("6_SD\tsession=",session,"\tsc=",scenarios[i],"\tn_traj=",n_traj)
        u0,P=deepcopy(u0_0),deepcopy(P_0)
        P.M=0.8*P.M
        prob = KenyaCoV_screening.create_KenyaCoV_non_neg_prob(u0,(0.,2*365.),P)
        @time results = KenyaCoV_screening.run_βrand(P,prob,n_traj;interventions = CallbackSet())
        @save folder*"SD_sc"*string(scenarios[i])*"_"*string(n_traj)*"sims.jld2" results
    end
println();println();

session=02  ;    scenarios=[1]   ;
    folder="./Screening/output/6_SD/session"*string(session)*"_"*string(n_traj)*"sims/";if !isdir(folder)   mkdir(folder)   end
    @time for i=1:size(scenarios,1)
        println("6_SD\tsession=",session,"\tsc=",scenarios[i],"\tn_traj=",n_traj)
        u0,P=deepcopy(u0_0),deepcopy(P_0)
        P.M=0.6*P.M
        prob = KenyaCoV_screening.create_KenyaCoV_non_neg_prob(u0,(0.,2*365.),P)
        @time results = KenyaCoV_screening.run_βrand(P,prob,n_traj;interventions = CallbackSet())
        @save folder*"SD_sc"*string(scenarios[i])*"_"*string(n_traj)*"sims.jld2" results
    end
println();println();

######## Intervention 7 : Contact tracing of hospitalized only with Social distancing
session=01  ;    scenarios=[1]   ;
    folder="./Screening/output/7_CTHSD/session"*string(session)*"_"*string(n_traj)*"sims/";if !isdir(folder)   mkdir(folder)   end
    @time for i=1:size(scenarios,1)
        println("Intervention=7_CTHSD\tsession=",session,"\tsc=",scenarios[i],"\tn_traj=",n_traj)
        u0,P=deepcopy(u0_0),deepcopy(P_0)
        P.CT_n=30 ########################################################################################
        P.M=0.8*P.M
        cb=CallbackSet(callback_selection_and_mvt_probabilities,callback_detect_H,callback_CT)
        prob = KenyaCoV_screening.create_KenyaCoV_non_neg_prob(u0,(0.,2*365.),P)
        @time results = KenyaCoV_screening.run_βrand(P,prob,n_traj;interventions = cb)
        @save folder*"CTHSD_sc"*string(scenarios[i])*"_"*string(n_traj)*"sims.jld2" results
    end
println();println();

session=02  ;    scenarios=[1]   ;
    folder="./Screening/output/7_CTHSD/session"*string(session)*"_"*string(n_traj)*"sims/";if !isdir(folder)   mkdir(folder)   end
    @time for i=1:size(scenarios,1)
        println("Intervention=7_CTHSD\tsession=",session,"\tsc=",scenarios[i],"\tn_traj=",n_traj)
        u0,P=deepcopy(u0_0),deepcopy(P_0)
        P.CT_n=30 ########################################################################################
        P.M=0.6*P.M
        cb=CallbackSet(callback_selection_and_mvt_probabilities,callback_detect_H,callback_CT)
        prob = KenyaCoV_screening.create_KenyaCoV_non_neg_prob(u0,(0.,2*365.),P)
        @time results = KenyaCoV_screening.run_βrand(P,prob,n_traj;interventions = cb)
        @save folder*"CTHSD_sc"*string(scenarios[i])*"_"*string(n_traj)*"sims.jld2" results
    end
println();println();

session=03  ;    scenarios=[1]   ;
    folder="./Screening/output/7_CTHSD/session"*string(session)*"_"*string(n_traj)*"sims/";if !isdir(folder)   mkdir(folder)   end
    @time for i=1:size(scenarios,1)
        println("Intervention=CTHSD\tsession=",session,"\tsc=",scenarios[i],"\tn_traj=",n_traj)
        u0,P=deepcopy(u0_0),deepcopy(P_0)
        P.CT_n=100 ########################################################################################
        P.M=0.8*P.M
        cb=CallbackSet(callback_selection_and_mvt_probabilities,callback_detect_H,callback_CT)
        prob = KenyaCoV_screening.create_KenyaCoV_non_neg_prob(u0,(0.,2*365.),P)
        @time results = KenyaCoV_screening.run_βrand(P,prob,n_traj;interventions = cb)
        @save folder*"CTHSD_sc"*string(scenarios[i])*"_"*string(n_traj)*"sims.jld2" results
    end
println();println();

session=04  ;    scenarios=[1]   ;
    folder="./Screening/output/7_CTHSD/session"*string(session)*"_"*string(n_traj)*"sims/";if !isdir(folder)   mkdir(folder)   end
    @time for i=1:size(scenarios,1)
        println("Intervention=CTHSD\tsession=",session,"\tsc=",scenarios[i],"\tn_traj=",n_traj)
        u0,P=deepcopy(u0_0),deepcopy(P_0)
        P.CT_n=100 ########################################################################################
        P.M=0.6*P.M
        cb=CallbackSet(callback_selection_and_mvt_probabilities,callback_detect_H,callback_CT)
        prob = KenyaCoV_screening.create_KenyaCoV_non_neg_prob(u0,(0.,2*365.),P)
        @time results = KenyaCoV_screening.run_βrand(P,prob,n_traj;interventions = cb)
        @save folder*"CTHSD_sc"*string(scenarios[i])*"_"*string(n_traj)*"sims.jld2" results
    end
println();println();
