using Distributed
@everywhere push!(LOAD_PATH, "./src")
    @everywhere push!(LOAD_PATH, "./Screening")
    @everywhere using Parameters,Distributions,DifferentialEquations,JLD2,DataFrames,FileIO,MAT,RecursiveArrayTools,CSV,Revise
@everywhere begin
        import KenyaCoV_screening
        using LinearAlgebra:eigen
        using Statistics: median, quantile
        include("./Screening/callbacks_screening.jl")
###### PARAMS TO MODIFY
    code_folder="/home/gemvi/rabia_aziza_midasnetwork_us/julia_projects/2020-07-22_KenyaCoV/Screening/simFiles_4/rand_K_1000tests/"
        code_folder="./Screening/simFiles_4/rand_K_1000tests/"
    output_folder="./Screening/output_rand_K_1000tests/"
    n_traj=1000
    ϵ=1.
    strategy_nb_tests=1000
    #nb_months=1
    is_distributed=true
    R0_is_random=true
    scenarios=[i for i=1:5]

######
    u0_0,P_0=KenyaCoV_screening.baseline_config(ϵ)
    if R0_is_random
        P_0.β = rand(KenyaCoV_screening.d_R₀) #Choose R₀ randomly from 2-3 range
    else
        P_0.β = 2
    end

    P_0.screening_delay=1;    P_0.clear_quarantine=14
    P_0.test_sensitivity=.9;    P_0.test_specificity=.98
    P_0.CT_dur=14;    P_0.CT_n=30

    function get_strategy(nb_months,strategy_nb_tests)
        strategies=[zeros(Int64,370)  for i=1:8];
        for i=1:8
            for day=(i-1)*30+1:(i-1)*30+30*nb_months
                strategies[i][day]=strategy_nb_tests;
            end
        end
        return strategies
    end


    function CT_nₜ_fixed()
        return [30 for t=1:365]
    end
        function CT_nₜ_decaying_1month(sc)#CT_n full force CT (fixed to 30) for 1 month starting at beginning of sc_th month, then decaying over 3 months
            daily_dec=1 #daily decrease of the capacity of tracing
            V=[30 for t=1:sc*30]
            for t=sc*30+1:3:sc*30+3*30
                push!(V,V[end]);push!(V,V[end]);push!(V,V[end]-1)
            end
            return V
        end
        function CT_nₜ_decaying_3months(sc)#CT_n full force CT (fixed to 30) for 3 months starting at beginning of sc_th month, then decaying over 3 months
            #daily_dec=1 #daily decrease of the capacity of tracing
            V=[30 for t=1:(sc+2)*30]
            for t=(sc+2)*30+1:3:(sc+2)*30+3*30
                push!(V,V[end]);push!(V,V[end]);push!(V,V[end]-1)
            end
            return V
        end
        #=v=CT_nₜ_decaying_2(3)
            plot(v,xticks=[i for i=0:30:270])=#

    if !isdir(output_folder)   mkdir(output_folder)   end
        if !isdir(output_folder*"I0_NoInterv/")   mkdir(output_folder*"I0_NoInterv/")   end
        if !isdir(output_folder*"I1_CTH/")   mkdir(output_folder*"I1_CTH/")   end
        if !isdir(output_folder*"I2_SympS/")   mkdir(output_folder*"I2_SympS/")   end
        if !isdir(output_folder*"I3_SympSCT/")   mkdir(output_folder*"I3_SympSCT/")   end
        if !isdir(output_folder*"I4_MS/")   mkdir(output_folder*"I4_MS/")   end
        if !isdir(output_folder*"I5_MSCT/")   mkdir(output_folder*"I5_MSCT/")   end



    function run_intervention_session(session,scenarios,intervention_label,nb_months,cb,n_traj,CT_n_decaying=0)
        folder=output_folder*intervention_label*"/session"*string(session)*"_"*string(n_traj)*"sims/";if !isdir(folder)   mkdir(folder)   end
        @time for i ∈ scenarios
            println("Intervention="*intervention_label*"\tsession=",session,"\tsc=",i,"\tn_traj=",n_traj)
            u0,P=deepcopy(u0_0),deepcopy(P_0)
            if CT_n_decaying==0
                P.CT_nₜ=CT_nₜ_fixed()
            elseif CT_n_decaying==1
                P.CT_nₜ=CT_nₜ_decaying_1month(i)
            elseif CT_n_decaying==3
                P.CT_nₜ=CT_nₜ_decaying_3months(i)
            end
    		strategy=get_strategy(nb_months,strategy_nb_tests)
            for r=1:KenyaCoV_screening.n,day=1:370     #Same strategy in all counties
            #for r∈[28,30],day=1:370     #Strategy in Nairobi and Mombasa
                P.strategy[r,day]=strategy[i][day]
            end
            prob = KenyaCoV_screening.create_KenyaCoV_non_neg_prob(u0,(0.,1*365.),P)
            @time results = KenyaCoV_screening.run_simulations(P,prob,n_traj,R0_is_random,is_distributed;interventions=cb)
            @save folder*"sc"*string(i)*".jld2" results
        end
        println();println();
    end

end
