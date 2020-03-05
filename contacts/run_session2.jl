push!(LOAD_PATH, "./src")
push!(LOAD_PATH, "./contacts")
using Plots,Parameters,Distributions,DifferentialEquations,JLD2,DataFrames,StatsPlots,FileIO,MAT,RecursiveArrayTools
import KenyaCoV_contacts
using LinearAlgebra:eigen
using Statistics: median, quantile

#include("../output/forecast_functions.jl")
gr()#Plotting frontend


u0_0,P_0,P_dest_0 = KenyaCoV_contacts.model_ingredients_from_data("data/data_for_age_structuredmodel.jld2","data/flight_numbers.csv","data/projected_global_prevelance.csv")
#u0[KenyaCoV_contacts.ind_nairobi_as,5,4] = 5#Five initial infecteds in Nairobi in the 20-24 age group
P=P_0

plt_Nairobi=plot();
for j=1:10
    #global u0_0,P_0,P_dest_0
    u0,P,P_dest = u0_0,P_0,P_dest_0
    u0[KenyaCoV_contacts.ind_nairobi_as,5,4] = 5#Five initial infecteds in Nairobi in the 20-24 age group

    d_incubation = LogNormal(log(4.8),0.25) #Liu et al
    mean(d_incubation);
    (quantile(d_incubation,0.025),median(d_incubation),quantile(d_incubation,0.975));
    d_R₀ = Gamma(100,2.92/100); ##Liu et al
    mean(d_R₀);
    (quantile(d_R₀,0.025),median(d_R₀),quantile(d_R₀,0.975));

    P.dt = 0.5;
    P.ext_inf_rate = 0.;
    P.ϵ = rand(Uniform(0.,0.5))
    P.δ = 0.2
    P.γ = 1/2.5
    P.σ = 1/rand(d_incubation)
    P.β = rand(d_R₀)*P.γ/(P.δ + P.ϵ*(1-P.δ))

    P.τₚ=0.05
    P.κ=5
    P.κₘ=3
    P.Δₜ=3
    P.κ_per_event4=30
    P.Κ_max_capacity=1e3

    for wa=1:KenyaCoV_contacts.n_a, a=1:KenyaCoV_contacts.n_a   #**** #Calculate P.Mₚ = Age mixing pobabilities matrix
        P.Mₚ[wa,a]=P.M[wa,a]/sum(P.M[wa,:])
    end

    #P.τₚ=0.1
    prob = KenyaCoV_contacts.create_KenyaCoV_non_neg_prob(u0,(0.,365.),P)
    print("simulation with proba tau=",P.τₚ,"           ")
    @time sol=solve(prob,FunctionMap(),dt = P.dt)

    times = 0:1:365;    I_area = zeros(Int64,20,length(sol.t))
    for i = 1:20,(ind,t) in enumerate(sol.t)    I_area[i,ind] = sum(sol(t)[i,:,3:4] .+ sol(t)[i,:,10])   end
    plt = plot(sol.t,I_area[1,:], lab = 1);
    for i = 2:20    plot!(plt,sol.t,I_area[i,:],lab = i);   end
    #display(plt);
    savefig(plt,"./contacts/results_plots_sessions/session2/detection_tests_"*string(P.τₚ)*"_v"*string(j)*".png")

    plot!(plt_Nairobi,sol.t,I_area[4,:],lab = "Nairobi"*string(j));
end
#display(plt_Nairobi);
savefig(plt_Nairobi,"./contacts/results_plots_sessions/session2/detection_tests_"*string(P.τₚ)*"_NAIROBI.png")
