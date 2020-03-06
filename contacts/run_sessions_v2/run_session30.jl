push!(LOAD_PATH, "./src")
push!(LOAD_PATH, "./contacts")
using Plots,Parameters,Distributions,DifferentialEquations,JLD2,DataFrames,StatsPlots,FileIO,MAT,RecursiveArrayTools
import KenyaCoV_contacts
using LinearAlgebra:eigen
using Statistics: median, quantile


session_number="30"
results_folder="./contacts/result_sessions_v2/session"*session_number*"/"
mkdir(results_folder)
plt_Nairobi=plot();
solutions=[]
for j=1:100
    global P
    u0,P,P_dest = KenyaCoV_contacts.model_ingredients_from_data("data/data_for_age_structuredmodel.jld2","data/flight_numbers.csv","data/projected_global_prevelance.csv")
    u0[KenyaCoV_contacts.ind_nairobi_as,5,4] = 5#Five initial infecteds in Nairobi in the 20-24 age group

    for wa=1:KenyaCoV_contacts.n_a, a=1:KenyaCoV_contacts.n_a   #**** #Calculate P.Mₚ = Age mixing pobabilities matrix
        P.Mₚ[wa,a]=P.M[wa,a]/sum(P.M[wa,:])
    end

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
    P.τ=1/3.
    P.τₚ=0.8
    P.κ=5
    P.κₘ=7
    P.Δₜ=7
    P.κ_per_event4=30
    P.Κ_max_capacity=1e4

    prob = KenyaCoV_contacts.create_KenyaCoV_non_neg_prob(u0,(0.,365.),P)
    print("simulation "*string(j)*" with proba tau=",P.τₚ,"        ")
    @time sol=solve(prob,FunctionMap(),dt = P.dt)
    push!(solutions,sol)

    times = 0:1:365;    I_area = zeros(Int64,20,length(sol.t))
    for i = 1:20,(ind,t) in enumerate(sol.t)    I_area[i,ind] = sum(sol(t)[i,:,3:4] .+ sol(t)[i,:,10])   end
    plt = plot(sol.t,I_area[1,:], lab = 1);
    for i = 2:20    plot!(plt,sol.t,I_area[i,:],lab = i);     #display(plt);
    savefig(plt,results_folder*"session"*session_number*"_taup"*string(P.τₚ)*"_sim"*string(j)*".png")
    plot!(plt_Nairobi,sol.t,I_area[4,:],lab = "Nairobi"*string(j)); end
end

savefig(plt_Nairobi,results_folder*"detection_tests_"*string(P.τₚ)*"_NAIROBI.png");  #display(plt_Nairobi);
@save results_folder*"solutions.jld2" solutions
