
push!(LOAD_PATH, "./src")
    push!(LOAD_PATH, "./data")
    push!(LOAD_PATH, "./Vaccination")

    using DifferentialEquations, Plots,CSV,JLD2,SparseArrays,Revise, Parameters
    using LinearAlgebra
    import KenyaCoV_vaccination

# Load data for model --- and create an initially susceptible population
u0,P,transport_matrix = KenyaCoV_vaccination.model_ingredients_from_data_vacc();
    sol_ode = solve(ODEProblem(KenyaCoV_vaccination.ode_model,convert.(Float64,u0),(0.,365.),P),Tsit5());

#Plot
I_ode = [sum(sol_ode(t)[:,:,4:6]) for t in 0.:1.:sol_ode.t[end]]
    plt = plot(0.:1.:sol_ode.t[end],I_ode .+ 1,color = :red, ylabel = "Infecteds",legend=false)


## Solve vaccination ODE
v_mode=1
u0,P,transport_matrix = KenyaCoV_vaccination.model_ingredients_from_data_vacc();
if v_mode==1
    P.v=1/300.   #vaccination rate
    cb=Callba
elseif v_mode==2
    P.vₐ .= 1/300.
else
    ##Set vaccination parameters depending on time
    timing_daily(u,t,integrator) = t>=1 && t==floor(t)
    function affect_set_daily_vaccination_strategy!(integrator)
        @unpack vᵗ,vₐ = integrator.p
        vₐ[:,:] .= vᵗ[Int(integrator.t),:,:]
    end
    callback_set_daily_vaccination_strategy = DiscreteCallback(timing_daily,affect_set_daily_vaccination_strategy!)

end


    ##Solve vacc using the vaccination rate v
    u0,P,transport_matrix = KenyaCoV_vaccination.model_ingredients_from_data_vacc();
    P.v_mode = 1
    P.vᵗ[80:100,:,11:end] .= 10000000 #Vaccinating sonly the over-50s in all Kenya
    sol_ode_v = solve(ODEProblem(KenyaCoV_vaccination.ode_model_vacc,convert.(Float64,u0),(0.,365.),P),Tsit5(),callback=callback_set_daily_vaccination_strategy,tstops=[1:1:365;]);




    P.vᵗ[80:100,:,11:end] .= 10000000 #Vaccinating sonly the over-50s in all Kenya
    sol_ode_v = solve(ODEProblem(KenyaCoV_vaccination.ode_model_vacc,convert.(Float64,u0),(0.,365.),P),Tsit5(),callback=callback_set_daily_vaccination_strategy,tstops=[1:1:365;]);

    plt = plot(0.:1.:sol_ode.t[end],I_ode .+ 1,color = :red, ylabel = "Infecteds",legend=false)
    I_ode_v = [sum(sol_ode_v(t)[:,:,4:6]) for t in 0.:1.:sol_ode_v.t[end]]
    plot!(plt,0.:1.:sol_ode_v.t[end],I_ode_v .+ 1,color=:blue)#, ylabel = "Infecteds",legend=false)
