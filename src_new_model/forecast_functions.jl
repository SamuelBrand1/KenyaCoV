using NamedArrays

function get_rescaled_contact_matrices()
        M_Data = FileIO.load("data/agemixingmatrix_Kenya_all_types.jld2")
        M_Kenya = M_Data["M_Kenya"]
        M_h = M_Data["M_Kenya_ho"]
        M_o = M_Data["M_Kenya_other"]
        M_w = M_Data["M_Kenya_work"]
        M_s = M_Data["M_Kenya_school"]
        Ronescaler = 1/Real(eigvals(M_Kenya)[end])
        M_h .= M_h.*Ronescaler
        M_o .= M_o.*Ronescaler
        M_w .= M_w.*Ronescaler
        M_s .= M_s.*Ronescaler
        return M_h,M_o,M_s,M_w
end

function get_population_size_matrix()
        df_pop = DataFrame(CSV.File("data/2019_census_age_pyramids_counties.csv"))
        N = NamedArray([Float64(df_pop[i,a+1][1]) for a = 1:17,i = 1:length(df_pop.county)])
        setnames!(N,vcat([string((a-1)*5)*"-"*string(a*5 - 1) for a = 1:16],["80+"]),1)
        setnames!(N,df_pop.county,2)
        setdimnames!(N,["Age","County"])
        return Array(N)
end


function model_ingredients_from_data()
    @load("data/baseline_fitted_parameters.jld2")
    #p_baseline.

    M_h,M_o,M_s,M_w = get_rescaled_contact_matrices()
    N_kenya = permutedims(get_population_size_matrix())

    ## @load agestructuredata_filename N_region_age M_Kenya movements_per_person P_dest ρ T σ rel_detection_rates M_h hosp_rate_by_age ICU_rate_by_age_cond_hosp
    @load "data/data_for_age_structuredmodel_with_counties.jld2"
    n,n_a = size(N_region_age)

    #Population state array
    suspop_kenya = zeros(Int64,n,n_a,n_s) #Array by area, age group and disease state
    for i = 1:n,j=1:n_a
        suspop_kenya[i,j,1] = N_region_age[i,j]
    end

    #Define the change matrix
    dc = sparse(zeros(Int64,n*n_a*n_s,n*n_a*n_ta))
    change_matrix(dc)
    #println("   ***   ",sum(dc))

    #Estimate effective population size in each area after mobility
    #N̂ = T*N_region_age
    #N̂[immobile_age_indices] .= N_region_age[immobile_age_indices]

    #Parameter definition
    P = CoVParameters3(β₀=p_baseline.β₀, α=p_baseline.α, ϵ=p_baseline.ϵ, αₚ=p_baseline.αP,
                        γA=p_baseline.γA, γM=p_baseline.γM, δ=p_baseline.δ, σ=p_baseline.σ, υ=p_baseline.υ, γV=p_baseline.γV,
                        β_home=p_baseline.β_home, β_school=p_baseline.β_school, β_other=p_baseline.β_other, β_work=p_baseline.β_work,
                        M_h=M_h, M_w=M_w, M_s=M_s, M_o=M_o,
                        N=N_kenya,
                        dc = dc)#,  ct=ct)
    return suspop_kenya,P#,P_dest
end


function remake_prob(prob,i,repeat) #Remember to rescale susceptibility by the inverse leading eigenvalue
    _P = deepcopy(prob.p)
    #_P.β = rand(d_R₀)
    return remake(prob,p=_P)
end

function output_simulation_data(sol,i)
    times = sol.prob.tspan[1]:1:sol.prob.tspan[end]

    #=total_cases_A_by_area_and_age=sol(sol.prob.tspan[end])[:,:,9]   ##
    total_cases_M_by_area_and_age=sol(sol.prob.tspan[end])[:,:,10]  ##
    total_cases_V_by_area_and_age=sol(sol.prob.tspan[end])[:,:,11]  ##
    total_q_by_area_and_age=sol(sol.prob.tspan[end])[:,:,16]        ##
    total_qs_by_area_and_age=sol(sol.prob.tspan[end])[:,:,17]       ##
    =#
    incidence_A = diff([sum(sol(t)[:,:,9],dims=2)[:]  for t in times])
    incidence_M = diff([sum(sol(t)[:,:,10],dims=2)[:]  for t in times])
    incidence_V = diff([sum(sol(t)[:,:,11],dims=2)[:]  for t in times])
    incidence_I = incidence_A + incidence_M + incidence_V           ##
    #=incidence_H_by_area = diff([sum(sol(t)[:,:,12],dims=2)[:]  for t in times])
    incidence_H_by_area_and_age = diff([sol(t)[:,:,12]  for t in times])
    T = length(incidence_H_by_area_and_age)
    nc,na = size(incidence_H_by_area_and_age[1])
    total_hosp_occup = zeros(nc,T)
    total_ICU_occup = zeros(nc,T)
    total_new_ICU = zeros(nc,T)
    total_death_incidence = zeros(nc,T)
    hosp_by_area_and_age = zeros(nc,na)
    ICU_by_area_and_age = zeros(nc,na)
    deaths_by_area_and_age = zeros(nc,na)

    for cn in 1:nc, a in 1:na
        hosp_by_area_and_age[cn,a] = sum(sol(sol.prob.tspan[end])[cn,a,12])
    end

    for cn in 1:nc, a in 1:na
        hosp_occup, ICU_occup, new_ICU,death_incidence = generate_hospitalisation_outcomes([inc_h[cn,a] for inc_h in incidence_H_by_area_and_age],a,cn)
        total_hosp_occup[cn,:] .+= hosp_occup
        total_ICU_occup[cn,:] .+= ICU_occup
        total_new_ICU[cn,:] .+= new_ICU
        total_death_incidence[cn,:] .+= death_incidence
        ICU_by_area_and_age[cn,a] = sum(new_ICU)
        deaths_by_area_and_age[cn,a] = sum(death_incidence)
    end
    =#
    return (incidence_A=VectorOfArray(incidence_A)[:,:],
            incidence_M=VectorOfArray(incidence_M)[:,:],
            incidence_V=VectorOfArray(incidence_V)[:,:],
            incidence_I=VectorOfArray(incidence_I)[:,:])
            #=,                        ##
            incidence_H=VectorOfArray(incidence_H_by_area)[:,:],
            hosp_occup_by_area_ts = total_hosp_occup,
            ICU_occup_by_area_ts = total_ICU_occup,
            incidence_ICU_by_area_ts = total_new_ICU,
            death_incidence_by_area_ts = total_death_incidence,
            total_hosp_by_area_and_age = hosp_by_area_and_age,
            total_ICU_by_area_and_age = ICU_by_area_and_age,
            total_deaths_by_area_and_age = deaths_by_area_and_age,
            total_cases_A_by_area_and_age=total_cases_A_by_area_and_age,        ##
            total_cases_M_by_area_and_age=total_cases_M_by_area_and_age,        ##
            total_cases_V_by_area_and_age=total_cases_V_by_area_and_age,        ##
            total_q_by_area_and_age=total_q_by_area_and_age,
            total_qs_by_area_and_age=total_qs_by_area_and_age),false =#           ##
end

function run_sims()
    u0,P#=,P_dest=# = model_ingredients_from_data()
    u0[4,8,5] = 10 # initial Mild symptomatics in Nairobi
    ensemble_prob = EnsembleProblem(DiscreteProblem(nonneg_tauleap,u0,(0.,1*200.),P),
                                    prob_func = remake_prob,
                                    output_func = output_simulation_data)

    return solve(ensemble_prob,EnsembleDistributed(),dt = P.dt, trajectories = 1)
end
