"""
    function get_flightdata(filename)

Method for collecting flight data from relevant .csv file
"""
function get_flightdata(filename)
    flights=readtable(filename)
    into_mombassa=convert(Vector,flights[:,1])
    into_nairobi=convert(Vector,flights[:,2])
    return into_mombassa,into_nairobi
end
"""
    function get_prevdata(filename)

Method for collecting global prevalence from relevant .csv file
"""
function get_prevdata(filename)
    prev_table=readtable(filename)
    global_prev=convert(Vector,prev_table[:,1])
    return global_prev
end


function transportstructure_params!(P::CoVParameters,ρ::Vector{Float64},transport_matrix)
    d1,d2 = size(transport_matrix)
   #Put in the correct location matrix
   for i = 1:d1,j = 1:d2
       if i != j
           P.T[i,j] = P.ρ[j]*transport_matrix[i,j]
       else
           P.T[i,j] = 1-P.ρ[j]
       end
   end
   return nothing
end

function transportstructure_params!(P::CoVParameters_AS,ρ::Vector{Float64},transport_matrix)
    d1,d2 = size(transport_matrix)
    P.ρ = ρ
   #Put in the correct location matrix
   for i = 1:d1,j = 1:d2
       if i != j
           P.T[i,j] = P.ρ[j]*transport_matrix[i,j]
       else
           P.T[i,j] = 1-P.ρ[j]
       end
   end
   return nothing
end

function transportstructure_params!(P::CoVParameters,ρ::Float64,transport_matrix)
    P.ρ = ρ
    #Put in the correct location matrix
    for i = 1:n,j = 1:n
        if i != j
            P.T[i,j] = P.ρ*transport_matrix[i,j]
        else
            P.T[i,j] = 1-P.ρ
        end
    end
    return nothing
end



function model_ingredients_from_data(datatablename,mixingmatrixname,travelmatrixname,flight_filename,prev_filename)
    KenyaTbl = readtable(datatablename)
    n_data, = size(KenyaTbl)
    @load mixingmatrixname T_opt ρ_county
    @load travelmatrixname P_opt
    if n_data != n
        println("ERROR: difference between number of areas in data, and definition.")
    end
    #Population state array --- read number of areas from T_opt from
    d1,d2 = size(T_opt)
    suspop_kenya = zeros(Int64,d1,n_s,2) #Array by area, state and urban vs rural
    for i = 1:d1
        suspop_kenya[i,1,1] = KenyaTbl.Urban[i]
        if !ismissing(KenyaTbl.Rural[i])
            suspop_kenya[i,1,2] = KenyaTbl.Rural[i]
        end
    end
    N_urb = [sum(suspop_kenya[i,:,1]) for i = 1:d1]
    N_rural = [sum(suspop_kenya[i,:,2]) for i = 1:d1]
    N̂ = T_opt*N_urb + N_rural

    global_prev = get_prevdata(prev_filename)
    into_mom, into_nai = get_flightdata(flight_filename)
    change_matrix(dc,1,1,1,1)
    # Now get flight numbers and global prevalence
    #Parameter definition
    P = CoVParameters(T = T_opt,ρ = ρ_county,
                        into_mom = into_mom, into_nai = into_nai,
                        global_prev = global_prev,
                        N̂=N̂,dc = sparse(dc))
    #Matrix for in-place tau-leap calculations
    return suspop_kenya,P,P_opt

end

function model_ingredients_from_data(agestructuredata_filename,flight_filename,prev_filename)
    @load agestructuredata_filename N_region_age M_Kenya movements_per_person P_dest ρ T σ rel_detection_rates M_Kenya_ho hosp_rate_by_age
    n_wa,n_a = size(N_region_age)
    #Population state array
    suspop_kenya = zeros(Int64,n_wa,n_a,n_s) #Array by area, age group and disease state
    for i = 1:n_wa,j=1:n_a
        suspop_kenya[i,j,1] = N_region_age[i,j]
    end

    #Estimate effective population size in each area after mobility
    # N_mobile = N_region_age[:,mobile_age_indices]
    # N_immobile = N_region_age[:,immobile_age_indices]
    N̂ = T*N_region_age
    N̂[immobile_age_indices] .= N_region_age[immobile_age_indices]

    # Now get flight numbers and global prevalence
    global_prev = get_prevdata(prev_filename)
    into_mom, into_nai = get_flightdata(flight_filename)

    #Define the change matrix
    dc = sparse(zeros(Int64,n_wa*n_a*n_s,n_wa*n_a*n_ta))
    change_matrix(dc)

    #Parameter definition
    P = CoVParameters_AS(T = T,ρ = ρ,χ = σ,rel_detection_rate = rel_detection_rates[:,2], #This assumes that the detection rates match ϵ = 0.1
                        into_mom = into_mom, into_nai = into_nai,
                        global_prev = global_prev,
                        M = M_Kenya,
                        M_ho = M_Kenya_ho,
                        hₐ = hosp_rate_by_age,
                        dc=dc,
                        N̂=N̂)
    #Matrix for in-place tau-leap calculations
    return suspop_kenya,P,P_dest

end




function create_KenyaCoV_prob(u0,tspan,P::CoVParameters)
    u0_vec = u0[:]
    prob_tl = DiscreteProblem(u0_vec,tspan,P)
    return JumpProblem(prob_tl,Direct(),reg_jumps_forKenyaCoV)
end

function solve_KenyaCoV_prob(u0,tspan,P::CoVParameters,dt=1.)
    jump_prob_tl = create_KenyaCoV_prob(u0,tspan,P)
    return solve(jump_prob_tl,SimpleTauLeaping();dt = dt)
end

function create_KenyaCoV_non_neg_prob(u0,tspan,P::CoVParameters)
    u0_vec = u0[:]
    return prob = DiscreteProblem(nonneg_tauleap,u0_vec,tspan,P)
end

function create_KenyaCoV_non_neg_prob(u0,tspan,P::CoVParameters_AS)
    return prob = DiscreteProblem(nonneg_tauleap,u0,tspan,P)
end

function create_KenyaCoV_ode_prob(u0,tspan,P::CoVParameters)
    _u0 = convert.(Float64,u0)
    return prob = ODEProblem(ode_model,_u0,tspan,P)
end

function create_KenyaCoV_ode_prob(u0,tspan,P::CoVParameters_AS)
    _u0 = convert.(Float64,u0)
    return prob = ODEProblem(ode_model,_u0,tspan,P)
end

function calculate_R₀_scale(P::CoVParameters_AS)
    sus_matrix = repeat(P.χ,1,KenyaCoV.n_a)
    trans_matrix = repeat(P.δ*P.rel_detection_rate' + P.ϵ*(1 .- P.δ*P.rel_detection_rate'),KenyaCoV.n_a,1)
    eigs, = eigen(sus_matrix.*M_China.*trans_matrix)
    return Real(eigs[end])
end
function calculate_R₀(P::CoVParameters_AS)
    sus_matrix = repeat(P.χ,1,KenyaCoV.n_a)
    trans_matrix = repeat(P.δ*P.rel_detection_rate' + P.ϵ*(1 .- P.δ*P.rel_detection_rate'),KenyaCoV.n_a,1)
    K = sus_matrix.*P.M.*trans_matrix*P.β/P.γ
    v₀ = zeros(n_a)
    v₀[6] = 1
    eigs,evects = eigen(K)
    return Real(eigs[end]), LinearAlgebra.normalize!(abs.(evects[:,end]),1),(K^5)*v₀
end

function calculate_R₀_homeonly(P::CoVParameters_AS)
    sus_matrix = repeat(P.χ,1,KenyaCoV.n_a)
    trans_matrix = repeat(P.δ*P.rel_detection_rate' + P.ϵ*(1 .- P.δ*P.rel_detection_rate'),KenyaCoV.n_a,1)
    K = sus_matrix.*P.M.*trans_matrix*P.β/P.γ
    v₀ = zeros(n_a)
    v₀[6] = 1
    eigs,evects = eigen(K)
    return Real(eigs[end]), LinearAlgebra.normalize!(abs.(evects[:,end]),1),(K^2)*v₀
end
