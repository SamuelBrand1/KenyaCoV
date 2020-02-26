
function transportstructure_params!(P::CoVParameters,ρ::Vector{Float64},transport_matrix)
   #Put in the correct location matrix
   for i = 1:n,j = 1:n
       if i != j
           P.T[i,j] = P.ρ[j]*transport_matrix[i,j]
           P.T[i,j] = P.ρ*transport_matrix[i,j]
       else
           P.T[i,j] = 1-P.ρ[j]
           P.T[i,j] = 1-P.ρ
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
    #Population state array
    suspop_kenya = zeros(Int64,n,n_s,2) #Array by area, state and urban vs rural
    for i = 1:n
        suspop_kenya[i,1,1] = KenyaTbl.Urban[i]
        if !ismissing(KenyaTbl.Rural[i])
            suspop_kenya[i,1,2] = KenyaTbl.Rural[i]
        end
    end
    N_urb = [sum(suspop_kenya[i,:,1]) for i = 1:n]
    N_rural = [sum(suspop_kenya[i,:,2]) for i = 1:n]
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
    @load agestructuredata_filename N_region_age agemixingmatrix movements_per_person P_dest ρ T age_specific_sus

    #Population state array
    suspop_kenya = zeros(Int64,n_wa,n_a,n_s) #Array by area, age group and disease state
    for i = 1:n_wa,j=1:n_a
        suspop_kenya[i,j,1] = N_region_age[i,j]
    end

    #Estimate effective population size in each area after
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
    P = CoVParameters_AS(T = T,ρ = ρ,χ = age_specific_sus,
                        into_mom = into_mom, into_nai = into_nai,
                        global_prev = global_prev,
                        M = agemixingmatrix,
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
