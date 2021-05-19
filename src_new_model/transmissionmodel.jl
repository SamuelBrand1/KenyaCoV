"""
    function get_flightdata(filename)

Method for collecting flight data from relevant .csv file
"""
function get_flightdata(filename)
    flights=CSV.read(filename)
    into_mombassa=convert(Vector,flights[:,1])
    into_nairobi=convert(Vector,flights[:,2])
    return into_mombassa,into_nairobi
end
"""
    function get_prevdata(filename)

Method for collecting global prevalence from relevant .csv file
"""
function get_prevdata(filename)
    prev_table=CSV.read(filename)
    global_prev=convert(Vector,prev_table[:,1])
    return global_prev
end


# function transportstructure_params!(P::CoVParameters,ρ::Vector{Float64},transport_matrix)
#     d1,d2 = size(transport_matrix)
#    #Put in the correct location matrix
#    for i = 1:d1,j = 1:d2
#        if i != j
#            P.T[i,j] = P.ρ[j]*transport_matrix[i,j]
#        else
#            P.T[i,j] = 1-P.ρ[j]
#        end
#    end
#    return nothing
# end
#
# function transportstructure_params!(P::CoVParameters_AS,ρ::Vector{Float64},transport_matrix)
#     d1,d2 = size(transport_matrix)
#     P.ρ = ρ
#    #Put in the correct location matrix
#    for i = 1:d1,j = 1:d2
#        if i != j
#            P.T[i,j] = P.ρ[j]*transport_matrix[i,j]
#        else
#            P.T[i,j] = 1-P.ρ[j]
#        end
#    end
#    return nothing
# end
#
# function transportstructure_params!(P::CoVParameters,ρ::Float64,transport_matrix)
#     P.ρ = ρ
#     #Put in the correct location matrix
#     for i = 1:n,j = 1:n
#         if i != j
#             P.T[i,j] = P.ρ*transport_matrix[i,j]
#         else
#             P.T[i,j] = 1-P.ρ
#         end
#     end
#     return nothing
# end


"""
    function model_ingredients_from_data(agestructuredata_filename,flight_filename,prev_filename)

model_ingredients_from_data reads in a JLD2 file which includes all the relevant information required for simulation.
"""
function model_ingredients_from_data(agestructuredata_filename,flight_filename,prev_filename)
    @load agestructuredata_filename N_region_age M_Kenya movements_per_person P_dest ρ T σ rel_detection_rates M_Kenya_ho hosp_rate_by_age ICU_rate_by_age_cond_hosp
    n_wa,n_a = size(N_region_age)

    #Population state array
    suspop_kenya = zeros(Int64,n_wa,n_a,n_s) #Array by area, age group and disease state
    for i = 1:n_wa,j=1:n_a
        suspop_kenya[i,j,1] = N_region_age[i,j]
    end

    #Estimate effective population size in each area after mobility
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
                        ICUₐ = ICU_rate_by_age_cond_hosp,
                        dc=dc,
                        N̂=N̂)
    #Matrix for in-place tau-leap calculations
    return suspop_kenya,P,P_dest

end

"""
    function create_KenyaCoV_non_neg_prob(u0,tspan,P::CoVParameters_AS)

This generates a DEProblem which can be solved using the DifferentialEquations package FunctionMap solver.
"""
function create_KenyaCoV_non_neg_prob(u0,tspan,P::CoVParameters_AS)
    return prob = DiscreteProblem(nonneg_tauleap,u0,tspan,P)
end

"""
    function create_KenyaCoV_ode_prob(u0,tspan,P::CoVParameters_AS)

This generates an ODEProblem which can be solved using a DifferentialEquations package ODE solver.
"""
function create_KenyaCoV_ode_prob(u0,tspan,P::CoVParameters_AS)
    _u0 = convert.(Float64,u0)
    return prob = ODEProblem(ode_model,_u0,tspan,P)
end

# function calculate_R₀_scale(P::CoVParameters_AS)
#     sus_matrix = repeat(P.χ,1,KenyaCoV.n_a)
#     trans_matrix = repeat(P.δ*P.rel_detection_rate' + P.ϵ*(1 .- P.δ*P.rel_detection_rate'),KenyaCoV.n_a,1)
#     eigs, = eigen(sus_matrix.*M_China.*trans_matrix)
#     return Real(eigs[end])
# end

"""
    function calculate_R₀(P::CoVParameters_AS)

This outputs the R₀ and corresponding eigenvector associated with the current parameters.
"""
function calculate_R₀(P::CoVParameters_AS)
    sus_matrix = repeat(P.χ,1,KenyaCoV.n_a)
    trans_matrix = repeat(P.δ*P.rel_detection_rate' + P.ϵ*(1 .- P.δ*P.rel_detection_rate'),KenyaCoV.n_a,1)
    K = sus_matrix.*P.M.*trans_matrix*P.β/P.γ
    eigs,evects = eigen(K)
    return Real(eigs[end]), LinearAlgebra.normalize!(abs.(evects[:,end]),1)
end

# function calculate_R₀_homeonly(P::CoVParameters_AS)
#     sus_matrix = repeat(P.χ,1,KenyaCoV.n_a)
#     trans_matrix = repeat(P.δ*P.rel_detection_rate' + P.ϵ*(1 .- P.δ*P.rel_detection_rate'),KenyaCoV.n_a,1)
#     K = sus_matrix.*P.M.*trans_matrix*P.β/P.γ
#     v₀ = zeros(n_a)
#     v₀[6] = 1
#     eigs,evects = eigen(K)
#     return Real(eigs[end]), LinearAlgebra.normalize!(abs.(evects[:,end]),1),(K^2)*v₀
# end
