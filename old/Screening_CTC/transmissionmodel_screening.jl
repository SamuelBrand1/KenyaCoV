
function model_ingredients_from_data_screening(agestructuredata_filename,flight_filename,prev_filename)
    @load agestructuredata_filename N_region_age M_Kenya movements_per_person P_dest ρ T σ rel_detection_rates M_Kenya_ho hosp_rate_by_age ICU_rate_by_age_cond_hosp
    n,n_a = size(N_region_age)

    #Population state array
    suspop_kenya = zeros(Int64,n,n_a,n_s) #Array by area, age group and disease state
    for i = 1:n,j=1:n_a
        suspop_kenya[i,j,1] = N_region_age[i,j]
    end

    #Estimate effective population size in each area after mobility
    N̂ = T*N_region_age
    N̂[immobile_age_indices] .= N_region_age[immobile_age_indices]

    # Now get flight numbers and global prevalence
    global_prev = get_prevdata(prev_filename)
    into_mom, into_nai = get_flightdata(flight_filename)

    #Define the change matrix
    dc = sparse(zeros(Int64,n*n_a*n_s,n*n_a*n_ta))
    change_matrix(dc)

    #Parameter definition
    P = CoVParameters_Screening(T = T,ρ = ρ,χ = σ,rel_detection_rate = rel_detection_rates[:,2], #This assumes that the detection rates match ϵ = 0.1
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


function create_KenyaCoV_non_neg_prob(u0,tspan,P::CoVParameters_Screening)
    return prob = DiscreteProblem(nonneg_tauleap,u0,tspan,P)
end
