
"""
    function model_ingredients_from_data(agestructuredata_filename,flight_filename,prev_filename)

model_ingredients_from_data reads in a JLD2 file which includes all the relevant information required for simulation.
"""
function model_ingredients_from_data_vacc()
    agestructuredata_filename="./data/data_for_age_structuredmodel_with_counties.jld2"
    flight_filename="./data/flight_numbers.csv"
    prev_filename="./data/projected_global_prevelance.csv"
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
    P = CoVParameters_Vaccination(T = T,ρ = ρ,χ = σ,rel_detection_rate = rel_detection_rates[:,2], #This assumes that the detection rates match ϵ = 0.1
                        into_mom = into_mom, into_nai = into_nai,
                        global_prev = global_prev,
                        M = M_Kenya,
                        M_ho = M_Kenya_ho,
                        hₐ = hosp_rate_by_age,
                        ICUₐ = ICU_rate_by_age_cond_hosp,
                        dc=dc,
                        N̂=N̂)



    u0,P,transport_matrix = suspop_kenya,P,P_dest
        # #Change some parameters away from default
        for (i,p) in enumerate(P.global_prev)
            P.global_prev[i] = 0.
        end
        P.ϵ = 1.
        P.χ = ones(17)

        counties = CSV.read("./data/2019_census_age_pyramids_counties.csv")
        Nairobi_index = findfirst(counties.county .== "Nairobi")
        Mombassa_index = findfirst(counties.county .== "Mombasa")
        Kwale_index = findfirst(counties.county .== "Kwale")
        Kilifi_index = findfirst(counties.county .== "Kilifi")
        Mandera_index  = findfirst(counties.county .== "Mandera")

        #Rescale by maximum eigenvalue of the age mixing matrix rescaled by susceptibility
        sus_matrix = repeat(P.χ,1,17)
        R_A = P.ϵ*((1/P.σ₂) + (1/P.γ) )
        R_M = (P.ϵ/P.σ₂) + (P.ϵ_D/P.γ)
        R_V = (P.ϵ/P.σ₂) + (P.ϵ_V/P.τ)
        R_vector = [(1-P.rel_detection_rate[a])*R_A + P.rel_detection_rate[a]*(1-P.hₐ[a])*R_M + P.rel_detection_rate[a]*P.hₐ[a]*R_V for a = 1:17]
        inf_matrix = repeat(R_vector',17,1)

        eigs, = eigen(sus_matrix.*P.M.*inf_matrix)
        max_eigval = Real(eigs[end])


        P.β = 2.5/max_eigval# Choosing the transmission rate that matches R₀ = 2.5
        P.dt = 0.25; #<---- KenyaCoV is a discrete time simulation, timestep here is 1/10th of a day

        u0[Nairobi_index,8,3] = 1000 #10 initial pre-symptomatics in Nairobi
        u0[Mombassa_index,8,3] = 5000 #10 initial pre-symptomatics in Mombasa
        u0[Mandera_index,8,3] = 2000 #5 initial pre-symptomatics in Mandera
    #Matrix for in-place tau-leap calculations
    return u0,P,transport_matrix

end
