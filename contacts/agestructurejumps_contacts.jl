"""
Basic representation of the state of the age structured model:
u[wa_index,age_group,disease state]

This is row major unpacked so
first 1,...,n_wa entries are 0-4 year old susceptibles in wider areas 1,...,n_wa
then n_wa+1,...,2n_wa are 5-9 year old susceptibles in wider areas 1,...,n_wa
.
.
.
then n_wa*n_s + 1, ...,  n_wa*n_s + n_wa entries are 0-4 year old exposed in wider areas 1,...,n_wa
.
.
.

States:
1 -> S
2 -> E
3 -> I_subclinical
4 -> I_diseased                                         #****  Iᴰ  Not to be hospitalised
5 -> Q(uarantined)                                      #***! We renamed H to Q (in ALL CURRENT FILE)
6 -> Recovered
7 -> Cumulative I_sub
8 -> Cumulative I_dis
9 -> Cumulative Dead
10 -> I_d(iseased_to_be_)h(ospitalised)                 #****  IQ

Events for each wider area and age group:

1-> Transmission
2-> Incubation into asymptomatic E->A
3-> Incubation into diseased E->D                       #****
4-> Diseased become hospitalised/treated                #****  Iᴰ becomes hospitalised : was Iᴰ->Q BECOMES IQ->Q
5-> Quarantined recover Q->R
6-> Diseased recover D->R
7-> Asymptomatics recover A->R
8-> Quarantined->death
9-> Incubation into diseased to be detected  E->IQ    #****
10->Diseased to be Quarantined who recover DQ->R       #**** IQ->R

"""

dc_age = zeros(Int64, n_wa * n_a * n_s, n_ta * n * n_a)

function calculate_infection_rates!(u, p::CoVParameters_AS, t)
    Iᴬ = @view u[:, :, 3]
    Iᴰ = @view(u[:, :, 4]) .+ @view(u[:, :, 10])                                                 #**** also summing for IQ
    mul!(p.Î, p.T, p.ϵ * Iᴬ .+ Iᴰ)#Local infecteds **if** everyone moved around
    p.Î[:, immobile_age_indices] .= p.ϵ * Iᴬ[:, immobile_age_indices] .+ Iᴰ[:, immobile_age_indices]#This corrects for immobility
    mul!(p.λ_loc, p.β .* (p.Î ./ p.N̂), p.M')#Local force of infection due to age-mixing
    mul!(p.λ, p.T', p.λ_loc)#this accounts for mobile susceptibles contracting away from home
    p.λ[:, immobile_age_indices] .= p.λ_loc[:, immobile_age_indices]#This corrects for immobility
    p.λ[ind_mombasa_as, :] .+= p.ext_inf_rate * import_rate_mom(t, p.into_mom, p.global_prev)
    p.λ[ind_nairobi_as, :] .+= p.ext_inf_rate * import_rate_nai(t, p.into_nai, p.global_prev)
    return nothing
end


function rates(out, u, p::CoVParameters_AS, t)
    @unpack λ, γ, σ, δ, τ, μ₁, χ, τₚ = p                  #**** added τₚ
    #τₚ=τ/(τ+γ)                                      #****
    calculate_infection_rates!(u, p, t)
    for k = 1:length(out)
        i, a, eventtype = Tuple(index_as_events[k])
        if eventtype == 1
            out[k] = χ[a] * λ[i, a] * u[i, a, 1] #Transmission
        end
        if eventtype == 2
            out[k] = (1 - δ) * σ * u[i, a, 2] # E->A
        end
        if eventtype == 3
            out[k] = δ * σ * (1 - τₚ) * u[i, a, 2] # E->D         #**** added τₚ
        end
        if eventtype == 4
            #out[k] = τ*u[i,a,4] # D->Q
            out[k] = 0#τ*u[i,a,4]                   #**** was D->Q became IQ->Q, rate =0 because we do this manually in the function update_l_IQ()
        end
        if eventtype == 5
            out[k] = γ * u[i, a, 5] # Q->R
        end
        if eventtype == 6
            out[k] = γ * u[i, a, 4] # D->R
        end
        if eventtype == 7
            out[k] = γ * u[i, a, 3] # A->R
        end
        if eventtype == 8
            out[k] = μ₁ * u[i, a, 5] # Q->death
        end
        if eventtype == 9                            #**** adding the event for E -> IQ
            out[k] = δ * σ * τₚ * u[i, a, 2]
        end
        if eventtype == 10                           #**** adding the event IQ->R
            out[k] = γ * u[i, a, 10] # IQ->R
        end
    end
end

function change_matrix(dc)
    d1, d2 = size(dc)
    for k = 1:d2
        i, a, eventtype = Tuple(index_as_events[k])
        if eventtype == 1 #Transmission
            ind_S = linear_as[i, a, 1]
            ind_E = linear_as[i, a, 2]
            dc[ind_S, k] = -1
            dc[ind_E, k] = 1
        end
        if eventtype == 2 # E->A
            ind_E = linear_as[i, a, 2]
            ind_A = linear_as[i, a, 3]
            ind_cumA = linear_as[i, a, 7]
            dc[ind_E, k] = -1
            dc[ind_A, k] = 1
            dc[ind_cumA, k] = 1
        end
        if eventtype == 3 # E->D
            ind_E = linear_as[i, a, 2]
            ind_D = linear_as[i, a, 4]
            ind_cumD = linear_as[i, a, 8]
            dc[ind_E, k] = -1
            dc[ind_D, k] = 1
            dc[ind_cumD, k] = 1
        end
        if eventtype == 4                              #**** was Iᴰ->Q   BECOMES IQ->Q
            ind_DQ = linear_as[i, a, 10]
            ind_Q = linear_as[i, a, 5]
            dc[ind_DQ, k] = -1
            dc[ind_Q, k] = 1
        end
        if eventtype == 5 # Q->R
            ind_Q = linear_as[i, a, 5]
            ind_R = linear_as[i, a, 6]
            dc[ind_Q, k] = -1
            dc[ind_R, k] = 1
        end
        if eventtype == 6 # D->R
            ind_D = linear_as[i, a, 4]
            ind_R = linear_as[i, a, 6]
            dc[ind_D, k] = -1
            dc[ind_R, k] = 1
        end
        if eventtype == 7 # A->R
            ind_A = linear_as[i, a, 3]
            ind_R = linear_as[i, a, 6]
            dc[ind_A, k] = -1
            dc[ind_R, k] = 1
        end
        if eventtype == 8 # Q->death
            ind_Q = linear_as[i, a, 5]
            ind_cumdead = linear_as[i, a, 9]
            dc[ind_Q, k] = -1
            dc[ind_cumdead, k] = 1
        end
        if eventtype == 9            #**** adding the event for E -> IQ
            ind_E = linear_as[i, a, 2]
            ind_DQ = linear_as[i, a, 10]
            ind_cumD = linear_as[i, a, 8]
            dc[ind_E, k] = -1
            dc[ind_DQ, k] = 1
            dc[ind_cumD] = 1
        end
        if eventtype == 10           #**** IQ->R
            ind_DQ = linear_as[i, a, 10]
            ind_R = linear_as[i, a, 6]
            dc[ind_DQ, k] = -1
            dc[ind_R, k] = 1
        end
    end
end

#This caps the Poisson processes at causing no more transitions than the group being effected
function max_change(out, u, p::CoVParameters_AS)
    for wa = 1:n_wa, a = 1:n_a, s = 1:n_s
        if u[wa, a, s] < 0
            u[wa, a, s] = 0
        end
        if out[linear_as_events[wa, a, s]] < 0
            out[linear_as_events[wa, a, s]] = 0
        end
    end
    for i = 1:n_wa, a = 1:n_a
        ind_trans = linear_as_events[i, a, 1]
        ind_EA = linear_as_events[i, a, 2]
        ind_ED = linear_as_events[i, a, 3]
        ind_IqQ = linear_as_events[i, a, 4]               #**** ind_DH became ind_IqQ
        ind_QR = linear_as_events[i, a, 5]
        ind_DR = linear_as_events[i, a, 6]
        ind_AR = linear_as_events[i, a, 7]
        ind_Qdeath = linear_as_events[i, a, 8]
        ind_EIq = linear_as_events[i, a, 9]                #****
        ind_IqR = linear_as_events[i, a, 10]               #****

        out[ind_trans] = min(out[ind_trans], u[i, a, 1])
        out[ind_QR] = min(out[ind_QR], u[i, a, 5])
        out[ind_AR] = min(out[ind_AR], u[i, a, 3])

        out[ind_DR] = min(out[ind_DR], u[i, a, 4])         #**** No more Iᴰ->R than actual Iᴰ
        out[ind_IqR] = min(out[ind_IqR], u[i, a, 10])       #**** No more IQ->R than actual IQ
        out[ind_IqQ] = min(out[ind_IqQ], u[i, a, 10])       #**** No more IQ->Q than actual IQ

        #spliting events using binomial sampling
        if out[ind_EA] + out[ind_ED] > u[i, a, 2] && u[i, a, 2] >= 0  #More incubations than actual rural E population
            #out[ind_ED] = rand(Binomial(u[i,a,2],p.δ)) #Binomially distributed the rural incubations between Asympotomatic and symptomatic
            #out[ind_EA] =  u[i,a,2] - out[ind_ED]
            D = rand(Binomial(u[i, a, 2], p.δ)) #Binomially distributed the incubations between Asympotomatic and symptomatic   #****
            out[ind_EA] = u[i, a, 2] - D
            out[ind_EIq] = rand(Binomial(D, p.τₚ)) # Binomially distribute symptomatics between DQ and D
            out[ind_ED] = D - out[ind_EIq]
        end
        #=if out[ind_DR] + out[ind_DQ] > u[i,a,4] # More end of infection events than diseased/symptomatic infecteds
            out[ind_DQ] = rand(Binomial(u[i,a,4],p.τ/(p.τ + p.γ))) #Binomially distributed the end of infections between hospitalisation and recovery
            out[ind_DR] = u[i,a,4] - out[ind_DQ]
        end=#
        if out[ind_IqQ] + out[ind_IqR] > u[i, a, 10]      #**** More IQ->Q + IQ->R than actual IQ
            out[ind_IqQ] = rand(Binomial(u[i, a, 10], p.τₚ))
            out[ind_IqR] = u[i, a, 10] - out[ind_IqQ]
        end
    end

    for wa = 1:n_wa, a = 1:n_a, s = 1:n_s
        if out[linear_as_events[wa, a, s]] < 0
            out[linear_as_events[wa, a, s]] = 0
        end
    end
end


function nonneg_tauleap(du, u, p::CoVParameters_AS, t)
    @unpack dc, dN, poi_rates, du_linear, τₚ, Κ_current, Κ_max_capacity = p
    rates(poi_rates, u, p, t)       #calculate rates of underlying Poisson processes
    PP_drivers(dN, poi_rates, p)    #Generate Poisson rvs with rates scaled by time step dt
    max_change(dN, u, p)            #Cap the size of the Poisson rvs to maintain non-negativity

    ##Contacts
    if τₚ != 0 && Κ_current < Κ_max_capacity    ## if the detection probability is not zero`
        #IQ_make_contacts(u,p,t)
        IQ_make_contacts_Nairobi(u, p, t)
        update_l_IQ(dN, u, p, t)                ##All IQ->R(10) are removed from l_IQ + All IQ with tau<=0 are added to dN and removed from l_IQ + decrease all tau + All E->IQ in dN are added to l_IQ:
        max_change(dN, u, p)
        update_contact_states(dN, u, p, t)      #updates the contact states, except for contacts made during this timestep
    elseif Κ_current >= Κ_max_capacity
        println("Tracing stopped at ", t)
    end

    mul!(du_linear, dc, dN)#Calculates the effect on the state in the inplace du vector
    du .= reshape(du_linear, n_wa, n_a, n_s)
    du .+= u #Calculates how the state should change
end

function ode_model(du, u, p::CoVParameters_AS, t)
    @unpack λ, γ, σ, δ, τ, μ₁, χ = p
    calculate_infection_rates!(u, p, t)
    for i = 1:n_wa, a = 1:n_a
        du[i, a, 1] = (-1) * χ[a] * λ[i, a] * u[i, a, 1]
        du[i, a, 2] = χ[a] * λ[i, a] * u[i, a, 1] - σ * u[i, a, 2]
        du[i, a, 3] = (1 - δ) * σ * u[i, a, 2] - γ * u[i, a, 3]
        du[i, a, 4] = δ * σ * u[i, a, 2] - γ * u[i, a, 4] - τ * u[i, a, 4]
        du[i, a, 5] = τ * u[i, a, 4] - γ * u[i, a, 5] - μ₁ * u[i, a, 5]
        du[i, a, 6] = γ * (u[i, a, 5] + u[i, a, 4] + u[i, a, 3])
        du[i, a, 7] = (1 - δ) * σ * u[i, a, 2]
        du[i, a, 8] = δ * σ * u[i, a, 2]
        du[i, a, 9] = μ₁ * u[i, a, 5]
    end
end

using PoissonRandom, SimpleRandom
function IQ_make_contacts(u, p::CoVParameters_AS, t::Float64)                     #****
    @unpack l_IQ, κ, Mₚ, dt, uₚ = p
    #### Recalculate uₚ=zeros(n_wa,n_a,n_s)   : Matrix of probabilities: when contacting someone with a specific wa and a, what is the chance of him being S, E, IQ, Iᴰ,... WE DO NOT MEET Q!
    for wa = 1:n_wa, a = 1:n_a, s = 1:n_s   ## Recalculate uₚ
        if s ∉ [5, 7, 8, 9] ## we only include the states 1(S) 2(E) 3(Iᴬ) 4(Iᴰ) 6(R) and 10(IQ). We do not include 5(Q) nor the cumulatives.
            uₚ[wa, a, s] = u[wa, a, s] / (sum(u[wa, a, 1:4]) + u[wa, a, 6] + u[wa, a, 10]) ## we only include the states 1(S) 2(E) 3(Iᴬ) 4(Iᴰ) 6(R) and 10(IQ). We do not include 5(Q) nor the cumulatives.
        else  # we have no chance of meeting an Q, and we remove the cumulative states
            uₚ[wa, a, s] = 0  # we have no chance of meeting an Q, and we remove the cumulative states
        end  # we have no chance of meeting an Q, and we remove the cumulative states
    end

    ## make contacts. We suppose all contacts are made in the same area. Contacts depend on the age mixing matrix M[a_i,a_j], we suppose a_i contacts a_j
    for i = 1:size(l_IQ, 1)
        n_contacts = pois_rand(κ * dt)  #number of contacts scaled by dt
        #i_IQ=l_IQ[i] ## is the current individual who's making the contacts
        append!(
            l_IQ[i].contacts,
            [Contact(
                wa = l_IQ[i].wa,
                a = random_choice(Mₚ[l_IQ[i].a, :]),
                s = random_choice(uₚ[l_IQ[i].wa, l_IQ[i].a, :]),
                contact_t = t,
                contact_dur = 0,
            ) for c = 1:n_contacts],
        )
    end
end

function IQ_make_contacts_Nairobi(u, p::CoVParameters_AS, t::Float64)                     #****
    @unpack l_IQ, κ, Mₚ, dt, uₚ = p
    #### Recalculate uₚ=zeros(n_wa,n_a,n_s)   : Matrix of probabilities: when contacting someone with a specific wa and a, what is the chance of him being S, E, IQ, Iᴰ,... WE DO NOT MEET Q!
    for wa in 4, a = 1:n_a, s = 1:n_s   ## Recalculate uₚ
        if s ∉ [5, 7, 8, 9] ## we only include the states 1(S) 2(E) 3(Iᴬ) 4(Iᴰ) 6(R) and 10(IQ). We do not include 5(Q) nor the cumulatives.  We have no chance of meeting an Q, and we remove the cumulative states
            uₚ[wa, a, s] = u[wa, a, s] / (sum(u[wa, a, 1:4]) + u[wa, a, 6] + u[wa, a, 10])
        else
            uₚ[wa, a, s] = 0
        end
    end

    ## make contacts. We suppose all contacts are made in the same area. Contacts depend on the age mixing matrix M[a_i,a_j], we suppose a_i contacts a_j
    for i = 1:size(l_IQ, 1)
        if (l_IQ[i].wa == 4)
            n_contacts = pois_rand(κ * dt)  #number of contacts scaled by dt
            #i_IQ=l_IQ[i] ## is the current individual who's making the contacts
            append!(
                l_IQ[i].contacts,
                [Contact(
                    wa = l_IQ[i].wa,
                    a = random_choice(Mₚ[l_IQ[i].a, :]),
                    s = random_choice(uₚ[l_IQ[i].wa, l_IQ[i].a, :]),
                    contact_t = t,
                    contact_dur = 0,
                ) for c = 1:n_contacts],
            )
        end
    end
end

function update_l_IQ(dN::Vector{Int64}, u, p::CoVParameters_AS, t)                 #****
    @unpack l_IQ, τ = p

    ##All events (10) IQ->R mean that an IQ recovered before being detected. Those need to be removed from l_IQ
    #IQ_to_R_events = dN[linear_as_events[:, :, 10]]
    for wa = 1:n_wa, a = 1:n_a
        n_Iq_to_R = dN[linear_as_events[wa, a, 10]]    ##the number of events IQ->R
        for i = 1:n_Iq_to_R
            selected_IQ_list = filter(x -> (x.wa == wa && x.a == a), l_IQ)
            IQ_indices = findall(x -> x == selected_IQ_list[rand(1:end)], l_IQ)  #draw one random IQ from l_IQ with same wa and a, and get its indices in l_IQ. This should be one element in an array
            if size(IQ_indices, 1) > 0
                deleteat!(l_IQ, IQ_indices[1])
            end       #delete the chosen element from l_IQ
            #if size(IQ_indices,1)>n_Iq_to_R      IQ_indices=IQ_indices[1:n_IIQ_to_R];    end    ## FORGOT why I wrote this line
        end
    end

    ## Update dN: All IQ with tau<=0 will be moved to Q with the event 4 (IQ->Q) and remove them from the list
    i = 1
    while i <= size(l_IQ, 1)
        if l_IQ[i].detection_dur <= 0
            dN[linear_as_events[Int(l_IQ[i].wa), Int(l_IQ[i].a), 4]] += 1          ## move IQ to Q, which is the event number 4
            intervention_trace_contacts(i, u, p, t)                                ## CONTACT TRACING INTERVENTION
            deleteat!(l_IQ, i)                                                   ## delete from the l_IQ list
        end
        i += 1
    end

    ## Decrease dt of all exponential durations of old IQ
    for i = 1:size(l_IQ, 1)
        if l_IQ[i].detection_dur - p.dt >= 0
            l_IQ[i].detection_dur -= p.dt
        else
            l_IQ[i].detection_dur = 0
        end
    end

    ## Add all new IQ to the list l_IQ l_IQ
    for wa = 1:n_wa, a = 1:n_a
        if dN[linear_as_events[wa, a, 9]] != 0
            for exp = 1:dN[linear_as_events[wa, a, 9]]     ##need to rescale by dt??
                push!(l_IQ, IQ_Person(wa = wa, a = a, detection_dur = randexp() / τ))     ##need to rescale by dt??
            end
        end
    end
end

function update_contact_states(dN::Vector{Int64}, u, p::CoVParameters_AS, t)       #**** Update the sates of all previously made contacts
    @unpack l_IQ, κₘ, dt = p

    for i = 1:size(l_IQ, 1)                         ## for all current IQ
        j = 1
        while j <= size(l_IQ[i].contacts, 1)        ## for all contacts
            e = l_IQ[i].contacts[j]                 # e is the contact

            if t - e.contact_t >= κₘ / dt           ## Forget the old contacts
                deleteat!(l_IQ[i].contacts, j)      ##forgetting the contact
            elseif e.contact_t != t                 ## if the contact hasnt been made during the current timestep
                if u[e.wa,e.a,1]!=0 && e.s == 1 && rand(Binomial(1,dN[linear_as_events[e.wa,e.a,1]]/u[e.wa,e.a,1]))==1 ## If contact is S(1) and event S->E(1) occured
                    e.s = 2
                end
                if u[e.wa, e.a, 2] != 0 && e.s == 2 ## if contact is E(2) then he can have E->Iᴬ(event 2), E->Iᴰ(3) or E->IQ(9)
                    p_Iᴬ = dN[linear_as_events[e.wa, e.a, 2]] / u[e.wa, e.a, 2]  # probability of becoming Iᴬ
                    p_Iᴰ = dN[linear_as_events[e.wa, e.a, 3]] / u[e.wa, e.a, 2]  # probability of becoming Iᴰ
                    p_IQ = dN[linear_as_events[e.wa, e.a, 9]] / u[e.wa, e.a, 2]  # probability of becoming IQ
                    p_E = 1 - p_Iᴬ - p_Iᴰ - p_IQ                                 #probability of staying E
                    new_ps = random_choice([p_E, p_Iᴬ, p_Iᴰ, p_IQ])
                    if new_ps == 2
                        e.s = 3
                    elseif new_ps == 3
                        e.s = 4
                    elseif new_ps == 4
                        e.s = 10
                    end
                end

                if u[e.wa, e.a, 3] != 0 && e.s == 3 && rand(Binomial(1, dN[linear_as_events[e.wa, e.a, 7]] / u[e.wa, e.a, 3])) == 1 ## If contact is Iᴬ(3) and event Iᴬ->R(7) occured
                    e.s = 6 ## If contact is Iᴬ(3) and event Iᴬ->R(7) occured
                end
                if u[e.wa, e.a, 4] != 0 && e.s == 4 && rand(Binomial(1, dN[linear_as_events[e.wa, e.a, 6]] / u[e.wa, e.a, 4])) == 1 ## If contact is Iᴰ(4) and event Iᴰ->R(6) occured
                    e.s = 6 ## If contact is Iᴰ(4) and event Iᴰ->R(6) occured
                end
                if u[e.wa, e.a, 5] != 0 && e.s == 5 && rand(Binomial(1, dN[linear_as_events[e.wa, e.a, 5]] / u[e.wa, e.a, 5])) == 1 ## If contact is Q(5) and event Q->R(5) occured
                    e.s = 6 ## If contact is Q(5) and event Q->R(5) occured
                end
                ## in the last line, we might never add a Q to the contactsm but an S or E or IQ might progress to Q with time
                ## if e.s==6 means e is an R, so he stays R
                if u[e.wa, e.a, 10] != 0 && e.s == 10 ## if contact is IQ then he can have IQ->Q(event 4) or IQ->R(event 10)
                    p_Q = dN[linear_as_events[e.wa, e.a, 4]] / u[e.wa, e.a, 10] # probability of becoming Q
                    p_R = dN[linear_as_events[e.wa, e.a, 10]] / u[e.wa, e.a, 10] # probability of becoming R
                    p_IQ = 1 - p_Q - p_R    # probability of staying IQ
                    if p_Q >= 0 && p_R >= 0 && p_IQ >= 0
                        new_ps = random_choice([p_IQ, p_Q, p_R])
                        if new_ps == 2
                            e.s = 5
                        elseif new_ps == 3
                            e.s = 6
                        end
                    end
                end
            end
            j += 1
        end
    end
end


function intervention_trace_contacts(traced_index::Int, u, p::CoVParameters_AS, t)   #****
    @unpack l_IQ, Δₜ, κ_per_event4, Κ_current, Κ_max_capacity, dt = p
    contacteds = [c for c in l_IQ[traced_index].contacts if t - c.contact_t <= Δₜ/dt]
    if contacteds != []
        shuffle!(contacteds)
        contacteds = @view contacteds[1:min(κ_per_event4, size(contacteds, 1))]      #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        for c in contacteds                                                     ## for each contact, remove the S, E, Iᴬ, Iᴰ or IQ to Q
            if c.s in [1, 2, 3, 4, 10] && u[c.wa, c.a, c.s] > 0 && Κ_current < Κ_max_capacity
                u[c.wa, c.a, c.s] -= 1
                u[c.wa, c.a, 5] += 1
                Κ_current += 1                                                  ## Increment the total number of traceds
            end
        end
    end
end

function check_negativity(u, dN::Vector{Int64}, t)
    for wa = 1:n_wa, a = 1:n_a, s = 1:n_s
        if u[wa, a, s] < 0
            u[wa, a, s] = 0
        end
        if dN[linear_as_events[wa, a, s]] < 0
            dN[linear_as_events[wa, a, s]] = 0
        end
    end
end
