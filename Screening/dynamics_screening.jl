"""
Basic representation of the state of the age structured model:
u[wa_index,age_group,disease state]

This is row major unpacked so
first 1,...,n entries are 0-4 year old susceptibles in wider areas 1,...,n
then n+1,...,2n are 5-9 year old susceptibles in wider areas 1,...,n
.
.
.
then n*n_s + 1, ...,  n*n_s + n entries are 0-4 year old exposed in wider areas 1,...,n
.
.
.

States:
1 -> S(usceptibles)
2 -> E(xposed/latent infected)
3 -> P(re-symptomatic infected)
4 -> A(symptomatic)
5 -> M(ild) symptomatics
6 -> first mild then eventually (se)V(ere) symptomatics
7 -> H(ospitalised)
8 -> Recovered
9 -> Cumulative P->A
10-> Cumulative P->M
11-> Cumulative P->V
12 -> Cumulative V->H
13 -> Quarantine
14 -> Q(uarantined) (se)V(ere) Qᵥ
15 -> Quarantined Susceptibles Qₛ
16 -> Cumulative Q (Q+Qᵥ+Qₛ)


Events for each wider area and age group:

1-> S to E
2-> E to P
3-> P to A
4-> P to M
5-> P to V
6-> V to H
7-> M to R
8-> A to R
9-> Q to R
10-> Qₛ to S
11-> Qᵥ to H
12-> E to Q
13->P to Q
14->A to Q
15->M to Q
16->V to Qᵥ
17-> R to Q
18-> S to Qₛ

"""

#dc_age = zeros(Int64,n_wa*n_a*n_s,n_ta*n*n_a)

function import_rate_mom(t,into_mom,global_prev)
    if t+1>min(length(into_mom),length(global_prev))
        t_int=min(length(into_mom),length(global_prev))
    else
        t_int=Int(floor(t))+1
    end
    return into_mom[t_int]*global_prev[t_int]
end

function import_rate_nai(t,into_nai,global_prev)
    if t+1>min(length(into_nai),length(global_prev))
        t_int=min(length(into_nai),length(global_prev))
    else
        t_int=Int(floor(t))+1
    end
    return into_nai[t_int]*global_prev[t_int]
end

"""
    function calculate_infection_rates!(u,p::CoVParameters_Screening,t)

Inplace method for calculating the force of infection on susceptibles in each spatial and age group.
"""
function calculate_infection_rates!(u,p::CoVParameters_Screening,t)
    I_P = @view u[:,:,3]
    I_A = @view u[:,:,4]
    I_M = @view u[:,:,5]
    I_V = @view u[:,:,6]
    mul!(p.Î,p.T,p.ϵ*I_P .+ p.ϵ*I_A .+ p.ϵ_D*I_M .+ p.ϵ_V*I_V)  #Local infecteds **if** everyone moved around
    p.Î[:,immobile_age_indices] .= p.ϵ*I_P[:,immobile_age_indices] .+ p.ϵ*I_A[:,immobile_age_indices] .+ p.ϵ_D*I_M[:,immobile_age_indices] .+ p.ϵ_V*I_V[:,immobile_age_indices]#This corrects for immobility
    mul!(p.λ_loc,p.β*p.c_t(t).*(p.Î ./p.N̂),p.M)#Local force of infection due to age-mixing --- M is in to (row), from (col) format
    mul!(p.λ,p.T',p.λ_loc)#this accounts for mobile susceptibles contracting away from home
    p.λ[:,immobile_age_indices] .= p.λ_loc[:,immobile_age_indices]#This corrects for immobility of susceptibles
    p.λ[ind_mombasa_as,:] .+= p.ext_inf_rate*import_rate_mom(t,p.into_mom,p.global_prev)
    p.λ[ind_nairobi_as,:] .+= p.ext_inf_rate*import_rate_nai(t,p.into_nai,p.global_prev)
    return nothing
end

"""
    function rates(out,u,p::CoVParameters_Screening,t)

Inplace method for calculating the rate of each of the 8 events per spatial and age group.
    1-> S to E
    2-> E to P
    3-> P to A
    4-> P to M
    5-> P to V
    6-> V to H
    7-> M to R
    8-> A to R
The out vector is in linear index form.
"""
function rates(out,u,p::CoVParameters_Screening,t)
    @unpack λ,γ,σ₁,σ₂,δ,τ,μ₁,χ,rel_detection_rate,clear_quarantine,hₐ = p
    calculate_infection_rates!(u,p,t)
    for k = 1:length(out)
        i,a,eventtype = Tuple(index_as_events[k])
        if eventtype ==1
            out[k] = χ[a]*λ[i,a]*u[i,a,1] #Transmission #S->E
        end
        if eventtype ==2
            out[k] = σ₁*u[i,a,2] #E->P
        end
        if eventtype ==3
            out[k] = σ₂*(1-(δ*rel_detection_rate[a]))*u[i,a,3] #P->A
        end
        if eventtype ==4
            out[k] = σ₂*δ*(1-hₐ[a])*rel_detection_rate[a]*u[i,a,3] #P->M
        end
        if eventtype ==5
            out[k] = σ₂*δ*hₐ[a]*rel_detection_rate[a]*u[i,a,3] #P->V
        end
        if eventtype ==6
            out[k] = τ*u[i,a,6] # V->H
        end
        if eventtype ==7
            out[k] = γ*u[i,a,5] # M->R
        end
        if eventtype ==8
            out[k] = γ*u[i,a,4] # A->R
        end
        if eventtype ==9
            out[k] = 1/clear_quarantine*u[i,a,13] # Q->R
        end
        if eventtype ==10
            out[k] = 1/clear_quarantine*u[i,a,15] # Qₛ->S
        end
        if eventtype ==11
            out[k] = τ*u[i,a,14] # Qᵥ to H
        end
    end
end

"""
    function change_matrix(dc)

Inplace method for creating a change matrix, that is a matrix encoding the linear map between frequency of events
    occuring and the change this causes in the state of the epidemic. Both the list of events and the list of states
    are encoded in a linear index.
"""
function change_matrix(dc)
    d1,d2 = size(dc)
    for k = 1:d2 #loop over
        i,a,eventtype = Tuple(index_as_events[k]) #
        if eventtype ==1 #S->E
            ind_S = linear_as[i,a,1]
            ind_E = linear_as[i,a,2]
            dc[ind_S,k] = -1
            dc[ind_E,k] = 1
        end
        if eventtype ==2 #E->P
            ind_E = linear_as[i,a,2]
            ind_P = linear_as[i,a,3]
            dc[ind_E,k] = -1
            dc[ind_P,k] = 1
        end
        if eventtype ==3 # P->A
            ind_P = linear_as[i,a,3]
            ind_A = linear_as[i,a,4]
            ind_cumA = linear_as[i,a,9]
            dc[ind_P,k] = -1
            dc[ind_A,k] = 1
            dc[ind_cumA,k] = 1
        end
        if eventtype ==4 # P->M
            ind_P = linear_as[i,a,3]
            ind_M = linear_as[i,a,5]
            ind_cumM = linear_as[i,a,10]
            dc[ind_P,k] = -1
            dc[ind_M,k] = 1
            dc[ind_cumM,k] = 1
        end
        if eventtype ==5 # P->V
            ind_P = linear_as[i,a,3]
            ind_V = linear_as[i,a,6]
            ind_cumV = linear_as[i,a,11]
            dc[ind_P,k] = -1
            dc[ind_V,k] = 1
            dc[ind_cumV,k] = 1
        end
        if eventtype ==6# V->H
            ind_V = linear_as[i,a,6]
            ind_H = linear_as[i,a,7]
            ind_cumVH = linear_as[i,a,12]
            dc[ind_V,k] = -1
            dc[ind_H,k] = 1
            dc[ind_cumVH,k] = 1
        end
        if eventtype ==7 # M->R
            ind_M = linear_as[i,a,5]
            ind_R = linear_as[i,a,8]
            dc[ind_M,k] = -1
            dc[ind_R,k] = 1
        end
        if eventtype ==8 # A->R
            ind_A = linear_as[i,a,4]
            ind_R = linear_as[i,a,8]
            dc[ind_A,k] = -1
            dc[ind_R,k] = 1
        end
        if eventtype ==9 # Q->R
            ind_Q = linear_as[i,a,13]
            ind_R = linear_as[i,a,8]
            dc[ind_Q,k] = -1
            dc[ind_R,k] = 1
        end
        if eventtype ==10 # Qₛ->S
            ind_Qₛ = linear_as[i,a,15]
            ind_S = linear_as[i,a,1]
            dc[ind_Qₛ,k] = -1
            dc[ind_S,k] = 1
        end
        if eventtype ==11# # Qᵥ -> H
            ind_Qᵥ = linear_as[i,a,14]
            ind_H = linear_as[i,a,7]
            ind_cumVH = linear_as[i,a,12]
            dc[ind_Qᵥ,k] = -1
            dc[ind_H,k] = 1
            dc[ind_cumVH,k] = 1
        end
        #Screening AND/OR Contact contracting
        if eventtype ==12# E->Q via screening or CT
            ind_E = linear_as[i,a,2]
            ind_Q = linear_as[i,a,13]
            ind_cumQ = linear_as[i,a,16]
            dc[ind_E,k] = -1
            dc[ind_Q,k] = 1
            dc[ind_cumQ,k] = 1
        end
        if eventtype ==13# P->Q via screening or CT
            ind_P = linear_as[i,a,3]
            ind_Q = linear_as[i,a,13]
            ind_cumQ = linear_as[i,a,16]
            dc[ind_P,k] = -1
            dc[ind_Q,k] = 1
            dc[ind_cumQ,k] = 1
        end
        if eventtype ==14# A->Q via screening or CT
            ind_A = linear_as[i,a,4]
            ind_Q = linear_as[i,a,13]
            ind_cumQ = linear_as[i,a,16]
            dc[ind_A,k] = -1
            dc[ind_Q,k] = 1
            dc[ind_cumQ,k] = 1
        end
        if eventtype ==15# M->Q via screening or CT
            ind_M = linear_as[i,a,5]
            ind_Q = linear_as[i,a,13]
            ind_cumQ = linear_as[i,a,16]
            dc[ind_M,k] = -1
            dc[ind_Q,k] = 1
            dc[ind_cumQ,k] = 1
        end
        if eventtype ==16# V->Qᵥ via screening or CT
            ind_V = linear_as[i,a,6]
            ind_Qᵥ = linear_as[i,a,14]
            ind_cumQ = linear_as[i,a,16]
            dc[ind_V,k] = -1
            dc[ind_Qᵥ,k] = 1
            dc[ind_cumQ,k] = 1
        end
        if eventtype ==17# R->Q via CT
            ind_R = linear_as[i,a,8]
            ind_Q = linear_as[i,a,13]
            ind_cumQ = linear_as[i,a,16]
            dc[ind_R,k] = -1
            dc[ind_Q,k] = 1
            dc[ind_cumQ,k] = 1
        end
        if eventtype ==18# S->Qₛ via CT
            ind_S = linear_as[i,a,1]
            ind_Qₛ = linear_as[i,a,15]
            ind_cumQ = linear_as[i,a,16]
            dc[ind_S,k] = -1
            dc[ind_Qₛ,k] = 1
            dc[ind_cumQ,k] = 1
        end
    end
end

"""
    function PP_drivers(dN::Vector{Int64},rates,p)
This method in-place generates the crude number of each type of event proposed by a Poisson process with
rates calculated by rates(out,u,p::CoVParameters_Screening,t).
"""
function PP_drivers(dN::Vector{Int64},rates,p)
    for i = 1:length(dN)
        if rates[i] >= 0.
                dN[i] = rand(Poisson(p.dt*rates[i]))
        else
            dN[i] = 0
        end
    end
end

"""
    function max_change(out,u,p::CoVParameters_Screening)
This method in-place modifies the number of each type of event proposed by the Poisson process
    so that non-negativity is respected.
"""
function max_change(out,u,p::CoVParameters_Screening)
    @unpack δ,rel_detection_rate,hₐ = p
    for i = 1:n,a = 1:n_a
        ind_trans = linear_as_events[i,a,1]
        ind_EP = linear_as_events[i,a,2]
        ind_PA = linear_as_events[i,a,3]
        ind_PM = linear_as_events[i,a,4]
        ind_PV = linear_as_events[i,a,5]
        ind_VH = linear_as_events[i,a,6]
        ind_MR = linear_as_events[i,a,7]
        ind_AR = linear_as_events[i,a,8]
        #For interventions:
        ind_QR = linear_as_events[i,a,9]
        ind_QₛS = linear_as_events[i,a,10]
        ind_QᵥH = linear_as_events[i,a,11]
        ind_EQ = linear_as_events[i,a,12]
        ind_PQ = linear_as_events[i,a,13]
        ind_AQ = linear_as_events[i,a,14]
        ind_MQ = linear_as_events[i,a,15]
        ind_VQᵥ = linear_as_events[i,a,16]
        ind_RQ = linear_as_events[i,a,17]
        ind_SQₛ = linear_as_events[i,a,18]

        #out[ind_trans] = min(out[ind_trans],u[i,a,1])
        if out[ind_trans] + out[ind_SQₛ] > u[i,a,1]     #priority to event SQₛ
            out[ind_SQₛ]=min(out[ind_SQₛ],u[i,a,1])
            out[ind_trans]=u[i,a,1]-out[ind_SQₛ]
        end
        #out[ind_EP] = min(out[ind_EP],u[i,a,2])
        if out[ind_EP] + out[ind_EQ] > u[i,a,2]     #priority to event EQ
            out[ind_EQ]=min(out[ind_EQ],u[i,a,2])
            out[ind_EP]=u[i,a,2]-out[ind_EQ]
        end
        #out[ind_VH] = min(out[ind_VH],u[i,a,6])
        if out[ind_VH] + out[ind_VQᵥ] > u[i,a,6]     #priority to event VQᵥ
            out[ind_VQᵥ]=min(out[ind_VQᵥ],u[i,a,6])
            out[ind_VH]=u[i,a,6]-out[ind_VQᵥ]
        end
        #out[ind_MR] = min(out[ind_MR],u[i,a,5])
        if out[ind_MR] + out[ind_MQ] > u[i,a,5]     #priority to event MQ
            out[ind_MQ]=min(out[ind_MQ],u[i,a,5])
            out[ind_MR]=u[i,a,5]-out[ind_MQ]
        end
        #out[ind_AR] = min(out[ind_AR],u[i,a,4])
        if out[ind_AR] + out[ind_AQ] > u[i,a,4]     #priority to event AQ
            out[ind_AQ]=min(out[ind_AQ],u[i,a,4])
            out[ind_AR]=u[i,a,4]-out[ind_AQ]
        end
        out[ind_QR] = min(out[ind_QR],u[i,a,13])
        out[ind_QᵥH] = min(out[ind_QᵥH],u[i,a,14])
        out[ind_QₛS] = min(out[ind_QₛS],u[i,a,15])

        #splitting events using multinomial sampling
        if out[ind_PA] + out[ind_PM] + out[ind_PV] + out[ind_PQ] > u[i,a,3] #More transitions than actual P so all P individuals transition
            out[ind_PQ] = min(out[ind_PQ],u[i,a,3])     #priority to event PQ
            rel_rate_each_event = [1-(δ*rel_detection_rate[a]),δ*(1-hₐ[a])*rel_detection_rate[a],δ*hₐ[a]*rel_detection_rate[a]]
            D = rand(Multinomial(u[i, a, 3] - out[ind_PQ], LinearAlgebra.normalize!(rel_rate_each_event,1) )) #Draw the states of the P individuals after transition
            out[ind_PA] = D[1]
            out[ind_PM] = D[2]
            out[ind_PV] = D[3]
        end
    end
end

"""
    nonneg_tauleap(du,u,p::CoVParameters_Screening,t)

This performs one in-place simulation of a time step
"""
function nonneg_tauleap(du,u,p::CoVParameters_Screening,t)
    @unpack dc,dN,poi_rates,du_linear,dN_S = p
    if t%1==0 && t>0
        push_intervention(dN_S,p,t)
        max_change(dN_S,u,p)
        mul!(du_linear,dc,dN_S)#Calculates the effect on the state in the inplace du vector
        du .= reshape(du_linear,n,n_a,n_s)
        u .+= du #Calculates how the state should change
        fill!(dN_S,0)
    end

    rates(poi_rates,u,p,t) #calculate rates of underlying Poisson processes
    PP_drivers(dN,poi_rates,p)#Generate Poisson rvs with rates scaled by time step dt
    max_change(dN,u,p)#Cap the size of the Poisson rvs to maintain non-negativity
    mul!(du_linear,dc,dN)#Calculates the effect on the state in the inplace du vector
    du .= reshape(du_linear,n,n_a,n_s)
    du .+= u #Calculates how the state should change
end

"""
    function ode_model(du,u,p::CoVParameters_Screening,t)

This is the vector field of the related KenyaCoV ODE model
"""
function ode_model(du,u,p::CoVParameters_Screening,t)
    @unpack λ,γ,σ₁,σ₂,δ,τ,μ₁,χ,rel_detection_rate,hₐ = p
    calculate_infection_rates!(u,p,t)
    for i = 1:n,a = 1:n_a
        du[i,a,1] = (-1)*χ[a]*λ[i,a]*u[i,a,1]
        du[i,a,2] = χ[a]*λ[i,a]*u[i,a,1] - σ₁*u[i,a,2]
        du[i,a,3] = σ₁*u[i,a,2] - σ₂*u[i,a,3]
        du[i,a,4] = σ₂*(1-(δ*rel_detection_rate[a]))*u[i,a,3] - γ*u[i,a,4]
        du[i,a,5] = σ₂*δ*(1-hₐ[a])*rel_detection_rate[a]*u[i,a,3] - γ*u[i,a,5]
        du[i,a,6] = σ₂*δ*hₐ[a]*rel_detection_rate[a]*u[i,a,3] - τ*u[i,a,6]
        du[i,a,7] = τ*u[i,a,6]
        du[i,a,8] = γ*u[i,a,5] + γ*u[i,a,4]
        du[i,a,9] = σ₂*(1-(δ*rel_detection_rate[a]))*u[i,a,3]
        du[i,a,10] = σ₂*δ*(1-hₐ[a])*rel_detection_rate[a]*u[i,a,3]
        du[i,a,11] = σ₂*δ*hₐ[a]*rel_detection_rate[a]*u[i,a,3]
        du[i,a,12] = τ*u[i,a,6] #Same as being hosp. for this version
     end
end

############ Call intervention into nonneg tau leap
function push_intervention(dN,p,t)
    @unpack toQ=p
    if t%1==0 && t>0 && sum(toQ)!=0   #Contact tracing of hospitalized only  OR Mass screening OR
       for r=1:n,a=1:n_a
           if sum(toQ[r,a,:])!=0
               if toQ[r,a,1]!=0    dN[linear_as_events[r,a,18]]=toQ[r,a,1];    end
               if toQ[r,a,2]!=0    dN[linear_as_events[r,a,12]]=toQ[r,a,2];    end
               if toQ[r,a,3]!=0    dN[linear_as_events[r,a,13]]=toQ[r,a,3];    end
               if toQ[r,a,4]!=0    dN[linear_as_events[r,a,14]]=toQ[r,a,4];    end
               if toQ[r,a,5]!=0    dN[linear_as_events[r,a,15]]=toQ[r,a,5];    end
               if toQ[r,a,6]!=0    dN[linear_as_events[r,a,16]]=toQ[r,a,6];    end
               if toQ[r,a,8]!=0    dN[linear_as_events[r,a,17]]=toQ[r,a,8];    end
           end
       end
       fill!(toQ, 0)
    end
end
