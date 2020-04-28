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
3 -> P
4 -> A
5 -> M(ild) symptomatics
6 -> (se)V(ere) symptomatics
7 -> H(ospitalised)
8 -> Recovered
9 -> Cumulative A
10-> Cumulative M
11-> Cumulative V

Events for each wider area and age group:

1-> S to E
2-> E to P
3-> P to A
4-> P to M
5-> P to V
6-> V to H
7-> M to R
8-> A to R

"""

dc_age = zeros(Int64,n_wa*n_a*n_s,n_ta*n*n_a)
# asymp_indices = zeros(Bool,n_wa,n_a,n_s)
# asymp_indices[:,:,3] .= true;
# f_asymp_indices = findall(asymp_indices[:])
# diseased_indices = zeros(Bool,n_wa,n_a,n_s)
# diseased_indices[:,:,4] .= true;
# f_diseased_indices = findall(diseased_indices[:])
# asymp_indices = 0;#free memory
# diseased_indices = 0;

function calculate_infection_rates!(u,p::CoVParameters_AS,t)
    I_A = @view u[:,:,4]
    I_M = @view u[:,:,5]
    I_V = @view u[:,:,6]
    mul!(p.Î,p.T,p.ϵ*I_A .+ p.ϵ_D*I_M .+ p.ϵ*I_V)  #Local infecteds **if** everyone moved around
    p.Î[:,immobile_age_indices] .= p.ϵ*I_A[:,immobile_age_indices] .+ p.ϵ_D*I_M[:,immobile_age_indices] .+ p.ϵ_D*I_V[:,immobile_age_indices]#This corrects for immobility
    mul!(p.λ_loc,p.β*p.c_t(t).*(p.Î ./p.N̂),p.M)#Local force of infection due to age-mixing --- M is in to (row), from (col) format
    mul!(p.λ,p.T',p.λ_loc)#this accounts for mobile susceptibles contracting away from home
    p.λ[:,immobile_age_indices] .= p.λ_loc[:,immobile_age_indices]#This corrects for immobility
    p.λ[ind_mombasa_as,:] .+= p.ext_inf_rate*import_rate_mom(t,p.into_mom,p.global_prev)
    p.λ[ind_nairobi_as,:] .+= p.ext_inf_rate*import_rate_nai(t,p.into_nai,p.global_prev)
    return nothing
end


function rates(out,u,p::CoVParameters_AS,t)
    @unpack λ,γ,σ,δ,τ,μ₁,χ,rel_detection_rate,clear_quarantine,hₐ = p
    calculate_infection_rates!(u,p,t)
    for k = 1:length(out)
        i,a,eventtype = Tuple(index_as_events[k])
        if eventtype ==1
            out[k] = χ[a]*λ[i,a]*u[i,a,1] #Transmission #S->E
        end
        if eventtype ==2
            out[k] = σ*u[i,a,2] #E->P
        end
        if eventtype ==3
            #out[k] = (1-(δ*rel_detection_rate[a]))*u[i,a,3] #P->A
            out[k] = (1-δ)*u[i,a,3] #P->A
        end
        if eventtype ==4
            #out[k] = δ*(1-hₐ)*rel_detection_rate[a]*u[i,a,3] #P->M
            out[k] = δ*(1-hₐ)*u[i,a,3] #P->M
        end
        if eventtype ==5
            #out[k] = δ*hₐ*rel_detection_rate[a]*u[i,a,3] #P->V
            out[k] = δ*hₐ*u[i,a,3] #P->V
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
    end
end

function change_matrix(dc)
    d1,d2 = size(dc)
    for k = 1:d2
        i,a,eventtype = Tuple(index_as_events[k])
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
            dc[ind_V,k] = -1
            dc[ind_H,k] = 1
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
    end
end

#This caps the Poisson processes at causing no more transitions than the group being effected
function max_change(out,u,p::CoVParameters_AS)
    for i = 1:n_wa,a = 1:n_a
        ind_trans = linear_as_events[i,a,1]
        ind_EP = linear_as_events[i,a,2]
        ind_PA = linear_as_events[i,a,3]
        ind_PM = linear_as_events[i,a,4]
        ind_PV = linear_as_events[i,a,5]
        ind_VH = linear_as_events[i,a,6]
        ind_MR = linear_as_events[i,a,7]
        ind_AR = linear_as_events[i,a,8]


        out[ind_trans] = min(out[ind_trans],u[i,a,1])
        out[ind_EP] = min(out[ind_EP],u[i,a,2])
        out[ind_VH] = min(out[ind_VH],u[i,a,6])
        out[ind_MR] = min(out[ind_MR],u[i,a,5])
        out[ind_AR] = min(out[ind_AR],u[i,a,4])

        #spliting events using binomial sampling
        if out[ind_PA] + out[ind_PM] + out[ind_PV] > u[i,a,3] #More cases than actual P
            D = rand(Binomial(u[i, a, 3], p.δ))
            out[ind_PA] = u[i, a, 3] - D
            out[ind_PV] = rand(Binomial(D, p.hₐ))
            out[ind_PM] = D - out[ind_PM]
        end
    end
end


function nonneg_tauleap(du,u,p::CoVParameters_AS,t)
    @unpack dc,dN,poi_rates,du_linear = p
    rates(poi_rates,u,p,t) #calculate rates of underlying Poisson processes
    PP_drivers(dN,poi_rates,p)#Generate Poisson rvs with rates scaled by time step dt
    max_change(dN,u,p)#Cap the size of the Poisson rvs to maintain non-negativity
    mul!(du_linear,dc,dN)#Calculates the effect on the state in the inplace du vector
    du .= reshape(du_linear,n_wa,n_a,n_s)
    du .+= u #Calculates how the state should change
end

function ode_model(du,u,p::CoVParameters_AS,t)
    @unpack λ,γ,σ,δ,τ,μ₁,χ,rel_detection_rate,clear_quarantine,hₐ = p
    calculate_infection_rates!(u,p,t)
    for i = 1:n_wa,a = 1:n_a
        du[i,a,1] = (-1)*χ[a]*λ[i,a]*u[i,a,1]
        du[i,a,2] = χ[a]*λ[i,a]*u[i,a,1] - σ*u[i,a,2]
        du[i,a,3] = σ*u[i,a,2] - rel_detection_rate[a]*u[i,a,3]
        du[i,a,4] = (1-(δ*rel_detection_rate[a]))*u[i,a,3] - γ*u[i,a,4]
        du[i,a,5] = δ*(1-hₐ)*rel_detection_rate[a]*u[i,a,3] - γ*u[i,a,5]
        du[i,a,6] = δ*hₐ*rel_detection_rate[a]*u[i,a,3] - τ*u[i,a,6]
        du[i,a,7] = τ*u[i,a,6]
        du[i,a,8] = γ*u[i,a,5] + γ*u[i,a,4]
     end
end
