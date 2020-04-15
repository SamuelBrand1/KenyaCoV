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
2 -> E                  #HH becomes Eᶠ (Exposed from first introduction into a HH)
3 -> I_subclinical      #HH becomes IAᶠ
4 -> I_diseased         #HH becomes IDᶠ
5 -> H(ospitalised)
6 -> Recovered
7 -> Cumulative I_sub
8 -> Cumulative I_dis
9 -> Cumulative Dead
10-> Eˢ                 #HH Eˢ (Exposed from subsequent infection)
11-> IAˢ                #HH
12-> IDˢ                #HH

Events for each wider area and age group:

1-> Transmission                                #HH becomes S->Eᶠ
2-> Incubation into asymptomatic E->A           #HH becomes Eᶠ->A
3-> Incubation into diseased E->D               #HH becomes Eᶠ->D
4-> Diseased become hospitalised/treated D->H
5-> Hospitalised/treated recover
6-> Diseased recover D->R
7-> Asymptomatics recover A->R
8-> Hospitalised/treated->death
9-> Transmission S->Eˢ                          #HH
10-> Eˢ->A                                      #HH
11-> Eˢ->D                                      #HH

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
    I_Aᶠ = @view u[:,:,3]
    I_Dᶠ = @view u[:,:,4]
    I_Aˢ = @view u[:,:,12]  #HH
    I_Dˢ = @view u[:,:,11]  #HH
    #=mul!(p.Î,p.T,p.ϵ*I_A .+ p.ϵ_D*I_D)
    p.Î[:,immobile_age_indices] .= p.ϵ*I_A[:,immobile_age_indices] .+ p.ϵ_D*I_D[:,immobile_age_indices]
    mul!(p.λ_loc,p.β*p.c_t(t).*(p.Î ./p.N̂),p.M)#Local force of infection due to age-mixing --- M is in to (row), from (col) format
    mul!(p.λ,p.T',p.λ_loc)#this accounts for mobile susceptibles contracting away from home
    p.λ[:,immobile_age_indices] .= p.λ_loc[:,immobile_age_indices]#This corrects for immobility
    p.λ[ind_mombasa_as,:] .+= p.ext_inf_rate*import_rate_mom(t,p.into_mom,p.global_prev)
    p.λ[ind_nairobi_as,:] .+= p.ext_inf_rate*import_rate_nai(t,p.into_nai,p.global_prev)=#

    #HH:
    mul!(p.Îᶠ,p.T,p.ϵ*I_Aᶠ .+ p.ϵ_D*I_Dᶠ .+ p.ϵ*I_Aˢ .+ p.ϵ_D*I_Dˢ)
    p.Îᶠ[:,immobile_age_indices] .= p.ϵ*I_Aᶠ[:,immobile_age_indices] .+ p.ϵ_D*I_Dᶠ[:,immobile_age_indices] .+ p.ϵ*I_Aˢ[:,immobile_age_indices] .+ p.ϵ_D*I_Dˢ[:,immobile_age_indices]
    mul!(p.λᶠ_loc,p.βᶠ*p.c_t(t).*(p.Îᶠ ./p.N̂),p.M)
    mul!(p.λᶠ,p.T',p.λᶠ_loc)
    p.λᶠ[:,immobile_age_indices] .= p.λᶠ_loc[:,immobile_age_indices]

    mul!(p.Îˢ,p.T,p.ϵ*I_Aˢ .+ p.ϵ_D*I_Dˢ)
    p.Îˢ[:,immobile_age_indices] .= p.ϵ*I_Aˢ[:,immobile_age_indices] .+ p.ϵ_D*I_Dˢ[:,immobile_age_indices]
    mul!(p.λˢ_loc,p.βˢ*p.c_t(t).*(p.Îˢ ./p.N̂),p.M)
    mul!(p.λˢ,p.T',p.λˢ_loc)
    p.λˢ[:,immobile_age_indices] .= p.λˢ_loc[:,immobile_age_indices]
    return nothing
end


function rates(out,u,p::CoVParameters_AS,t)
    @unpack #=λ=#λᶠ,λˢ,γ,σ,δ,τ,μ₁,χ,rel_detection_rate,clear_quarantine = p
    calculate_infection_rates!(u,p,t)
    for k = 1:length(out)
        i,a,eventtype = Tuple(index_as_events[k])
        if eventtype ==1
            out[k] = χ[a]*λᶠ[i,a]*u[i,a,1] #Transmission    #HH  S->Eᶠ
        end
        if eventtype ==2
            out[k] = (1-(δ*rel_detection_rate[a]))*σ*u[i,a,2] #HH becomes Eᶠ->A (without modifs here)
        end
        if eventtype ==3
            out[k] = δ*rel_detection_rate[a]*σ*u[i,a,2]     #HH becomes Eᶠ->D (without modifs here)
        end
        if eventtype ==4
            out[k] = τ*u[i,a,4] # D->H                      #HH Dᶠ->Q
        end
        if eventtype ==5
            out[k] = clear_quarantine*u[i,a,5] # H->R
        end
        if eventtype ==6
            out[k] = γ*u[i,a,4] # D->R                      #HH Dᶠ->R
        end
        if eventtype ==7
            out[k] = γ*u[i,a,3] # A->R                      #HH Aᶠ->R
        end
        if eventtype ==8
            out[k] = μ₁*u[i,a,5] # H->death
        end
        if eventtype ==9
            out[k] = χ[a]*λˢ[i,a]*u[i,a,1] #Transmission   #HH  S->Eˢ
        end
        if eventtype ==10
            out[k] = (1-(δ*rel_detection_rate[a]))*σ*u[i,a,10]         #HH  Eˢ->A
        end
        if eventtype ==11
            out[k] = δ*rel_detection_rate[a]*σ*u[i,a,10]    #HH  Eˢ->D
        end
        if eventtype ==12
            out[k] = τ*u[i,a,11]                           #HH Dˢ->Q
        end
        if eventtype ==13
            out[k] = γ*u[i,a,11]                           #HH Dˢ->R
        end
        if eventtype ==14
            out[k] = γ*u[i,a,12]                           #HH Aˢ->R
        end
    end
end

function change_matrix(dc)
    d1,d2 = size(dc)
    for k = 1:d2
        i,a,eventtype = Tuple(index_as_events[k])
        if eventtype ==1 #Transmission
            ind_S = linear_as[i,a,1]
            ind_E = linear_as[i,a,2]
            dc[ind_S,k] = -1
            dc[ind_E,k] = 1
        end
        if eventtype ==2 # E->A
            ind_E = linear_as[i,a,2]
            ind_A = linear_as[i,a,3]
            ind_cumA = linear_as[i,a,7]
            dc[ind_E,k] = -1
            dc[ind_A,k] = 1
            dc[ind_cumA,k] = 1
        end
        if eventtype ==3 # E->D
            ind_E = linear_as[i,a,2]
            ind_D = linear_as[i,a,4]
            ind_cumD = linear_as[i,a,8]
            dc[ind_E,k] = -1
            dc[ind_D,k] = 1
            dc[ind_cumD,k] = 1
        end
        if eventtype ==4# D->H
            ind_D = linear_as[i,a,4]
            ind_H = linear_as[i,a,5]
            dc[ind_D,k] = -1
            dc[ind_H,k] = 1
        end
        if eventtype ==5 # H->R
            ind_H = linear_as[i,a,5]
            ind_R = linear_as[i,a,6]
            dc[ind_H,k] = -1
            dc[ind_R,k] = 1
        end
        if eventtype ==6 # D->R
            ind_D = linear_as[i,a,4]
            ind_R = linear_as[i,a,6]
            dc[ind_D,k] = -1
            dc[ind_R,k] = 1
        end
        if eventtype ==7 # A->R
            ind_A = linear_as[i,a,3]
            ind_R = linear_as[i,a,6]
            dc[ind_A,k] = -1
            dc[ind_R,k] = 1
        end
        if eventtype ==8 # H->death
            ind_H = linear_as[i,a,5]
            ind_cumdead = linear_as[i,a,9]
            dc[ind_H,k] = -1
            dc[ind_cumdead,k] = 1
        end
        if eventtype ==9                        #HH  S->Eˢ
            ind_S = linear_as[i,a,1]
            ind_E = linear_as[i,a,10]           #HH ind_E here is Eˢ
            dc[ind_S,k] = -1
            dc[ind_E,k] = 1
        end
        if eventtype ==10                       #HH  Eˢ->A
            ind_E = linear_as[i,a,10]           #HH ind_E here is Eˢ
            ind_A = linear_as[i,a,3]
            ind_cumA = linear_as[i,a,7]
            dc[ind_E,k] = -1
            dc[ind_A,k] = 1
            dc[ind_cumA,k] = 1
        end
        if eventtype ==11                       #HH  Eˢ->D
            ind_E = linear_as[i,a,10]           #HH ind_E here is Eˢ
            ind_D = linear_as[i,a,4]
            ind_cumD = linear_as[i,a,8]
            dc[ind_E,k] = -1
            dc[ind_D,k] = 1
            dc[ind_cumD,k] = 1
        end
        if eventtype ==12                        #HH Dˢ->Q
            ind_D = linear_as[i,a,11]           #HH ind_D here is Dˢ
            ind_H = linear_as[i,a,5]
            dc[ind_D,k] = -1
            dc[ind_H,k] = 1
        end
        if eventtype ==13                       #HH Dˢ->R
            ind_D = linear_as[i,a,11]
            ind_R = linear_as[i,a,6]
            dc[ind_D,k] = -1
            dc[ind_R,k] = 1
        end
        if eventtype ==14                        #HH Aˢ->R
            ind_A = linear_as[i,a,12]
            ind_R = linear_as[i,a,6]
            dc[ind_A,k] = -1
            dc[ind_R,k] = 1
        end
    end
end

#This caps the Poisson processes at causing no more transitions than the group being effected
function max_change(out,u,p::CoVParameters_AS)
    for i = 1:n_wa,a = 1:n_a
        ind_trans = linear_as_events[i,a,1]
        ind_EA = linear_as_events[i,a,2]
        ind_ED = linear_as_events[i,a,3]
        ind_DH = linear_as_events[i,a,4]
        ind_HR = linear_as_events[i,a,5]
        ind_DR = linear_as_events[i,a,6]
        ind_AR = linear_as_events[i,a,7]
        ind_Hdeath = linear_as_events[i,a,8]
        #HH:
        ind_trans2 = linear_as_events[i,a,9]
        ind_EA2 = linear_as_events[i,a,10]
        ind_ED2 = linear_as_events[i,a,11]
        ind_DH2 = linear_as_events[i,a,12]
        ind_DR2 = linear_as_events[i,a,13]
        ind_AR2 = linear_as_events[i,a,14]


        #out[ind_trans] = min(out[ind_trans],u[i,a,1])  #HH (to be split binomially)
        out[ind_HR] = min(out[ind_HR],u[i,a,5])
        out[ind_AR] = min(out[ind_AR],u[i,a,3])

        #spliting events using binomial sampling
        if out[ind_EA] + out[ind_ED] > u[i,a,2] #More incubations than actual rural E population
            out[ind_ED] = rand(Binomial(u[i,a,2],p.δ)) #Binomially distributed the rural incubations between Asympotomatic and symptomatic
            out[ind_EA] =  u[i,a,2] - out[ind_ED]
        end
        if out[ind_DR] + out[ind_DH] > u[i,a,4] # More end of infection events than diseased/symptomatic infecteds
            out[ind_DH] = rand(Binomial(u[i,a,4],p.τ/(p.τ + p.γ))) #Binomially distributed the end of infections between hospitalisation and recovery
            out[ind_DR] = u[i,a,4] - out[ind_DH]
        end

        #HH:
        out[ind_AR2] = min(out[ind_AR2],u[i,a,12])
        if out[ind_trans] + out[ind_trans2] > u[i,a,1] #More transmissions than actual S population
            out[ind_trans] = rand(Binomial(u[i,a,1],p.λᶠ))
            out[ind_trans2] =  u[i,a,1] - out[ind_trans2]
        end
        if out[ind_EA2] + out[ind_ED2] > u[i,a,10] #More SUBSEQUENT incubations than actual Eˢ population
            out[ind_ED2] = rand(Binomial(u[i,a,10],p.δ)) #Binomially distributed the rural incubations between Asympotomatic and symptomatic
            out[ind_EA2] =  u[i,a,10] - out[ind_ED2]
        end
        if out[ind_DR2] + out[ind_DH2] > u[i,a,11] # More end of infection events than diseased/symptomatic infecteds
            out[ind_DH2] = rand(Binomial(u[i,a,11],p.τ/(p.τ + p.γ))) #Binomially distributed the end of infections between hospitalisation and recovery
            out[ind_DR2] = u[i,a,11] - out[ind_DH2]
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
    @unpack λ,γ,σ,δ,τ,μ₁,χ,rel_detection_rate,clear_quarantine = p
    calculate_infection_rates!(u,p,t)
    for i = 1:n_wa,a = 1:n_a
        du[i,a,1] = (-1)*χ[a]*λ[i,a]*u[i,a,1]
        du[i,a,2] = χ[a]*λ[i,a]*u[i,a,1] - σ*u[i,a,2]
        du[i,a,3] = (1-(δ*rel_detection_rate[a]))*σ*u[i,a,2] - γ*u[i,a,3]
        du[i,a,4] = δ*rel_detection_rate[a]*σ*u[i,a,2] - γ*u[i,a,4] - τ*u[i,a,4]
        du[i,a,5] = τ*u[i,a,4] - clear_quarantine*u[i,a,5] - μ₁*u[i,a,5]
        du[i,a,6] = clear_quarantine*u[i,a,5] + γ*(u[i,a,4] + u[i,a,3])
        du[i,a,7] = (1-δ)*σ*u[i,a,2]
        du[i,a,8] = δ*σ*u[i,a,2]
        du[i,a,9] = μ₁*u[i,a,5]
     end
end
