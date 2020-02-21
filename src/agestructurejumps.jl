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
4 -> I_diseased
5 -> H(ospitalised)
6 -> Recovered
7 -> Cumulative I_sub
8 -> Cumulative I_dis
9 -> Cumulative Dead

Events for each wider area and age group:

1-> Transmission
2-> Incubation into asymptomatic E->A
3-> Incubation into diseased E->D
4-> Diseased become hospitalised/treated D->H
5-> Hospitalised/treated recover
6-> Diseased recover D->R
7-> Asymptomatics recover A->R
8->Hospitalised/treated->death

"""

dc_age = zeros(Int64,n_wa*n_a*n_s,n_ta*n*n_a)
asymp_indices = zeros(Bool,n_wa,n_a,n_s)
asymp_indices[:,:,3] .= true;
f_asymp_indices = findall(asymp_indices[:])
diseased_indices = zeros(Bool,n_wa,n_a,n_s)
diseased_indices[:,:,4] .= true;
f_diseased_indices = findall(diseased_indices[:])
asymp_indices = 0;#free memory
diseased_indices = 0;

mobile_age_indices = 5:11; #This assumes that 16-49 year ols move around and others don't
immobile_age_indices = [1,2,3,4,12,13,14,15,16]
function calculate_infection_rates!(u,p::CoVParameters_AS,t)
    I_A = @view u[f_asymp_indices]
    I_D = @view u[f_diseased_indices]
    I_A = reshape(I_A,n_wa,n_a)
    I_D = reshape(I_D,n_wa,n_a)
    mul!(p.Î,p.T,p.ϵ*I_A .+ I_D)#Local infecteds **if** everyone moved around
    p.Î[:,immobile_age_indices] .= p.ϵ*I_A[:,immobile_age_indices] .+ I_D[:,immobile_age_indices]#This corrects for immobility
    mul!(p.λ_loc,p.M,p.β .*(p.Î ./p.N̂))#Local force of infection due to age-mixing
    mul!(p.λ,p.T',p.λ_loc)#this accounts for mobile susceptibles contracting away from home
    p.λ[:,immobile_age_indices] .= p.λ_loc[:,immobile_age_indices]#This corrects for immobility
    p.λ[ind_mombasa_as] += p.ext_inf_rate*import_rate_mom(t,p.into_mom,p.global_prev)
    p.λ[ind_nairobi_as] += p.ext_inf_rate*import_rate_nai(t,p.into_nai,p.global_prev)
    return nothing
end
