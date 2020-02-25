"""
For this code to work the 3-dim representation should
be unpacked as a vector

This is row major unpacked so
first 1,...,n entries are susceptibles in urban areas 1,...,n
then n+1,...,2n are exposeds in urban areas 1,...,n
.
.
.
then n*num_states + 1, ..., n*num_states + n entries are susceptibles in rural areas 1,...,n
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

Events for each location/county:

1-> Urban transmission
2-> Rural transmission
3-> #urban E->A
4-> rural E->A
5-> urban E->D
6-> rural E->D
7-> urban D->H
8-> rural D->H
9-> urban H->R
10-> rural H->R
11-> urban D->R
12-> rural D->R
13-> urban A->R
14-> rural A->R
15-> urban H->death
16-> rural H->death

"""
dc = zeros(Int64,n*n_s*2,n_t*n)


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

function calculate_infection_rates!(u,p::CoVParameters,t)
    I_urb_A = @view u[((3-1)*n + 1):((3-1)*n + n)]
    I_urb_D = @view u[((4-1)*n + 1):((4-1)*n + n)]
    I_rur_A = @view u[((3-1)*n + n_s*n + 1):((3-1)*n + n_s*n + n)]
    I_rur_D = @view u[((4-1)*n + n_s*n + 1):((4-1)*n + n_s*n + n)]
    mul!(p.Î,p.T,p.ϵ*I_urb_A .+ I_urb_D  )
    p.Î .+=  p.ϵ*I_rur_A .+ I_rur_D
    mul!(p.λ_urb,p.T',p.β .*(p.Î ./p.N̂))
    p.λ_urb[ind_mombasa] += p.ext_inf_rate*import_rate_mom(t,p.into_mom,p.global_prev)
    p.λ_urb[ind_nairobi] += p.ext_inf_rate*import_rate_nai(t,p.into_nai,p.global_prev)
    p.λ_rur .= p.β .*(p.Î ./p.N̂)
    return nothing
end

function rates(out,u,p::CoVParameters,t)
    @unpack λ_urb,λ_rur,β,γ,σ,δ,τ,μ₁,into_mom,into_nai,global_prev = p
    calculate_infection_rates!(u,p,t)
    for i = 1:n
        out[(i-1)*n_t+1] = λ_urb[i]*u[(1-1)*n + i] #urban transmission
        out[(i-1)*n_t+2] =  λ_rur[i]*u[(1-1)*n + n_s*n + i] #rural transmission
        out[(i-1)*n_t+3] = (1-δ)*σ*u[(2-1)*n + i] #urban E->A
        out[(i-1)*n_t+4] =  (1-δ)*σ*u[(2-1)*n + n_s*n + i] #rural E->A
        out[(i-1)*n_t+5] = δ*σ*u[(2-1)*n + i] #urban E->D
        out[(i-1)*n_t+6] =  δ*σ*u[(2-1)*n + n_s*n + i] #rural E->D
        out[(i-1)*n_t+7] = τ*u[(4-1)*n + i] #urban D->H
        out[(i-1)*n_t+8] =  τ*u[(4-1)*n + n_s*n + i] #rural D->H
        out[(i-1)*n_t+9] = γ*u[(5-1)*n + i] #urban H->R
        out[(i-1)*n_t+10] =  γ*u[(5-1)*n + n_s*n + i] #rural H->R
        out[(i-1)*n_t+11] = γ*u[(4-1)*n + i] #urban D->R
        out[(i-1)*n_t+12] =  γ*u[(4-1)*n + n_s*n + i] #rural D->R
        out[(i-1)*n_t+13] = γ*u[(3-1)*n + i] #urban A->R
        out[(i-1)*n_t+14] =  γ*u[(3-1)*n + n_s*n + i] #rural A->R
        out[(i-1)*n_t+15] = μ₁*u[(5-1)*n + i] #urban H->death
        out[(i-1)*n_t+16] =  μ₁*u[(5-1)*n + n_s*n + i] #rural H->death
    end
end



function change_matrix(dc,u,p,t,mark)
    for i = 1:n
        dc[(1-1)*n + i,(i-1)*n_t+1] = -1
        dc[(2-1)*n + i,(i-1)*n_t+1] = 1 #change due to urban transmission
        dc[(1-1)*n + n_s*n + i,(i-1)*n_t+2] = -1
        dc[(2-1)*n + n_s*n + i,(i-1)*n_t+2] = 1 #change due to rural transmission
        dc[(2-1)*n + i,(i-1)*n_t+3] = -1
        dc[(3-1)*n + i,(i-1)*n_t+3] = 1
        dc[(7-1)*n + i,(i-1)*n_t+3] = 1#change due to urban E->A
        dc[(2-1)*n + n_s*n + i,(i-1)*n_t+4] = -1
        dc[(3-1)*n + n_s*n + i,(i-1)*n_t+4] = 1
        dc[(7-1)*n + n_s*n + i,(i-1)*n_t+4] = 1 #change due to rural E->A
        dc[(2-1)*n + i,(i-1)*n_t+5] = -1
        dc[(4-1)*n + i,(i-1)*n_t+5] = 1
        dc[(8-1)*n + i,(i-1)*n_t+5] = 1 #change due to urban E->D
        dc[(2-1)*n + n_s*n + i,(i-1)*n_t+6] = -1
        dc[(4-1)*n + n_s*n + i,(i-1)*n_t+6] = 1
        dc[(8-1)*n + n_s*n + i,(i-1)*n_t+6] = 1 #change due to rural E->D
        dc[(4-1)*n + i,(i-1)*n_t+7] = -1
        dc[(5-1)*n + i,(i-1)*n_t+7] = 1 #change due to urban D->H
        dc[(4-1)*n + n_s*n + i,(i-1)*n_t+8] = -1
        dc[(5-1)*n + n_s*n + i,(i-1)*n_t+8] = 1 #change due to rural D->H
        dc[(5-1)*n + i,(i-1)*n_t+9] = -1
        dc[(6-1)*n + i,(i-1)*n_t+9] = 1 #change due to urban H->R
        dc[(5-1)*n + n_s*n + i,(i-1)*n_t+10] = -1
        dc[(6-1)*n + n_s*n + i,(i-1)*n_t+10] = 1 #change due to rural H->R
        dc[(4-1)*n + i,(i-1)*n_t+11] = -1
        dc[(6-1)*n + i,(i-1)*n_t+11] = 1 #change due to urban D->R
        dc[(4-1)*n + n_s*n + i,(i-1)*n_t+12] = -1
        dc[(6-1)*n + n_s*n + i,(i-1)*n_t+12] = 1 #change due to rural D->R
        dc[(3-1)*n + i,(i-1)*n_t+13] = -1
        dc[(6-1)*n + i,(i-1)*n_t+13] = 1 #change due to urban A->R
        dc[(3-1)*n + n_s*n + i,(i-1)*n_t+14] = -1
        dc[(6-1)*n + n_s*n + i,(i-1)*n_t+14] = 1 #change due to rural A->R
        dc[(5-1)*n + i,(i-1)*n_t+15] = -1
        dc[(9-1)*n + i,(i-1)*n_t+15] = 1#change due to urban H->death
        dc[(5-1)*n + n_s*n + i,(i-1)*n_t+16] = -1
        dc[(9-1)*n + n_s*n + i,(i-1)*n_t+16] = 1#change due to rural H->death
    end
end

reg_jumps_forKenyaCoV = RegularJump(rates,change_matrix,dc;constant_c=true)

#This generates the underlying Poisson drivers according to the rate function

function PP_drivers(dN::Vector{Int64},rates,p)
    for i = 1:length(dN)
        if rates[i] >= 0.
            dN[i] = rand(Poisson(p.dt*rates[i]))
        else
            dN[i] = 0
        end
    end
end

#This caps the Poisson processes at causing no more transitions than the group being effected
function max_change(out,u,p::CoVParameters)
    for i = 1:n
        out[(i-1)*n_t+1] = min(out[(i-1)*n_t+1],u[(1-1)*n + i] )#urban transmission
        out[(i-1)*n_t+2] =  min(out[(i-1)*n_t+2],u[(1-1)*n + n_s*n + i] )#rural transmission
        # out[(i-1)*n_t+3] = min(out[(i-1)*n_t+3],u[(2-1)*n + i] )#urban E->A
        # out[(i-1)*n_t+4] =  min(out[(i-1)*n_t+4],u[(2-1)*n + n_s*n + i] )#rural E->A
        # out[(i-1)*n_t+5] = min(out[(i-1)*n_t+5],u[(2-1)*n + i] )#urban E->D
        # out[(i-1)*n_t+6] =  min(out[(i-1)*n_t+6],u[(2-1)*n + n_s*n + i] )#rural E->D
        # out[(i-1)*n_t+7] = min(out[(i-1)*n_t+7],u[(4-1)*n + i] )#urban D->H
        # out[(i-1)*n_t+8] =  min(out[(i-1)*n_t+8],u[(4-1)*n + n_s*n + i] )#rural D->H
        out[(i-1)*n_t+9] = min(out[(i-1)*n_t+9],u[(5-1)*n + i] )#urban H->R
        out[(i-1)*n_t+10] = min(out[(i-1)*n_t+10],u[(5-1)*n + n_s*n + i] )#rural H->R
        # out[(i-1)*n_t+11] = min(out[(i-1)*n_t+11],u[(4-1)*n + i] )#urban D->R
        # out[(i-1)*n_t+12] =  min(out[(i-1)*n_t+12],u[(4-1)*n + n_s*n + i] )#rural D->R
        out[(i-1)*n_t+13] = min(out[(i-1)*n_t+13],u[(3-1)*n + i] )#urban A->R
        out[(i-1)*n_t+14] =  min(out[(i-1)*n_t+14],u[(3-1)*n + n_s*n + i] )#rural A->R
        out[(i-1)*n_t+15] = min(out[(i-1)*n_t+15],u[(5-1)*n + i] )#urban H->death
        out[(i-1)*n_t+16] =  min(out[(i-1)*n_t+16],u[(5-1)*n + n_s*n + i] )#rural H->death
        if out[(i-1)*n_t+3] + out[(i-1)*n_t+5] > u[(2-1)*n + i] #More incubations than actual urban E population
            out[(i-1)*n_t+3] = rand(Binomial(u[(2-1)*n + i],p.δ)) #Binomially distributed the urban incubations between Asympotomatic and symptomatic
            out[(i-1)*n_t+5] = u[(2-1)*n + i] - out[(i-1)*n_t+3]
        end
        if out[(i-1)*n_t+4] + out[(i-1)*n_t+6] > u[(2-1)*n + n_s*n + i]  #More incubations than actual rural E population
            out[(i-1)*n_t+4] = rand(Binomial(u[(2-1)*n + n_s*n + i],p.δ)) #Binomially distributed the rural incubations between Asympotomatic and symptomatic
            out[(i-1)*n_t+6] = u[(2-1)*n + n_s*n + i] - out[(i-1)*n_t+4]
        end
        if out[(i-1)*n_t+7] + out[(i-1)*n_t+11] > u[(4-1)*n + i] # More end of infection events than urban symptomatic infecteds
            out[(i-1)*n_t+7] = rand(Binomial(u[(4-1)*n + i],p.τ/(p.τ + p.γ))) #Binomially distributed the urban end of infections between hospitalisation and recovery
            out[(i-1)*n_t+11] = u[(4-1)*n + i] - out[(i-1)*n_t+7]
        end
        if out[(i-1)*n_t+8] + out[(i-1)*n_t+12] > u[(4-1)*n + n_s*n + i] # More end of infection events than rural symptomatic infecteds
            out[(i-1)*n_t+8] = rand(Binomial(u[(4-1)*n + n_s*n + i],p.τ/(p.τ + p.γ))) #Binomially distributed the rural end of infections between hospitalisation and recovery
            out[(i-1)*n_t+12] = u[(4-1)*n + n_s*n + i] - out[(i-1)*n_t+7]
        end
    end
end

function nonneg_tauleap(du,u,p::CoVParameters,t)
    @unpack dc,dN,poi_rates = p
    rates(poi_rates,u,p,t) #calculate rates of underlying Poisson processes
    PP_drivers(dN,poi_rates,p)#Generate Poisson rvs with rates scaled by time step dt
    max_change(dN,u,p)#Cap the size of the Poisson rvs to maintain non-negativity
    mul!(du,dc,dN)#Calculates the effect on the state in the inplace du vector
    du .+= u #Calculates how the state should change
end

function ode_model(du,u,p::CoVParameters,t)
    @unpack λ_urb,λ_rur,β,γ,σ,δ,τ,μ₁,into_mom,into_nai,global_prev = p
    calculate_infection_rates!(u,p,t)
    for i = 1:n
        du[i,1,1] = (-1)*λ_urb[i]*u[i,1,1]
        du[i,2,1] = λ_urb[i]*u[i,1,1] - σ*u[i,2,1]
        du[i,3,1] = (1-δ)*σ*u[i,2,1] - γ*u[i,3,1]
        du[i,4,1] = δ*σ*u[i,2,1] - γ*u[i,4,1] - τ*u[i,4,1]
        du[i,5,1] = τ*u[i,4,1] - γ*u[i,5,1] - μ₁*u[i,5,1]
        du[i,6,1] = γ*(u[i,3,1] + u[i,4,1] + u[i,5,1])
        du[i,7,1] = (1-δ)*σ*u[i,2,1]
        du[i,8,1] = δ*σ*u[i,2,1]
        du[i,9,1] = μ₁*u[i,5,1]
        du[i,1,2] = (-1)*λ_rur[i]*u[i,1,2]
        du[i,2,2] = λ_rur[i]*u[i,1,2] - σ*u[i,2,2]
        du[i,3,2] = (1-δ)*σ*u[i,2,2] - γ*u[i,3,2]
        du[i,4,2] = δ*σ*u[i,2,2] - γ*u[i,4,2] - τ*u[i,4,2]
        du[i,5,2] = τ*u[i,4,2] - γ*u[i,5,2] - μ₁*u[i,5,2]
        du[i,6,2] = γ*(u[i,3,2] + u[i,4,2] + u[i,5,2])
        du[i,7,2] = (1-δ)*σ*u[i,2,2]
        du[i,8,2] = δ*σ*u[i,2,2]
        du[i,9,2] = μ₁*u[i,5,2]
     end
end
