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
dc = zeros(n*n_s*2,n_t*n)


function calculate_infection_rates!(u,p::CoVParameters)
    I_urb_A = @view u[((3-1)*n + 1):((3-1)*n + n)]
    I_urb_D = @view u[((4-1)*n + 1):((4-1)*n + n)]
    I_rur_A = @view u[((3-1)*n + n_s*n + 1):((3-1)*n + n_s*n + n)]
    I_rur_D = @view u[((4-1)*n + n_s*n + 1):((4-1)*n + n_s*n + n)]
    mul!(p.Î,p.T,I_urb_A .+ I_urb_D  )
    p.Î .+=  I_rur_A .+ I_rur_D
    mul!(p.λ_urb,p.T',p.β .*(p.Î ./p.N̂))
    p.λ_rur .= p.β .*(p.Î ./p.N̂)
    return nothing
end

function rates(out,u,p::CoVParameters,t)
    @unpack λ_urb,λ_rur,β,γ,σ,δ,τ,μ₁,ϵ_mom,ϵ_nai = p
    calculate_infection_rates!(u,p)
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
    out[(29-1)*n_t+1] = ϵ_mom*u[(1-1)*n + 29]
    out[(31-1)*n_t+1] = ϵ_nai*u[(1-1)*n + 31]
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
