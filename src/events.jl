function unnormedCategorical(p̃::AbstractVector)
    normalize!(p̃,1)
    return rand(Categorical(p̃))
end


function urban_transmission(u,p,t)
    @unpack T,β = p
    S_urb = u[:,1,1]
    I_urb = u[:,3,1] .+ u[:,4,1]
    I_rural = u[:,3,2] .+ u[:,4,2]
    Î = T*I_urb + I_rural
    λ_urb = β*T'*(Î ./N̂)
    return sum(S_urb.*λ_urb)
end

function affect_urb_transmission!(integrator)
    @unpack T,β = integrator.p
    u = integrator.u
    S_urb = u[:,1,1]
    I_urb = u[:,3,1] .+ u[:,4,1]
    I_rural = u[:,3,2] .+ u[:,4,2]
    Î = T*I_urb + I_rural
    λ_urb = β*T'*(Î ./N̂)
    # inf_rate = S_urb.*λ_urb
    inf_rate = [S_urb[i]*λ_urb[i] for i = 1:n]
    area = unnormedCategorical(inf_rate)
    u[area,1,1] -= 1
    u[area,2,1] += 1
end
jump_urb_trans = ConstantRateJump(urban_transmission,affect_urb_transmission!)

function rural_transmission(u,p,t)
    @unpack T,β = p
    S_rural = u[:,1,2]
    I_urb = u[:,3,1] .+ u[:,4,1]
    I_rural = u[:,3,2] .+ u[:,4,2]
    Î = T*I_urb + I_rural
    λ_rural = β*(Î./N̂)
    return β*sum(S_rural.*λ_rural)
end
function affect_rural_transmission!(integrator)
    @unpack β,T = integrator.p
    u = integrator.u
    S_rural = u[:,1,2]
    I_urb = u[:,3,1] .+ u[:,4,1]
    I_rural = u[:,3,2] .+ u[:,4,2]
    Î = (T*I_urb) .+ I_rural
    λ_rural = β*(Î./N̂)
    inf_rate = [S_rural[i]*λ_rural[i] for i = 1:n]
    area = unnormedCategorical(inf_rate)
    u[area,1,2] -= 1
    u[area,2,2] += 1
end
jump_rural_trans = ConstantRateJump(rural_transmission,affect_rural_transmission!)

function incubation(u,p,t)
    return p.σ*sum(u[:,2,:])
end

function affect_incubation!(integrator)
    @unpack δ,σ = integrator.p
    u = integrator.u
    E_urb = u[:,2,1]
    E_rural = u[:,2,2]
    incubationrates = σ*(E_urb .+ E_rural)
    choose_disease = 3 + rand(Bernoulli(δ))
    area = unnormedCategorical(incubationrates)
    urb_vs_rural = 1+ rand(Bernoulli(E_rural[area]/(E_urb[area] + E_rural[area])))
    u[area,choose_disease,urb_vs_rural] += 1
    u[area,2,urb_vs_rural] -= 1
    u[area,choose_disease+4,urb_vs_rural] += 1#Add to cumulative cases
end
jump_incubation= ConstantRateJump(incubation,affect_incubation!)

function recovery(u,p,t)
    return p.γ*sum(u[:,3:5,:])
end

function affect_recovery!(integrator)
    @unpack γ = integrator.p
    u = integrator.u
    sub_I = sum(u, dims = 3)[:,3]
    diseased_I = sum(u, dims = 3)[:,4]
    H = sum(u, dims = 3)[:,5]
    area = unnormedCategorical(γ*(sub_I .+ diseased_I .+H))
    group = 2 + unnormedCategorical(γ*[sub_I[area],diseased_I[area],H[area]])
    urb_vs_rural = 1 + rand(Bernoulli(u[area,group,2]/(u[area,group,1]+u[area,group,2])))
    u[area,group,urb_vs_rural] -= 1
    u[area,6,urb_vs_rural] += 1
end
jump_recovery = ConstantRateJump(recovery,affect_recovery!)

function hospitalisation(u,p,t)
    return p.τ*sum(u[:,4,:])
end

function affect_hospitalisation!(integrator)
    @unpack τ = integrator.p
    u = integrator.u
    diseased_I_urb = u[:,4,1]
    diseased_I_rural = u[:,4,2]
    area = unnormedCategorical(τ*(diseased_I_urb .+ diseased_I_rural))
    urb_vs_rural = 1 + rand(Bernoulli(u[area,4,2]/(u[area,4,1]+u[area,4,2])))
    u[area,4,urb_vs_rural] -= 1
    u[area,5,urb_vs_rural] += 1
end
jump_hosp = ConstantRateJump(hospitalisation,affect_hospitalisation!)

function death(u,p,t)
    return p.μ₁*sum(u[:,5,:])
end

function affect_death!(integrator)
    @unpack μ₁ = integrator.p
    u = integrator.u
    H_urb = u[:,5,1]
    H_rural = u[:,5,2]
    area = unnormedCategorical(μ₁*(H_urb .+ H_rural))
    urb_vs_rural = 1 + rand(Bernoulli(u[area,5,2]/(u[area,5,1] +u[area,5,2])))
    u[area,5,urb_vs_rural] -=1
    u[area,9,urb_vs_rural] +=1
end
jump_death = ConstantRateJump(death,affect_death!)
