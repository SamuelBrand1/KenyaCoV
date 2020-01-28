function unnormedCategorical(p̃::AbstractVector)
    normalize!(p̃,1)
    return rand(Categorical(p̃))
end


function urban_to_urban_transmission(u,p,t)
    @unpack T,β = p
    S_urban = u[:,1,1]
    I_urban = u[:,3,1] .+ u[:,4,1]
    return sum(β*S_urban.*(T*I_urban) ./N_urban)
end

function affect_uu_transmission!(integrator)
    @unpack T,β = integrator.p
    u = integrator.u
    S_urban = u[:,1,1]
    I_urban = u[:,3,1] .+ u[:,4,1]
    inf_rate = β*S_urban.*(T*I_urban) ./N_urban
    area = unnormedCategorical(inf_rate)
    u[area,1,1] -= 1
    u[area,2,1] += 1
end
jump_uu_trans = ConstantRateJump(urban_to_urban_transmission,affect_uu_transmission!)

function within_county_transmission(u,p,t)
    @unpack ρ,β = p
    S_area = u[:,1,1] + u[:,1,2]
    I_area = (1-ρ)*(u[:,3,1] + u[:,4,1]) + u[:,3,2] + u[:,4,2]
    return β*sum(S_area.*I_area ./(N_rural .+ 1))
end
function affect_wct_transmission!(integrator)
    @unpack ρ,β = integrator.p
    u = integrator.u
    S_area = u[:,1,1] + u[:,1,2]
    I_area = (1-ρ)*(u[:,3,1] + u[:,4,1]) + u[:,3,2] + u[:,4,2]
    inf_rate = β*S_area.*I_area ./(N_rural + 1)
    area = unnormedCategorical(inf_rate)
    urb_vs_rural = 1 + rand(Bernoulli(u[area,1,2]/(u[area,1,1] +u[area,1,2])))
    u[area,1,urb_vs_rural] -= 1
    u[area,2,urb_vs_rural] += 1
end
jump_wct_trans = ConstantRateJump(within_county_transmission,affect_wct_transmission!)

function incubation(u,p,t)
    return p.σ*sum(u[:,2,:])
end

function affect_incubation!(integrator)
    @unpack δ,σ = integrator.p
    u = integrator.u
    E_urban = u[:,2,1]
    E_rural = u[:,2,2]
    incubationrates = σ*(E_urban .+ E_rural)
    choose_disease = 3 + rand(Bernoulli(δ))
    area = unnormedCategorical(incubationrates)
    urb_vs_rural = 1+ rand(Bernoulli(E_rural[area]/(E_urban[area] + E_rural[area])))
    u[area,choose_disease,urb_vs_rural] += 1
    u[area,2,urb_vs_rural] -= 1
    u[area,choose_disease+4,urb_vs_rural] += 1#Add to cumulative cases
end
jump_incubation= ConstantRateJump(incubation,affect_incubation!)

function recovery(u,p,t)
    return p.γ*sum(u[:,3:4,:])
end

function affect_recovery!(integrator)
    @unpack γ = integrator.p
    u = integrator.u
    sub_I = sum(u, dims = 3)[:,3]
    diseased_I = sum(u, dims = 3)[:,4]
    area = unnormedCategorical(γ*(sub_I .+ diseased_I))
    sub_vs_dis = 3 + rand(Bernoulli(diseased_I[area]/(sub_I[area] + diseased_I[area])))
    urb_vs_rural = 1+ rand(Bernoulli(u[area,sub_vs_dis,2]/(u[area,sub_vs_dis,1]+u[area,sub_vs_dis,2])))
    u[area,sub_vs_dis,urb_vs_rural] -= 1
    u[area,6,urb_vs_rural] += 1
end
jump_recovery = ConstantRateJump(recovery,affect_recovery!)

function hospitalisation(u,p,t)
    return p.τ*sum(u[:,4,:])
end

function affect_hospitalisation!(integrator)
    @unpack τ = integrator.p
    u = integrator.u
    diseased_I_urban = u[:,4,1]
    diseased_I_rural = u[:,4,2]
    area = unnormedCategorical(τ*(diseased_I_urban .+ diseased_I_rural))
    urb_vs_rural = 1+ rand(Bernoulli(u[area,4,2]/(u[area,4,1]+u[area,4,2])))
    u[area,4,urb_vs_rural] -= 1
    u[area,5,urb_vs_rural] += 1
end
jump_hosp = ConstantRateJump(hospitalisation,affect_hospitalisation!)

function leave_hospital(u,p,t)
    return (p.γ + p.μ)*sum(u[:,5,:])
end

function affect_leave_hospital!(integrator)
    @unpack γ,μ = integrator.p
    u = integrator.u
    diseased_H_urban = u[:,5,1]
    diseased_H_rural = u[:,5,2]
    area = unnormedCategorical((γ+μ)*(diseased_H_urban .+ diseased_H_rural))
    cure_vs_dead = rand(Bernoulli(γ/(γ+μ)))
    urb_vs_rural = 1+ rand(Bernoulli(u[area,5,2]/(u[area,5,1] +u[area,5,2])))
    if cure_vs_dead
        u[area,5,urb_vs_rural] -=1
        u[area,6,urb_vs_rural] +=1
    else
        u[area,5,urb_vs_rural] -=1
        u[area,9,urb_vs_rural] +=1
    end
end
jump_leave_hosp = ConstantRateJump(leave_hospital,affect_leave_hospital!)
