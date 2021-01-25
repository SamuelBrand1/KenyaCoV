push!(LOAD_PATH,joinpath(homedir(),"Github/KenyaCoV/src"))

using Revise
import KenyaCoV
using DifferentialEquations,ModelingToolkit,Latexify,SparseArrays,StaticArrays,LinearAlgebra,Plots
using DelimitedFiles,JLD2,BenchmarkTools,DataFrames,CSV,RecursiveArrayTools,Interpolations
using StaticArrays,SparseArrays

## Rescale contact matrices
JLD2.@load("data/agemixingmatrix_Kenya_all_types.jld2")
Ronescaler = 1/Real(eigvals(M_Kenya)[end])
M_Kenya_ho .= M_Kenya_ho.*Ronescaler
M_Kenya_other .= (M_Kenya_other .+ M_Kenya_work).*Ronescaler
M_Kenya_school .= M_Kenya_school.*Ronescaler


# JLD2.@load("notebooks/d_log_priors.jld2")
# @save("data/d_log_priors.jld2",d_log_priors,prior_param_names)
JLD2.@load("data/d_log_priors.jld2",d_log_priors,prior_param_names)
# prior_param_names = [:β₀,:α,:ϵ,:αP,:γA,:γM]
# append!(prior_param_names,[Symbol("δ_"*string(n)) for n = 1:16])
# append!(prior_param_names,[Symbol("σ_"*string(n)) for n = 1:15])

parameters = exp.(rand(d_log_priors))
β₀,α,ϵ,αP,γA,γM = parameters[1:6]
δ = vcat(parameters[7:(6+16)],parameters[6+16])
σ = vcat(parameters[(6+16+1):end],[1.,1.])


df_pop = DataFrame(CSV.File("data/2019_census_age_pyramids_counties.csv"))
N = [df_pop[df_pop.county .== "Nairobi",a+1][1] for a = 1:17]
bar(N,xticks = (1:17,vcat([string(a*5)*"-"*string(a*5 + 4) for a = 0:15],"80+")))

#Create function for contact rate change
@load("data/projected_contact_data_20112020.jld2")
contact_data = projected_contactrate_nairobi.contactrate
xs = 1:294

function build_ct_and_grad_ct(xs,A)
        extrap = CubicSplineInterpolation(xs, A, extrapolation_bc = Line())
        return t -> extrap(t),t -> Interpolations.gradient(extrap,t)[1]
end

function ct_from_data(t,A)
        extrap = CubicSplineInterpolation(xs, A, extrapolation_bc = Line())
        return extrap(t)
end
function ∂ct∂t_from_data(t,A)
        extrap = CubicSplineInterpolation(xs, A, extrapolation_bc = Line())
        return Interpolations.gradient(extrap,t)[1]
end
ct_from_data(700,contact_data)
∂ct∂t_from_data(700,contact_data)


ct,∇ct = build_ct_and_grad_ct(xs,contact_data)


function kenyacov_ode(du,u,p,t)
        #Gather parameters
        β₀,α,ϵ,αP,γA,γM = p[1:6]
        δ = vcat(p[7:(6+16)],p[6+16])
        σ = vcat(p[(6+16+1):end],[1.,1.])
        υ = zeros(17)
        γV = γM
        #Gather state
        S = @view u[:,1]
        E = @view u[:,2]
        A = @view u[:,3]
        P = @view u[:,4]
        M = @view u[:,5]
        V = @view u[:,6]
        R = @view u[:,7]
        cum_incidence = @view u[:,8]
        #Force of infection
        λ = β₀.*((M_Kenya_ho .+ ct(t).*M_Kenya_other .+ 0*M_Kenya_school)*(ϵ.*A .+ P .+ M .+V))./N
        #RHS of vector field (in-place calculation)
        du[:,1] .= -S.*λ
        du[:,2] .= S.*λ .- α.*E
        du[:,3] .= α.*E.*(1 .- δ) .- γA.*A
        du[:,4] .= α.*E.*δ .- αP.*P
        du[:,5] .= αP.*P.*(1 .- υ) .- γM.*M
        du[:,6] .= αP.*P.*υ .- γV.*V
        du[:,7] .= γA.*A .+ γM.*M
        du[:,8] .= S.*λ

        return nothing
end
function kenyacov_ode_with_contactrates(du,u,p,t)
        #Gather parameters
        β₀,α,ϵ,αP,γA,γM = p[1:6]
        δ = vcat(p[7:(6+16)],p[6+16])
        σ = vcat(p[(6+16+1):37],[1.,1.])
        contactrate = @view p[38:end]

        υ = zeros(17)
        γV = γM
        #Gather state
        S = @view u[:,1]
        E = @view u[:,2]
        A = @view u[:,3]
        P = @view u[:,4]
        M = @view u[:,5]
        V = @view u[:,6]
        R = @view u[:,7]
        cum_incidence = @view u[:,8]
        #Force of infection
        λ = β₀.*((M_Kenya_ho .+ ct_from_data(t,contactrate).*M_Kenya_other .+ 0*M_Kenya_school)*(ϵ.*A .+ P .+ M .+V))./N
        #RHS of vector field (in-place calculation)
        du[:,1] .= -S.*λ
        du[:,2] .= S.*λ .- α.*E
        du[:,3] .= α.*E.*(1 .- δ) .- γA.*A
        du[:,4] .= α.*E.*δ .- αP.*P
        du[:,5] .= αP.*P.*(1 .- υ) .- γM.*M
        du[:,6] .= αP.*P.*υ .- γV.*V
        du[:,7] .= γA.*A .+ γM.*M
        du[:,8] .= S.*λ

        return nothing
end


function kenyacov_ode_tgrad(dT,u,p,t)
        #Gather parameters
        β₀,α,ϵ,αP,γA,γM = p[1:6]
        δ = vcat(p[7:(6+16)],p[6+16])
        σ = vcat(p[(6+16+1):end],[1.,1.])
        υ = zeros(17)
        γV = γM
        #Gather state
        S = @view u[:,1]
        E = @view u[:,2]
        A = @view u[:,3]
        P = @view u[:,4]
        M = @view u[:,5]
        V = @view u[:,6]
        R = @view u[:,7]
        cum_incidence = @view u[:,8]
        #Force of infection
        ∂λ∂t = β₀.*(∇ct(t).*M_Kenya_other)*(ϵ.*A .+ P .+ M .+V)./N
        #RHS of vector field (in-place calculation)
        dT[:,1] .= -S.*∂λ∂t
        du[:,2] .= S.*∂λ∂t
        du[:,3] .= 0.
        du[:,4] .= 0.
        du[:,5] .= 0.
        du[:,6] .= 0.
        du[:,7] .= 0.
        du[:,8] .= S.*∂λ∂t

        return nothing
end
#Method for caluclating jacobians without the time function problem

function kenyacov_ode_notrans(du,u,p,t)
        #Gather parameters
        β₀,α,ϵ,αP,γA,γM = p[1:6]
        δ = vcat(p[7:(6+16)],p[6+16])
        σ = vcat(p[(6+16+1):end],[1.,1.])
        υ = zeros(17)
        γV = γM
        #Gather state
        S = @view u[:,1]
        E = @view u[:,2]
        A = @view u[:,3]
        P = @view u[:,4]
        M = @view u[:,5]
        V = @view u[:,6]
        R = @view u[:,7]
        cum_incidence = @view u[:,8]
        #Force of infection
        # λ = β₀.*((M_Kenya_ho .+ ct(t).*M_Kenya_other .+ 0*M_Kenya_school)*(ϵ.*A .+ P .+ M .+V))./N
        #RHS of vector field (in-place calculation)
        du[:,1] .= 0
        du[:,2] .= 0 .- α.*E
        du[:,3] .= α.*E.*(1 .- δ) .- γA.*A
        du[:,4] .= α.*E.*δ .- αP.*P
        du[:,5] .= αP.*P.*(1 .- υ) .- γM.*M
        du[:,6] .= αP.*P.*υ .- γV.*V
        du[:,7] .= γA.*A .+ γM.*M
        du[:,8] .= 0

        return nothing
end

function kenyacov_ode_hometransonly(du,u,p,t)
        #Gather parameters
        β₀,α,ϵ,αP,γA,γM = p[1:6]
        δ = vcat(p[7:(6+16)],p[6+16])
        σ = vcat(p[(6+16+1):end],[1.,1.])
        υ = zeros(17)
        γV = γM
        #Gather state
        S = @view u[:,1]
        E = @view u[:,2]
        A = @view u[:,3]
        P = @view u[:,4]
        M = @view u[:,5]
        V = @view u[:,6]
        R = @view u[:,7]
        cum_incidence = @view u[:,8]
        #Force of infection
        λ = β₀.*((M_Kenya_ho)*(ϵ.*A .+ P .+ M .+V))./N
        #RHS of vector field (in-place calculation)
        du[:,1] .= -S.*λ
        du[:,2] .= S.*λ
        du[:,3] .= 0.
        du[:,4] .= 0.
        du[:,5] .= 0.
        du[:,6] .= 0.
        du[:,7] .= 0.
        du[:,8] .= S.*λ

        return nothing
end

function kenyacov_ode_othertransonly(du,u,p,t)
        #Gather parameters
        β₀,α,ϵ,αP,γA,γM = p[1:6]
        δ = vcat(p[7:(6+16)],p[6+16])
        σ = vcat(p[(6+16+1):end],[1.,1.])
        υ = zeros(17)
        γV = γM
        #Gather state
        S = @view u[:,1]
        E = @view u[:,2]
        A = @view u[:,3]
        P = @view u[:,4]
        M = @view u[:,5]
        V = @view u[:,6]
        R = @view u[:,7]
        cum_incidence = @view u[:,8]
        #Force of infection
        λ = β₀.*((M_Kenya_other)*(ϵ.*A .+ P .+ M .+V))./N
        #RHS of vector field (in-place calculation)
        du[:,1] .= -S.*λ
        du[:,2] .= S.*λ
        du[:,3] .= 0.
        du[:,4] .= 0.
        du[:,5] .= 0.
        du[:,6] .= 0.
        du[:,7] .= 0.
        du[:,8] .= S.*λ

        return nothing
end


function kenyacov_ode_schooltransonly(du,u,p,t)
        #Gather parameters
        β₀,α,ϵ,αP,γA,γM = p[1:6]
        δ = vcat(p[7:(6+16)],p[6+16])
        σ = vcat(p[(6+16+1):end],[1.,1.])
        υ = zeros(17)
        γV = γM
        #Gather state
        S = @view u[:,1]
        E = @view u[:,2]
        A = @view u[:,3]
        P = @view u[:,4]
        M = @view u[:,5]
        V = @view u[:,6]
        R = @view u[:,7]
        cum_incidence = @view u[:,8]
        #Force of infection
        λ = β₀.*((M_Kenya_school)*(ϵ.*A .+ P .+ M .+V))./N
        #RHS of vector field (in-place calculation)
        du[:,1] .= -S.*λ
        du[:,2] .= S.*λ
        du[:,3] .= 0.
        du[:,4] .= 0.
        du[:,5] .= 0.
        du[:,6] .= 0.
        du[:,7] .= 0.
        du[:,8] .= S.*λ

        return nothing
end


## Set up initial conditions and time span

I₀ = fill(10.,17)
u0 = zeros(17,8)
u0[:,1] = N
u0[:,2] = I₀
tspan = (0.,365*1)
p = copy(parameters)
p[1] = p[1]/3
p2 = copy(vcat(p,contact_data))

## Underlying sparse arrays for jac

prob_notrans = ODEProblem(kenyacov_ode_notrans,u0,tspan,p)
prob_hometransonly = ODEProblem(kenyacov_ode_hometransonly,u0,tspan,p)
prob_othertransonly = ODEProblem(kenyacov_ode_othertransonly,u0,tspan,p)
prob_schooltransonly = ODEProblem(kenyacov_ode_schooltransonly,u0,tspan,p)

sys_notrans = modelingtoolkitize(prob_notrans)
sys_hometransonly = modelingtoolkitize(prob_hometransonly)
sys_othertransonly = modelingtoolkitize(prob_othertransonly)
sys_schooltransonly = modelingtoolkitize(prob_schooltransonly)

write("src/jac_notrans.jl", string(ModelingToolkit.generate_jacobian(sys_notrans)[2]));
write("src/jac_hometransonly.jl", string(ModelingToolkit.generate_jacobian(sys_hometransonly)[2]));
write("src/jac_othertransonly.jl", string(ModelingToolkit.generate_jacobian(sys_othertransonly)[2]));
write("src/jac_schooltransonly.jl", string(ModelingToolkit.generate_jacobian(sys_schooltransonly)[2]));

jac_notrans =  eval(include("../src/jac_notrans.jl"))
jac_hometransonly =  eval(include("../src/jac_hometransonly.jl"))
jac_othertransonly =  eval(include("../src/jac_othertransonly.jl"))
jac_schooltransonly =  eval(include("../src/jac_schooltransonly.jl"))

u0[:]
J = zeros(136,136)
@time jac(J,u0,p,1.);



function jac(J,u,p,t)
        J_1 = similar(J)
        J_2 = similar(J)
        J_3 = similar(J)
        J_4 = similar(J)
        jac_notrans(J_1,u,p,t)
        jac_hometransonly(J_2,u,p,t)
        jac_othertransonly(J_3,u,p,t)
        jac_schooltransonly(J_4,u,p,t)
        @. J = J_1 + J_2 + ct(t)*J_3 + J_4
        return nothing
end




kenyacov_ode_jac = ODEFunction(kenyacov_ode;jac=jac)
kenyacov_ode_opt_tgrad = ODEFunction(kenyacov_ode;tgrad = kenyacov_ode_tgrad)

prob_nojac = ODEProblem(kenyacov_ode,u0,tspan,p)
prob_tgrad = ODEProblem(kenyacov_ode_opt_tgrad,u0,tspan,p)
prob_withcts = ODEProblem(kenyacov_ode_with_contactrates,u0,tspan,p2)

prob_jac = ODEProblem(kenyacov_ode_opt,u0,tspan,p)

# Time test benchmarks
@benchmark sol = solve(prob_nojac,Rosenbrock23(),p = p,saveat = 1.,u0=u0)
@benchmark sol = solve(prob_jac,Rosenbrock23(),p = p,saveat = 1.,u0=u0)
sol = solve(prob_jac,Rosenbrock23(),p = p,saveat = 1.,u0=u0)
@benchmark sol = solve(prob_nojac,BS3(),p = p,saveat = 1.,u0=u0)
@benchmark sol = solve(prob_withcts,BS3(),p = p2,saveat = 1.,u0=u0)

@benchmark sol = solve(prob_tgrad,BS3(),p = p,saveat = 1.,u0=u0)
@time sol = solve(prob_tgrad,BS3(),p = p,saveat = 1.,u0=u0)
@benchmark sol = solve(prob_nojac,Tsit5(),p = p,saveat = 1.,u0=u0)




# @benchmark sol = solve(prob,Rosenbrock23(),p = parameters,saveat = 1.)

function get_incidence_time_array(sol)
        Matrix(VectorOfArray(diff([u[:,8] for u in sol.u]))')
end
inc = get_incidence_time_array(sol)
plot!(inc)
bar!(sol.u[end][:,8]./N,lab = "")

## Using model toolkit functionality
@parameters t
@register ct(t)
sys = modelingtoolkitize(prob_nojac)
sys1 =


jac = eval(ModelingToolkit.generate_jacobian(sys,parallel=ModelingToolkit.MultithreadedForm())[2])
function tgrad(dT,u,p,t)

sp_pattern =ModelingToolkit.jacobian_sparsity(sys)

kenyacov_ode_sym = ODEFunction(kenyacov_ode,jac=jac,tgrad=tgrad,sparsity = sp_pattern)
prob_sym = ODEProblem(kenyacov_ode_sym,u0,tspan,p)

@benchmark sol2 = solve(prob_sym,Rosenbrock23(),p = parameters,saveat = 1.)
sys


@variables u[1:17,1:8]
@parameters t
@parameters p[1:37]

@variables x[1:100]
@variables A[1:100]

extrap = CubicSplineInterpolation(x, A, extrapolation_bc = Line())
extrap_func(a,A) = CubicSplineInterpolation(x, A, extrapolation_bc = Line())
ct_func(t) = extrap_func(a,A)(t)
@register ct_func(t)

du = simplify.(kenyacov_ode_notinplace(u,p,t))
ModelingToolkit.build_function(du,u,p,t)
fastf = eval(ModelingToolkit.build_function(du,u,p,t,
           parallel=ModelingToolkit.MultithreadedForm())[2])
