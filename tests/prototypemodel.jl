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

M_Kenya_ho_s = SMatrix{17,17}(M_Kenya_ho)
M_Kenya_other_s = SMatrix{17,17}(M_Kenya_other)
M_Kenya_school_s = SMatrix{17,17}(M_Kenya_school)


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
A = contact_data
extrap = CubicSplineInterpolation(xs, A, extrapolation_bc = Line())
extrap(700)
function build_ct_and_grad_ct(xs,A)
        extrap = CubicSplineInterpolation(xs, A, extrapolation_bc = Line())
        return t -> extrap(t),t -> Interpolations.gradient(extrap,t)[1]
end

ct,∇ct = build_ct_and_grad_ct(xs,A)


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

function kenyacov_jac(J,u,p,t)
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

        λ = β₀.*((M_Kenya_ho .+ ct(t).*M_Kenya_other .+ 0*M_Kenya_school)*(ϵ.*A .+ P .+ M .+V))./N
        # ∂(-Sₐ λₐ)/∂Sb = δ_ab λ_ₐ
        for a = 1:17
                J[a,a] = λ[a]
        end

        for a = 1:17,b = 1:17
                # ∂(-Sₐ λₐ)/∂Ab = -S_a ∂(λₐ)/∂Ab
                i = a
                j = b + 2*17
                J[i,j] = -S[a]*β₀.*(M_Kenya_ho[a,b] .+ ct(t).*M_Kenya_other[a,b] .+ 0*M_Kenya_school[a,b])*(ϵ)./N[a]
        end
end
length(u0[:])
J = spzeros(136,136)
@time kenyacov_jac(J,u0,parameters,0.)
## Fastest get incidence function

rem(17,17)
floor(Int64,10/17)


I₀ = fill(10.,17)
u0 = zeros(17,8)
u0[:,1] = N
u0[:,2] = I₀
tspan = (0.,365*1)
p = copy(parameters)
p[1] = p[1]/3
p[1] = 0.
@benchmark kenyacov_ode(similar(u0),u0,p,0.)


prob = ODEProblem(kenyacov_ode,u0,tspan,p)
sol = solve(prob,Tsit5(),p = parameters)
@benchmark sol = solve(prob,Rosenbrock23(),p = [3,N])
@benchmark sol = solve(prob,Tsit5(),p = p)
@benchmark sol = solve(prob,Tsit5(),p = p,saveat = 1.,u0=u0)
# @benchmark sol = solve(prob,Rosenbrock23(),p = parameters,saveat = 1.)

function get_incidence_time_array(sol)
        Matrix(VectorOfArray(diff([u[:,8] for u in sol.u]))')
end
inc = get_incidence_time_array(sol)
plot(inc)
bar(sol.u[end][:,8]./N,lab = "")

#Using model toolkit functionality
@parameters t
@register ct(t)
sys = modelingtoolkitize(prob)
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
du = simplify.(kenyacov_ode_n(u,p,t))
ModelingToolkit.build_function(du,u,p,t)
fastf = eval(ModelingToolkit.build_function(du,u,p,t,
           parallel=ModelingToolkit.MultithreadedForm())[2])


function kenyacov_ode_n(u,p,t)
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
        λ = β₀.*((M_Kenya_ho .+ M_Kenya_other .+ M_Kenya_school)*(ϵ.*A .+ P .+ M .+V))./N
        #RHS of vector field (in-place calculation)
        du = similar(u)
        du[:,1] .= -S.*λ
        du[:,2] .= S.*λ .- α.*E
        du[:,3] .= α.*E.*(1 .- δ) .- γA.*A
        du[:,4] .= α.*E.*δ .- αP.*P
        du[:,5] .= αP.*P.*(1 .- υ) .- γM.*M
        du[:,6] .= αP.*P.*υ .- γV.*V
        du[:,7] .= γA.*A .+ γM.*M
        du[:,8] .= S.*λ

        return du
end
