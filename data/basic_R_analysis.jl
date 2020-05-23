# Basic R₀ analysis
using Distributions,Plots

 case_dates = ["Feb 15","Feb 16","Feb 17","Feb 18","Feb 19","Feb 20","Feb 21","Feb 22","Feb 23","Feb 24","Feb 25","Feb 26","Feb 27","Feb 28","Feb 29","Mar 01","Mar 02","Mar 03","Mar 04","Mar 05","Mar 06","Mar 07","Mar 08","Mar 09","Mar 10","Mar 11","Mar 12","Mar 13","Mar 14","Mar 15","Mar 16","Mar 17","Mar 18","Mar 19","Mar 20","Mar 21","Mar 22","Mar 23","Mar 24","Mar 25","Mar 26","Mar 27","Mar 28","Mar 29","Mar 30","Mar 31","Apr 01","Apr 02","Apr 03","Apr 04","Apr 05","Apr 06","Apr 07","Apr08",
 "Apr 09","Apr 10","Apr 11","Apr 12","Apr 13","Apr 14","Apr 15","Apr 16","Apr 17","Apr 18","Apr 19","Apr 20","Apr 21","Apr 22","Apr 23","Apr 24","Apr 25","Apr 26","Apr 27","Apr 28","Apr 29","Apr 30","May 01","May 02","May 03","May 04","May 05","May 06","May 07","May 08","May 09","May 10","May 11","May 12","May 13","May 14","May 15","May 16","May 17","May 18","May 19","May 20","May 21"]
new_cases = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,1,3,0,0,0,8,1,9,3,3,0,7,4,8,9,22,29,12,4,16,16,14,7,5,5,2,6,11,8,9,9,12,16,8,11,15,7,17,16,7,12,8,11,10,12,15,24,30,25,45,47,25,14,28,23,28,15,22,21,23,49,57,25,51,66,80]
f = findfirst(case_dates .== "Mar 13")

generation_distribution = Weibull(2.826,)

dates_str = case_dates[f:end]
cases = new_cases[f:end]

scatter(cases,lab = "",xticks = (1:7:length(dates_str),dates_str[1:7:end]))


@with_kw struct CaseData
    data
end
P = CaseData(data=cases)

function (P::CaseData)(θ)
    @unpack r,α,c₀ = θ
    @unpack data = P

    log_likelihood = 0.
    for t in 1:length(data)
        d = data[t]
        μ = c₀*exp(r*t)
        σ² = μ + α*μ^2
        p_negbin = 1 - ((σ² - μ)/σ²)
        r_negbin = μ^2/(σ² - μ)
        log_likelihood += logpdf(NegativeBinomial(r_negbin,p_negbin),d)
    end
    return log_likelihood
end
@time P((r=0.1,α = 0.1,c₀ = .1))
trans = as((r = asℝ₊,α = asℝ₊,c₀ = asℝ₊))



@with_kw struct DynamicHMCPosterior
    "Algorithm for the ODE solver."
    algorithm = Tsit5()
    "An ODE problem definition (`DiffEqBase.DEProblem`)."
    problem
    "Time values at which the simulated path is compared to `data`."
    t
    "Data, as a matrix with each time value in a column."
    data
end

P = DynamicHMCPosterior(problem=prob,t = sol.t,data=data_agg)
function (P::DynamicHMCPosterior)(θ)
    @unpack parameters,α = θ
    @unpack algorithm, problem, data, t = P
    T = eltype(parameters)
    u0 = convert.(T,problem.u0)
    p = parameters
    _saveat = t
    sol = concrete_solve(problem, algorithm, u0, p; saveat = _saveat)

    log_likelihood = 0.
    for i = 1:length(sol.t)
        μ = max(1e-13,sol[1,i])
        p = 1 - ((σ² - μ)/σ²)
        r = μ^2/(σ² - μ)
        log_likelihood += logpdf(NegativeBinomial(r,p),data[i])
    end
    return log_likelihood
end
