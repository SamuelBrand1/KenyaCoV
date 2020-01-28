# @with_kw mutable struct CoVParameters
#     β::Float64 = 0.7676
#     γ::Float64 = 1/6
#     σ::Float64 = 1/5
#     ρ::Float64 = 0.1#Proportion of symptomatic/diseased vs non-symptomatic cases
#     τ::Float64 = 1/6. #treatment/isolation rate for symptomatics
#     μ₁::Float64 = 0.01#Excess mortality due to disease
#     ρ::Float64 = 0.01 #spatial coupling
#     T::Matrix{Float64} #transmission matrix
# end
@with_kw mutable struct CoVParameters
    β::Float64 = 0.7676
    γ::Float64 = 1/6
    σ::Float64 = 1/5
    δ::Float64 = 0.1#Proportion of symptomatic/diseased vs non-symptomatic cases
    τ::Float64 = 1/6. #treatment/isolation rate for symptomatics
    μ₁::Float64 = 0.01#Excess mortality due to disease
    ρ::Float64 = 0.01 #spatial coupling
    T::Matrix{Float64} #transmission matrix
end
