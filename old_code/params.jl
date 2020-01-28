#Some constant parameters
γ1 = 1/9
γ2 = 0.25
inf_2 = 0.5
sus_2 = 0.75
Δt = 2.0
ν = 1/180.

# Time varying mortality and crude birth rates based on 2014 demographic goals being uniformly successful
μ0 = 1/(57*365.25)
μ(t::Float64) = (1-(t/(21*365.25))) + (54./64.)*(t/(21*365.25))
BRScale2014 = 4.6*(1-(5/21)) + 2.6*(5/21)
BR(t::Float64) = (0.5/(57*365.25))*((4.6/BRScale2014)*(1-(t/(21*365.25))) + (2.6/BRScale2014)*(t/(21*365.25)) )

type VariableRSVParameters
    ρ::Float64
    β::Float64
    ξ::Float64
end

P = VariableRSVParameters(0.05,2.,0.15)

InitialKenyaEpidemic = zeros(Int64,n,6)
for i = 1:n
    InitialKenyaEpidemic[i,1] = round(Int64,KenyaTbl[:Total][i]*0.01)
    InitialKenyaEpidemic[i,4] = round(Int64,KenyaTbl[:Total][i]*0.99)
end
