#Model for transmission between Counties in Kenya
using DataFrames,Distributions,Plots#,MAT
gr()

include("ReadInCountyData.jl");
include("GravityModel.jl");
include("params.jl");
include("StochasticModelStep.jl");


#Choose some new parameters
P.ρ = 0.01
P.β = 1.5
P.ξ = 0.5
# Equilibriate
u = InitialKenyaEpidemic;
t = -730.0*Δt
TimeSeries = zeros(Int64,n,6,730)
TimeIndex = zeros(730)
#Add some infecteds
u[30,4] -= 20;
u[30,5] += 20;
#Equilibriate model
for timestep = 1:730
    TimeIndex[timestep] = t
    TimeSeries[:,:,timestep] = u
    u = EquilibriumEpidemicEvents!(t,u,P)
    t += Δt
end



CNum = 30 #Nairobi
# plot(Δt*collect(1:730)/365.25,InfTimeSeries[CNum,:,1]+InfTimeSeries[CNum,:,2],lab = KenyaTbl[:County][CNum],linewidth = 3)
plot(Δt*collect(1:730)/365.25,TimeSeries[CNum,2,:] + TimeSeries[CNum,5,:],lab = KenyaTbl[:County][CNum],linewidth = 3,title = "Number of infecteds")

CNum = 14 #Kilifi
plot!(Δt*collect(1:730)/365.25,TimeSeries[CNum,2,:]+ TimeSeries[CNum,5,:],lab = KenyaTbl[:County][CNum])

CNum = 28 #Mombasa
plot!(Δt*collect(1:730)/365.25,TimeSeries[CNum,2,:]+ TimeSeries[CNum,5,:],lab = KenyaTbl[:County][CNum])

CNum = 17 #Kisumu
plot!(Δt*collect(1:730)/365.25,TimeSeries[CNum,2,:]+ TimeSeries[CNum,5,:],lab = KenyaTbl[:County][CNum])

# file = matopen("KenyaEpidemic.mat","w")
# write(file,"T",TimeIndex)
# write(file,"Epidemic",TimeSeries)
# close(file)

savefig("BasicTimeSeries.pdf")
