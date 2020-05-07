push!(LOAD_PATH, "/Users/Sam/GitHub/KenyaCoV/src")
using DataFrames,CSV,MAT,Statistics,LinearAlgebra,Optim,Plots,JLD2,RData,Distances,DelimitedFiles
gr()


"""
Load conversion matrices, movement data and population data
"""
c_to_rr = matread("data/conversion_matrix_c_rr.mat")["conversion_matrix_c_to_rr"]
rr_to_c =  matread("data/conversion_matrix.mat")["conversion_matrix"]
@load "data/mv_matrix.jld2" mv_matrix;#In from-to format
county_populations = CSV.read("data/2019_census_age_pyramids_counties.csv")
N_pop = zeros(Int64,47,17)
for i = 1:47, a = 1:17
    N_pop[i,a] = county_populations[i,a+1]
end
mobile_pop_by_county = Vector(sum(N_pop[:,5:11],dims = 2)[:])
bar(mobile_pop_by_county,xticks = (1:47,county_populations.county))

"""
Get location data
"""

kenya_09_tbl = CSV.read("data/2009_National_Estimates_of_Rural_and_Urban_Populations_by_County.csv")
county_centre_points = []
for i = 1:47
    global county_centre_points
    x = join(split(kenya_09_tbl.Location[i],"("))
    x = join(split(x,")"))
    y = split(x,",")
    z = parse(Float64,y[2])*111.32,parse(Float64,y[1])*110.57
    push!(county_centre_points,z)
end

"""
Get distances between counties
"""

dist_matrix = zeros(47,47) #County level distances
for i = 1:47,j=1:47
    dist_matrix[i,j] = sqrt((county_centre_points[i][1] - county_centre_points[j][1])^2 +
                        (county_centre_points[i][2] - county_centre_points[j][2])^2)
end
"""
Convert movements into time spent in areas
"""
movements_per_person = zeros(20,20)
for i = 1:20,j=1:20
    movements_per_person[i,j] = mv_matrix[j,i]/1000#change to from col to row format for ease of multiplication with column vectors
end
heatmap(movements_per_person)
ρ = [sum(movements_per_person,dims = 1)[j]*5/30 for j = 1:20]
P_dest = zeros(20,20)#probability distribution of destination
for j = 1:20
    P_dest[:,j] = normalize(movements_per_person[:,j],1)
end
T = zeros(20,20) # combined location density matrix
for i = 1:20,j=1:20
    if i == j
        T[i,j] = 1-ρ[i]
    else
        T[i,j] = ρ[j]*P_dest[i,j]
    end
end
heatmap(T,clims = (0.,0.1))


"""

"""

function prediction_where_travelled_to(dist_matrix,N,α)
    prediction = similar(dist_matrix)
    d1,d2 = size(dist_matrix)
    for i = 1:d1,j = 1:d2
        if i != j
            prediction[i,j] = N[i]/(1 + dist_matrix[i,j]^α)
        else
            prediction[i,j] = 0.
        end
    end
    for j = 1:d2
        prediction[:,j] = normalize(prediction[:,j],1)
    end
    return prediction
end



function aggregate_movements_to_risk_regions(dist_matrix,N,α)
    prediction_county = prediction_where_travelled_to(dist_matrix,N,α)
    return c_to_rr*prediction_county*rr_to_c
end

function rel_entropy(q,p)#rel entropy of p wrt q
    sum = 0.
    for i=1:length(q)
        if p[i] > 0. && q[i] >0.
            sum += q[i]*log(q[i]/p[i])
        end
    end
    return sum
end

function gravity_model_rel_entropy(mv_data,dist_matrix,N,α)
    d1,d2 = size(mv_data) #<--- base distributions
    P = prediction_where_travelled_to(dist_matrix,N,α) #<--- target distributions
    P̂ = aggregate_movements_to_risk_regions(dist_matrix,N,α)
    sum_rel_entropy = 0.
    for j = 1:d1
        sum_rel_entropy += rel_entropy(P̂[:,j],mv_data[:,j])
    end
    return sum_rel_entropy
end
gravity_model_rel_entropy(normed_movements,dist_matrix,mobile_pop_by_county,2)

"""
Optimise the distance exponent
"""
function opt_func(θ)
    α = θ[1]
    error = gravity_model_rel_entropy(normed_movements,dist_matrix,mobile_pop_by_county,α)
end
fit = optimize(opt_func,[2.],Newton())
α_opt = Optim.minimizer(fit)
rel_min = Optim.minimum(fit)
P_opt = prediction_where_travelled_to(dist_matrix,mobile_pop_by_county,α_opt[1])

@save "data/optimal_movement_matrix.jld2" P_opt
P[28,30] - NtoMmoves
P[30,28] - MtoNmoves
"""
Plot the total relative entropy
"""
α_range = 0.:0.01:5.
rel_entropy_range = [opt_func([α]) for α in α_range]
plt = plot(α_range,rel_entropy_range,
    lab = "Total relative entropy",lw =3,
    xlabel = "distance exponent",
    title = "Fit of gravity model to aggregated movement data",
    ylabel = "Total rel. entropy")
scatter!(plt,[α_opt],[rel_min],ms=5, lab = "minimum at alpha = $(round(α_opt[1],digits = 2))")

savefig(plt,"data/rel_entropy_fit_for_movements.png")
"""
Estimate time spent away by county
"""
ρ_county = ones(47)
for i = 1:47
    rr_distribution = c_to_rr[:,i]
    ρ_county[i] = sum(ρ.*rr_distribution) #Mean value of rho over the likely risk region the county deweller is in
end

"""
Calculate optimal transition matrix
"""
T_opt = similar(P_opt)
d1,d2 = size(T_opt)
for i = 1:d1,j=1:d2
    if i == j
        T_opt[i,j] = 1-ρ_county[j]
    else
        T_opt[i,j] = ρ_county[j]*P_opt[i,j]
    end
end
heatmap(T_opt,clims = (0.,0.2))
@save "data/optimal_transition_matrix.jld2" T_opt ρ_county
