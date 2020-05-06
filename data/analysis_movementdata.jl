push!(LOAD_PATH, "/Users/Sam/GitHub/KenyaCoV/src")
using DataFrames,CSV,MAT,Statistics,LinearAlgebra,Optim,Plots,JLD2,RData,Distances
gr()



movement_data =  RData.load("data/movement_matrix_2020-03.rda")
movements = movement_data["incidence_matrix_dayOnly"]
heatmap(log10.(movements.+1))

sum(movements,dims = 2)
"""
Get data from 2009 and 2019 --- rescale and get movement data from mobile phone movements
"""
T_2009 = readtable("data/2009_population_estimates.csv")
T_2019 = readtable("data/2019_population_estimates.csv")
T = join(T_2009,T_2019,on = :County)
#Rescale the size of the rural and urban population to match 2019 totals
for i = 1:47
    if !ismissing(T.Rural[i])
        T.Rural[i] = round(Int64,T.Rural[i]*T[:total_size][i]/T[:Total][i])
    else
        T.Rural[i] = 0
    end
    if !ismissing(T.Urban[i])
        T.Urban[i] = round(Int64,T.Urban[i]*T[:total_size][i]/T[:Total][i])
    else
        T.Urban[i] = 0
    end
end
CSV.write("data/combined_population_estimates.csv",T)

#Get movement matrix data
# matread("data/movements.mat")
# mv_data = matopen("data/movements.mat")
# mv_matrix = read(mv_data)["movements"]
# close(mv_data)
# CSV.write("data/combined_population_estimates.csv",T)

@load "data/mv_matrix.jld2" mv_matrix;

"""
Idea:
1) Fix the per-capita flux Nairobi <-> Mombassa to 2009 estimates.
2) Construct a gravity model prediction for other moves, varying distance exponent
3) Fit to wider area data using relative entropy
"""
urban_pop_size = [u for u in T.Urban]
x_locations = [x for x in T.x_location]
y_locations = [y for y in T.y_location]
dist_matrix = zeros(47,47) #County level distances
for i = 1:47,j=1:47
    dist_matrix[i,j] = sqrt((x_locations[i] - x_locations[j])^2 + (y_locations[i] - y_locations[j])^2)
end
ind_mombasa = findfirst(T.County .== "Mombasa")
ind_nairobi = findfirst(T.County .== "Nairobi")
normed_movements = copy(mv_matrix) #Data at the wider spatial scale
for j = 1:20
    normed_movements[:,j] = LinearAlgebra.normalize(normed_movements[:,j],1)
end
NtoMmoves = normed_movements[12,4] #Nairobi to Mombasa
MtoNmoves = normed_movements[4,12] #Mombassa to Nairobi

function prediction_moves(dist_matrix,N,α)
    P = similar(dist_matrix)
    d1,d2 = size(dist_matrix)
    for j = 1:d2
        for i = 1:d1
            if j != i
                P[i,j] = N[i]/(dist_matrix[i,j]^α)
            else
                P[i,j] = 0.
            end
        end
        P[:,j] = LinearAlgebra.normalize(P[:,j],1)
        if j == ind_nairobi
            P[ind_mombasa,j] = NtoMmoves
            x = sum(P[:,j]) - NtoMmoves
            for i = 1:d1
                if i != ind_mombasa
                    P[i,j] = P[i,j]*(1-NtoMmoves)/x
                end
            end
        end
        if j == ind_mombasa
            P[ind_nairobi,j] = MtoNmoves
            x = sum(P[:,j]) - MtoNmoves
            for i = 1:d1
                if i != ind_nairobi
                    P[i,j] = P[i,j]*(1-MtoNmoves)/x
                end
            end
        end
    end
    return P
end



function aggregate_movements(P,T)
    n = length(unique(T.wider_area))
    P̂ = zeros(n,n)
    urban_population_size_wa = [sum(T.Urban[findall(T.wider_area .== i)]) for i = 1:20]
    for i = 1:n,j=1:n
        in_i = findall(T.wider_area .== i)
        in_j = findall(T.wider_area .== j)
        for l in in_i,k in in_j
            P̂[i,j] += (T.Urban[k]/urban_population_size_wa[j])*P[l,k]
        end
    end
    P̂[4,12] = MtoNmoves #Set to correct amount
    P̂[12,4] = NtoMmoves
    return P̂
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

function gravity_model_rel_entropy(mv_data,T,dist_matrix,N,α)
    d1,d2 = size(mv_data) #<--- base distributions
    P = prediction_moves(dist_matrix,N,α) #<--- target distributions
    P̂ = aggregate_movements(P,T)
    sum_rel_entropy = 0.
    for j = 1:d1
        sum_rel_entropy += rel_entropy(P̂[:,j],mv_data[:,j])
    end
    return sum_rel_entropy
end

"""
Optimise the distance exponent
"""
function opt_func(θ)
    α = θ[1]
    error = gravity_model_rel_entropy(normed_movements,T,dist_matrix,T.Urban,α)
end
fit = optimize(opt_func,[0.6],Newton())
α_opt = Optim.minimizer(fit)
rel_min = Optim.minimum(fit)
P_opt = prediction_moves(dist_matrix,urban_pop_size,α_opt[1])

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
    title = "Fit of gravity model to aggregated movement data")
scatter!(plt,[α_opt],[rel_min],ms=5, lab = "minimum at alpha = $(round(α_opt[1],digits = 2))")

savefig(plt,"data/rel_entropy.png")
"""
Estimate time spent away by county
"""
median_length_of_trip = 5.
ρ_widerarea = [sum(mv_matrix[:,j])*median_length_of_trip/(1000*365) for j = 1:20]
ρ_county = ones(47)
for i = 1:47
    wider_area = T.wider_area[i]
    ρ_county[i] = ρ_widerarea[wider_area]
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
names = T.County
D = Dict("movement_mat" => T_opt,"County_names"=>names)
matwrite("data/mixing_matrix.mat",D)

@save "data/optimal_transition_matrix.jld2" T_opt ρ_county
