#Conversion script for 2019 data, strings -> integers for total population
#And standardise names
using DataFrames,CSV
T_2009 = readtable("data/2009_National_Estimates_of_Rural_and_Urban_Populations_by_County.csv")
T_2009 = T_2009[1:47,:]
T_2019 = readtable("data/2019-population_census-report-per-county.csv")
total_size = [parse(Int64,join(split(T_2019.Total_Population19[i],","))) for i = 1:47]
T_2019[:total_size] = total_size
T_2019[:County][5] = "Elgeyo Marakwet"
T_2019[:County][41] = "Tharaka Nithi"
"""
Convert the lat/long positions into kilometres
Equator conversion 1 deg lat = 110.57km, 1 deg long = 111.32 km
"""
T_2009[:x_location] = zeros(47)
T_2009[:y_location] = zeros(47)

for i =1:47
    x = join(split(T_2009.Location_1[i],"("))
    x = join(split(x,")"))
    y = split(x,",")
    T_2009[:x_location][i] = parse(Float64,y[2])*111.32
    T_2009[:y_location][i] = parse(Float64,y[1])*110.57
end
CSV.write("data/2009_population_estimates.csv",T_2009)
CSV.write("data/2019_population_estimates.csv",T_2019)
