#Analysis of population distribution by sub-County
push!(LOAD_PATH, "/Users/Sam/GitHub/KenyaCoV/src")
using DataFrames,CSV,MAT,Statistics,LinearAlgebra,Optim,Plots,JLD2,RData,Distances
gr()
"""
Load names of districts
"""
names = readtable("data/admin_district_names.csv")
kenya_subcounty_names = lowercase.(names.ADM2_EN[names.ADM0_EN .== "Kenya"])
kenya_county_names = sort(unique(lowercase.(names.ADM1_EN[names.ADM0_EN .== "Kenya"])))

popsize_tbl = readtable("data/distribution-of-population-by-age-sex-county-and-sub-county-kenya-2019-census-volume-iii.csv")
"""
To avoid comparison issues filter out whitespaces in age group and change everything  them out e.g. "0 - 4" != "0-4" and "WEST POKOT" != "west pokot"
to lower case in the names
"""
age_categories = [ filter(x -> !isspace(x), a) for a in popsize_tbl.Age]
sub_county_name = [lowercase(name) for name in popsize_tbl.sub_county]
county_name = unique([lowercase(name) for name in popsize_tbl.county])

check_in_data_sc = [name in sub_county_name for name in kenya_subcounty_names]
check_in_data_county = [name in county_name for name in kenya_county_names]
findall(.!check_in_data_county  )
kenya_county_names[39]
county_name[39]


findall(popsize_tbl.Age .== "0 - 4")
