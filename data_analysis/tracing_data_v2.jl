
using CSV, DataFrames, GraphRecipes, Plots, ColorBrewer, Dates

function get_shp_wards()
        shp_wards=DataFrame(CSV.File("data_analysis/shp_wards.csv"))
        shp_wards.county .= uppercase.(shp_wards.county)
        shp_wards.subcounty .= uppercase.(shp_wards.subcounty)
        shp_wards.ward .= uppercase.(shp_wards.ward)
        shp_wards_non__NA=copy(shp_wards)
        for c in unique(shp_wards.county)
                for sc in unique(shp_wards.subcounty[shp_wards.county .== c])
                        push!(shp_wards,[c sc "__NA"])
                end
                push!(shp_wards,[c "__NA" "__NA"])
        end
        #push!(shp_subcounties,["_UNKNOWN" "_UNKNOWN" 0 0])
        sort!(shp_wards,[:county,:subcounty,:ward],rev=(false,false,false))
        return shp_wards,shp_wards_non__NA
end
        shp_wards,shp_wards_non__NA=get_shp_wards()=#

colors=permutedims(unique([ColorBrewer.palette("Paired",12)[2:end-2];ColorBrewer.palette("Accent",8)[1:3];ColorBrewer.palette("Accent",8)[5:end];ColorBrewer.palette("RdBu",11)[1:4];ColorBrewer.palette("RdBu",11)[8:11]]))

# For all counties
function get_county_colors()
    tracing_data= DataFrame(CSV.File("data/cleaned_linelist20210307.csv"))
    tracing_data = tracing_data[:,[:unique_id,:contact_of_new,:age_years_,:county_of_residence,:sub_county ,:date_of_lab_confirmation]]
    tracing_data = tracing_data[tracing_data.contact_of_new .!="NA",:]
    # add columns for parsed unique_id and contact_of_new into Int64
    tracing_data[:,:unique_id_int] = parse.(Int64,SubString.(tracing_data.unique_id,6))
    tracing_data[:,:contact_of_new_int] = parse.(Int64,SubString.(tracing_data.contact_of_new,6))
    # only keep logical'ish traced contacts, where the length of contact_of_new is <= unique_id
    tracing_data=tracing_data[length.(tracing_data.contact_of_new).<=length.(tracing_data.unique_id),:]
    # add column for the county of contact_of_new
    contact_of_new_counties=[linelist_data.county_of_residence[linelist_data.unique_id.==contact_of_new]     for contact_of_new ∈ tracing_data.contact_of_new  ]
    tracing_data[:,:contact_of_new_county] = [ length(contact_of_new_counties[i])>0 ? contact_of_new_counties[i][1] : tracing_data.county_of_residence[i]      for i=1:length(contact_of_new_counties)]

    # only keep infected and infectors from county c
        #tracing_data=tracing_data[in.(tracing_data.county_of_residence,[counties_to_consider]) .| in.(tracing_data.contact_of_new_county,[counties_to_consider]),:]
    # Make DataFrame of all cases involved in contact tracing (infectors and infecteds)
    id=sort(unique(vec(Matrix(tracing_data[:,[:unique_id_int,:contact_of_new_int]]))))
    cases = DataFrame(id=id, county=["" for i=1:length(id)])
    # Get the county of each involved case
    for i=1:length(cases.id)
        a = findfirst(x->x=="Case-0"*string(cases.id[i]),linelist_data.unique_id)
        if isnothing(a)
            a = findfirst(x->x=="Case-0"*string(cases.id[i]),linelist_data.contact_of_new)
        end
        if !isnothing(a)
            cases.county[i] = linelist_data.county_of_residence[a]
        end
    end
    # Prepare color per county DataFrame
    ## >>>>>>>>>>>>>>>> sorting colors based on counties' urban %
    county_colors = (county=unique(cases.county), color=colors[1:length(unique(cases.county))])

    return county_colors
end

## Plotting tracing_data for a set of counties
function get_tracing_data(counties_to_consider)
    tracing_data= DataFrame(CSV.File("data/cleaned_linelist20210307.csv"))
    tracing_data = tracing_data[:,[:unique_id,:contact_of_new,:age_years_,:county_of_residence,:sub_county ,:date_of_lab_confirmation]]
    tracing_data = tracing_data[tracing_data.contact_of_new .!="NA",:]
    # add columns for parsed unique_id and contact_of_new into Int64
    tracing_data[:,:unique_id_int] = parse.(Int64,SubString.(tracing_data.unique_id,6))
    tracing_data[:,:contact_of_new_int] = parse.(Int64,SubString.(tracing_data.contact_of_new,6))
    # only keep logical'ish traced contacts, where the length of contact_of_new is <= unique_id
    tracing_data=tracing_data[length.(tracing_data.contact_of_new).<=length.(tracing_data.unique_id),:]
    # add column for the county of contact_of_new
    contact_of_new_counties=[linelist_data.county_of_residence[linelist_data.unique_id.==contact_of_new]     for contact_of_new ∈ tracing_data.contact_of_new  ]
    tracing_data[:,:contact_of_new_county] = [ length(contact_of_new_counties[i])>0 ? contact_of_new_counties[i][1] : tracing_data.county_of_residence[i]      for i=1:length(contact_of_new_counties)]
    # only keep infected and infectors from county c
    tracing_data=tracing_data[in.(tracing_data.county_of_residence,[counties_to_consider]) .| in.(tracing_data.contact_of_new_county,[counties_to_consider]),:]

    # Make DataFrame of all cases involved in contact tracing (infectors and infecteds)
    id=sort(unique(vec(Matrix(tracing_data[:,[:unique_id_int,:contact_of_new_int]]))))
    cases = DataFrame(id=id, county=["" for i=1:length(id)], color=[RGB(0,0,0) for i=1:length(id)])
    # Get the county of each involved case
    for i=1:length(cases.id)
        a = findfirst(x->x=="Case-0"*string(cases.id[i]),linelist_data.unique_id)
        if isnothing(a)
            a = findfirst(x->x=="Case-0"*string(cases.id[i]),linelist_data.contact_of_new)
        end
        if !isnothing(a)
            cases.county[i] = linelist_data.county_of_residence[a]
        end
    end
    # Prepare color per county DataFrame
    ## >>>>>>>>>>>>>>>> sorting colors based on counties' urban %
    county_colors = get_county_colors()
    # Determine what color is each case depending on its coounty
    for i=1:length(cases.id)
        cases.color[i]=county_colors.color[county_colors.county.==cases.county[i]][1]
    end
    # Giving a unique id to each case
    tracing_data[:,:unique_id_node]=[findfirst(x->x==id,cases.id)  for id ∈ tracing_data.unique_id_int]
        tracing_data[:,:contact_of_new_node]=[findfirst(x->x==id,cases.id)  for id ∈ tracing_data.contact_of_new_int]
    # Add number of infecteds caused by each case
    cases[:,:nb_traced_inecteds]=zeros(Int64, length(cases.id))
    for i=1:length(cases.id)
        cases.nb_traced_inecteds[i]=length(findall(x->x==cases.id[i],tracing_data.contact_of_new_int))
    end
    # Add date of lab confirmation
    tracing_data.date_of_lab_confirmation.=replace.(tracing_data.date_of_lab_confirmation,"-"=>"/")
    cases[:,:date_of_lab_confirmation]=[Date("20/02/2020", dateformat"dd/mm/yyyy") for i=1:length(cases.id)]
    for i=1:length(cases.id)
        d = tracing_data.date_of_lab_confirmation[tracing_data.unique_id_int.==cases.id[i]]
        if length(d)>0
            cases.date_of_lab_confirmation[i]=Date(tracing_data.date_of_lab_confirmation[tracing_data.unique_id_int.==cases.id[i]][1], dateformat"dd/mm/yyyy")#length(findall(x->x==cases.id[i],tracing_data.contact_of_new_int))
        end
    end
    # Preparing transition matrix
    transition_matrix = zeros(Int64,length(cases.id),length(cases.id))
        for i=1:length(tracing_data.contact_of_new)
            transition_matrix[tracing_data.contact_of_new_node[i],tracing_data.unique_id_node[i]] = 1
        end
    return cases, transition_matrix
end
    cases, transition_matrix = get_tracing_data(["KAJIADO"])
Date(tracing_data.date_of_lab_confirmation[1])
function draw_graph(counties_to_consider)
    cases, transition_matrix = get_tracing_data(counties_to_consider)
        plt= graphplot(transition_matrix, method= :tree, curvature_scalar=0, nodeshape = :circle,
                    nodecolor=cases.color,
                    node_weights=cases.nb_traced_inecteds.+1,
                    title="Contact tracing for counties "*string(counties_to_consider)*"\n includes individuals from "*string(unique(cases.county)), titlefontsize=8,
                    #label=permutedims(SubString.(cases.county,1,2)), legend=:topleft
                    #names=SubString.(cases.county,1,2), fontsize=4 #string.(cases.id)
                    #names=string.(cases.nb_traced_inecteds), fontsize=7
                    )
        return cases, transition_matrix, plt
end
    cases, transition_matrix, plt = draw_graph(["KAJIADO"])
    display(plt)

for c ∈ counties.county
    cases, transition_matrix, plt = draw_graph([c])
    display(plt)
end

# color coding
scatter([(i,5) for i=1:17], ylims=(0,10),
        linewidth=false,legend=false,
        markersize=20,
        xticks=(1:1:17,county_colors.county), tickfontsize=3,
        yticks=false,
        color=county_colors.color,
        size=(400,100))


using ColorSchemes
get(colorschemes[:rainbow], collect(1/47:1/47:1))










#
