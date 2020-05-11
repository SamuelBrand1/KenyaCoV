push!(LOAD_PATH, joinpath(homedir(),"GitHub/KenyaCoV/src"))
using Parameters,Distributions,DifferentialEquations,JLD2,DataFrames,FileIO,DelimitedFiles,RecursiveArrayTools,Plots
using CSV,ExcelFiles,DataFrames
using Statistics: median, quantile
using LinearAlgebra: eigen


#Load data
# @load joinpath(homedir(),"Github/KenyaCoVOutputs/sims_consensus_baseline_vs2.jld2") sims_baseline
@load joinpath(homedir(),"Github/KenyaCoVOutputs/sims_consensus_end_lockdown.jld2") sims_end_regional_lockdown
# @load joinpath(homedir(),"Github/KenyaCoVOutputs/sims_consensus_open_schools_june.jld2") sims_open_schools_june
# @load joinpath(homedir(),"Github/KenyaCoVOutputs/sims_consensus_open_schools_august.jld2") sims_open_schools_august


## Functions --- to be used



include("hospitalisations.jl");
counties = CSV.read(joinpath(homedir(),"Github/KenyaCoV/data/2019_census_age_pyramids_counties.csv"))
names = counties.county
sims = sims_end_regional_lockdown

### Peak calculations by county



@time a,peak_A,peak_A_value = first_introduction_time_peak_peak_value(sims,1)
b,peak_M,peak_M_value = first_introduction_time_peak_peak_value(sims,2)
c,peak_V,peak_V_value = first_introduction_time_peak_peak_value(sims,3)

first_times = get_first_time(a,b,c)


df = print_onset_report(first_times,peak_A,peak_A_value,peak_M,peak_M_value,peak_V,peak_V_value,names,"KenyaCoVOutputs/peak_report_end_regional_lockdown.csv")

spare_capacity_H_by_county[2:end]
spare_capacity_ICU_by_county[30]

function print_hosp_report(sims)
    df = DataFrame(county = names)

    hosp_occup_per_sim,ICU_occup_per_sim,new_ICU_per_sim,death_incidence_per_sim = total_hospital_outcomes_per_sim(sims,30)
cum_death_per_sim = similar(death_incidence_per_sim)
cum_hosp_per_sim = similar(death_incidence_per_sim)
cum_ICU_per_sim = similar(death_incidence_per_sim)

for k = 1:1000
    cum_death_per_sim[k,:] .= cumsum(death_incidence_per_sim[k,:])
    cum_death_per_sim[k,:] .= cumsum(death_incidence_per_sim[k,:])

end

df = DataFrame(county = names)
hosp_maximums = [(hosp_max,hosp_max_time) = findmax(hosp_occup_per_sim[k,:]) for k = 1:1000]
ICU_maximums = [(ICU_max,ICU_max_time) = findmax(ICU_occup_per_sim[k,:]) for k = 1:1000]
perc_hosp_exceeds = [h[1]/spare_capacity_H_by_county[31] for h in hosp_maximums]
perc_ICU_exceeds = [h[1]/spare_capacity_ICU_by_county[31] for h in ICU_maximums]

chance_exceed_hosp = mean([p >= 1 for p in perc_hosp_exceeds])
day_hosp_exceeded = filter(!isnothing,[findfirst(hosp_occup_per_sim[k,:].>spare_capacity_H_by_county[31]) for k = 1:1000])
median_day_ICU_exceeded = median(filter(!isnothing,[findfirst(ICU_occup_per_sim[k,:].>spare_capacity_ICU_by_county[31]) for k = 1:1000]))

chance_exceed_ICU = mean([p >= 1 for p in perc_ICU_exceeds])







## Cumulative incidence
cum_incidence_A,cum_incidence_M,cum_incidence_V,cum_incidence_H = cum_incidence_for_each_sim_by_type(sims,30)
cum_incidence_total = cum_incidence_A.+cum_incidence_M.+cum_incidence_V
incidence_A = diff(cum_incidence_A,dims = 2)
incidence_M = diff(cum_incidence_M,dims = 2)
incidence_V = diff(cum_incidence_V,dims = 2)
incidence_H = diff(cum_incidence_H,dims = 2)
total_incidence = incidence_A.+incidence_M.+incidence_V




#Cumulative values for spreadsheet days 30,45,60,90,180,365(end)
#Total infections
median_total_infections = [ median(cum_incidence_total[:,t+1]) for t in [30,45,60,90,180,360,365] ]
lb_total_infections = [ quantile(cum_incidence_total[:,t+1],0.025) for t in [30,45,60,90,180,360,365] ]
ub_total_infections = [ quantile(cum_incidence_total[:,t+1],0.975) for t in [30,45,60,90,180,360,365] ]
total_inf_array = hcat(median_total_infections,lb_total_infections,ub_total_infections)
writedlm("output/save_total_inf_"*nametag*".csv",total_inf_array,',')
median_total_infections_at_peak = median([ cum_incidence_total[k,t+1] for (k,t) in enumerate(peak_days) ])
lb_total_infections_at_peak = quantile([ cum_incidence_total[k,t+1] for (k,t) in enumerate(peak_days) ],0.025)
ub_total_infections_at_peak = quantile([ cum_incidence_total[k,t+1] for (k,t) in enumerate(peak_days) ],0.975)
writedlm("output/save_total_inf_atpeak"*nametag*".csv",[median_total_infections_at_peak,lb_total_infections_at_peak,ub_total_infections_at_peak],',')

#Mild infection
median_mild_infections = [ median(cum_incidence_M[:,t+1]) for t in [30,45,60,90,180,360,365] ]
lb_mild_infections = [ quantile(cum_incidence_M[:,t+1],0.025) for t in [30,45,60,90,180,360,365] ]
ub_mild_infections = [ quantile(cum_incidence_M[:,t+1],0.975) for t in [30,45,60,90,180,360,365] ]
total_mild_array = hcat(median_mild_infections,lb_mild_infections,ub_mild_infections)
writedlm("output/save_total_mild_inf_"*nametag*".csv",total_mild_array,',')

median_mild_infections_at_peak = median([ cum_incidence_M[k,t+1] for (k,t) in enumerate(peak_days) ])
lb_mild_infections_at_peak = quantile([ cum_incidence_M[k,t+1] for (k,t) in enumerate(peak_days) ],0.025)
ub_mild_infections_at_peak = quantile([ cum_incidence_M[k,t+1] for (k,t) in enumerate(peak_days) ],0.975)
writedlm("output/save_total_mild_inf_atpeak"*nametag*".csv",[median_mild_infections_at_peak,lb_mild_infections_at_peak,ub_mild_infections_at_peak],',')

#Incidence
incidence_A = diff(cum_incidence_A,dims = 2)
incidence_M = diff(cum_incidence_M,dims = 2)
incidence_V = diff(cum_incidence_V,dims = 2)
incidence_H = diff(cum_incidence_H,dims = 2)
total_incidence = incidence_A.+incidence_M.+incidence_V

median_total_incidence = [ median(total_incidence[:,t]) for t in [30,45,60,90,180,360,365] ]
lb_total_incidence = [ quantile(total_incidence[:,t],0.025) for t in [30,45,60,90,180,360,365] ]
ub_total_incidence = [ quantile(total_incidence[:,t],0.975) for t in [30,45,60,90,180,360,365] ]
incidence_array = hcat(median_total_incidence,lb_total_incidence,ub_total_incidence)
writedlm("output/save_incidence_inf_"*nametag*".csv",incidence_array,',')
median_total_incidence_at_peak = median([ total_incidence[k,t] for (k,t) in enumerate(peak_days) ])
lb_total_incidence_at_peak = quantile([ total_incidence[k,t] for (k,t) in enumerate(peak_days) ],0.025)
ub_total_incidence_at_peak = quantile([ total_incidence[k,t] for (k,t) in enumerate(peak_days) ],0.975)
writedlm("output/save_incidence_inf_atpeak"*nametag*".csv",[median_total_incidence_at_peak,lb_total_incidence_at_peak,ub_total_incidence_at_peak],',')


#Mild infection
median_mild_incidence = [ median(incidence_M[:,t]) for t in [30,45,60,90,180,360,365] ]
lb_mild_incidence = [ quantile(incidence_M[:,t],0.025) for t in [30,45,60,90,180,360,365] ]
ub_mild_incidence = [ quantile(incidence_M[:,t],0.975) for t in [30,45,60,90,180,360,365] ]
mild_incidence_array = hcat(median_mild_incidence,lb_mild_incidence,ub_mild_incidence)
writedlm("output/save_mild_incidence_inf_"*nametag*".csv",mild_incidence_array,',')
median_mild_incidence_at_peak = median([ incidence_M[k,t] for (k,t) in enumerate(peak_days) ])
lb_mild_incidence_at_peak = quantile([ incidence_M[k,t] for (k,t) in enumerate(peak_days) ],0.025)
ub_mild_incidence_at_peak = quantile([ incidence_M[k,t] for (k,t) in enumerate(peak_days) ],0.975)
writedlm("output/save_mild_incidence_inf_atpeak"*nametag*".csv",[median_mild_incidence_at_peak,lb_mild_incidence_at_peak,ub_mild_incidence_at_peak],',')



#Severe infection
median_severe_incidence = [ median(incidence_H[:,t]) for t in [30,45,60,90,180,360,365] ]
lb_severe_incidence = [ quantile(incidence_H[:,t],0.025) for t in [30,45,60,90,180,360,365] ]
ub_severe_incidence = [ quantile(incidence_H[:,t],0.975) for t in [30,45,60,90,180,360,365] ]
median_severe_incidence_at_peak = median([ incidence_H[k,t] for (k,t) in enumerate(peak_days) ])
lb_severe_incidence_at_peak = quantile([ incidence_H[k,t] for (k,t) in enumerate(peak_days) ],0.025)
ub_severe_incidence_at_peak = quantile([ incidence_H[k,t] for (k,t) in enumerate(peak_days) ],0.975)


## Hospitalisation outcomes
hosp_occup_per_sim,ICU_occup_per_sim,new_ICU_per_sim,death_incidence_per_sim = total_hospital_outcomes_per_sim(sims,30)
cum_death_per_sim = similar(death_incidence_per_sim)
for k = 1:1000
    cum_death_per_sim[k,:] .= cumsum(death_incidence_per_sim[k,:])
end
#Peak quantities

death_incidence_peaks = find_peaks(death_incidence_per_sim)
D_ts_at_peak = [median(death_incidence_peaks),quantile(death_incidence_peaks,0.025),quantile(death_incidence_peaks,0.975)]
writedlm("output/save_death_incidence_atpeak"*nametag*".csv",D_ts_at_peak,',')
ICU_incidence_peaks = find_peaks(new_ICU_per_sim)
ICU_incidence_at_peak = [median(ICU_incidence_peaks),quantile(ICU_incidence_peaks,0.025),quantile(ICU_incidence_peaks,0.975)]
writedlm("output/save_ICU_incidence_atpeak"*nametag*".csv",ICU_incidence_at_peak,',')

H_incidence_peaks = find_peaks(incidence_H)
H_incidence_at_peak = [median(H_incidence_peaks),quantile(H_incidence_peaks,0.025),quantile(H_incidence_peaks,0.975)]
writedlm("output/save_hosp_incidence_atpeak"*nametag*".csv",H_incidence_at_peak,',')

H_occup_peaks = find_peaks(hosp_occup_per_sim)
H_occup_at_peak = [median(H_occup_peaks),quantile(H_occup_peaks,0.025),quantile(H_occup_peaks,0.975)]
writedlm("output/save_hosp_occup_atpeak"*nametag*".csv",H_occup_at_peak,',')

ICU_occup_peaks = find_peaks(ICU_occup_per_sim)
ICU_occup_at_peak = [median(ICU_occup_peaks),quantile(ICU_occup_peaks,0.025),quantile(ICU_occup_peaks,0.975)]
writedlm("output/save_ICU_occup_atpeak"*nametag*".csv",ICU_occup_at_peak,',')

median_hosp_occup = [ median(hosp_occup_per_sim[:,t]) for t in [30,45,60,90,180,360,365] ]
lb_hosp_occup  = [ quantile(hosp_occup_per_sim[:,t],0.025) for t in [30,45,60,90,180,360,365] ]
ub_hosp_occup  = [ quantile(hosp_occup_per_sim[:,t],0.975) for t in [30,45,60,90,180,360,365] ]
hosp_occup_array = hcat(median_hosp_occup,lb_hosp_occup,ub_hosp_occup)
writedlm("output/save_hosp_occup_"*nametag*".csv",hosp_occup_array,',')

median_H = [ median(incidence_H[:,t]) for t in [30,45,60,90,180,360,365]]
lb_H = [ quantile(incidence_H[:,t],0.025) for t in [30,45,60,90,180,360,365] ]
ub_H = [ quantile(incidence_H[:,t],0.975) for t in [30,45,60,90,180,360,365] ]
hosp_incidence_array = hcat(median_H,lb_H,ub_H)
writedlm("output/save_hosp_incidence_"*nametag*".csv",hosp_incidence_array',',')


median_ICU_incidence = [ median(new_ICU_per_sim[:,t]) for t in [30,45,60,90,180,360,365] ]
lb_ICU_incidence  = [ quantile(new_ICU_per_sim[:,t],0.025) for t in [30,45,60,90,180,360,365] ]
ub_ICU_incidence  = [ quantile(new_ICU_per_sim[:,t],0.975) for t in [30,45,60,90,180,360,365] ]
ICU_incidence_array = vcat(median_ICU_incidence',lb_ICU_incidence',ub_ICU_incidence')
writedlm("output/save_ICU_incidence_"*nametag*".csv",ICU_incidence_array,',')
median_ICU_incidence_ts = [ median(new_ICU_per_sim[:,t]) for t in 1:365 ]



median_ICU_occup = [ median(ICU_occup_per_sim[:,t]) for t in [30,45,60,90,180,360,365] ]
lb_ICU_occup  = [ quantile(ICU_occup_per_sim[:,t],0.025) for t in [30,45,60,90,180,360,365] ]
ub_ICU_occup  = [ quantile(ICU_occup_per_sim[:,t],0.975) for t in [30,45,60,90,180,360,365] ]
ICU_occup_array = vcat(median_ICU_occup',lb_ICU_occup',ub_ICU_occup')
writedlm("output/save_ICU_occup_"*nametag*".csv",ICU_occup_array,',')

median_death_incidence = [ median(death_incidence_per_sim[:,t]) for t in [30,45,60,90,180,360,365] ]
lb_death_incidence  = [ quantile(death_incidence_per_sim[:,t],0.025) for t in [30,45,60,90,180,360,365] ]
ub_death_incidence  = [ quantile(death_incidence_per_sim[:,t],0.975) for t in [30,45,60,90,180,360,365] ]
death_incidence_array = vcat(median_death_incidence',lb_death_incidence',ub_death_incidence')
writedlm("output/save_death_incidence_"*nametag*".csv",death_incidence_array,',')


median_death_cum = [ median(cum_death_per_sim[:,t]) for t in [30,45,60,90,180,360,365] ]
lb_cum_death_incidence  = [ quantile(cum_death_per_sim[:,t],0.025) for t in [30,45,60,90,180,360,365] ]
ub_cum_death_incidence  = [ quantile(cum_death_per_sim[:,t],0.975) for t in [30,45,60,90,180,360,365] ]

cum_death_incidence_array = vcat(median_death_cum',lb_cum_death_incidence',ub_cum_death_incidence')
writedlm("output/save_death_cum_"*nametag*".csv",cum_death_incidence_array,',')


## Plotting
using Dates
gr()
first_of_months = [dayofyear(2020,m,1) for m = 4:12] .- dayofyear(2020,3,13)
median_deaths_ts = [ median(death_incidence_per_sim[:,t]) for t in 1:365 ]
lb_deaths_ts = [ quantile(death_incidence_per_sim[:,t],0.025) for t in 1:365 ]
ub_deaths_ts = [ quantile(death_incidence_per_sim[:,t],0.975) for t in 1:365 ]



median_H_ts = [ median(incidence_H[:,t]) for t in 1:365 ]
lb_H_ts = [ quantile(incidence_H[:,t],0.025) for t in 1:365 ]
ub_H_ts = [ quantile(incidence_H[:,t],0.975) for t in 1:365 ]





median_Hocc_ts = [ median(hosp_occup_per_sim[:,t]) for t in 1:365 ]
lb_HOcc_ts = [ quantile(hosp_occup_per_sim[:,t],0.025) for t in 1:365 ]
ub_HOcc_ts = [ quantile(hosp_occup_per_sim[:,t],0.975) for t in 1:365 ]




median_ICUocc_ts = [ median(ICU_occup_per_sim[:,t]) for t in 1:365 ]
lb_ICUOcc_ts = [ quantile(ICU_occup_per_sim[:,t],0.025) for t in 1:365 ]
ub_ICUOcc_ts = [ quantile(ICU_occup_per_sim[:,t],0.975) for t in 1:365 ]


plt_incidence_nairobi = plot(0:364,medi,
                        lab = "Hospitalisations",
                        lw = 3,color = :red,
                        ribbon = (median_H_ts.-lb_H_ts,ub_H_ts .- median_H_ts),
                        fillalpha = 0.25,
                        xticks = (first_of_months,["Apr","May","June","July","Aug","Sept","Oct","Nov","Dec"]),
                        xlims = (0.,275.),
                        title = "Incidence of disease")
#
#  plot!(plt_incidence,0:364,median_deaths_ts,
#                         lab = "Deaths",
#                         lw = 3,color = :black,
#                         ribbon = (median_deaths_ts.-lb_deaths_ts,ub_deaths_ts .- median_deaths_ts),
#                         fillalpha = 0.5,
#                         xticks = (first_of_months,["Apr","May","June","July","Aug","Sept","Oct","Nov","Dec"]),
#                         xlims = (0.,275.),
#                         ylabel = "Daily Incidence" )

plot!(plt_incidence_both,0:364,median_H_ts,
                        lab = "",
                        lw = 3,color = :red,ls=:dot,
                        ribbon = (median_H_ts.-lb_H_ts,ub_H_ts .- median_H_ts),
                        fillalpha = 0.15,
                        xticks = (first_of_months,["Apr","May","June","July","Aug","Sept","Oct","Nov","Dec"]),
                        xlims = (0.,275.),
                        title = "Incidence of disease")
plot!(plt_incidence_both,0:364,median_deaths_ts,
                        lab = "",
                        lw = 3,color = :black,ls=:dot,
                        ribbon = (median_deaths_ts.-lb_deaths_ts,ub_deaths_ts .- median_deaths_ts),
                        fillalpha = 0.15,
                        xticks = (first_of_months,["Apr","May","June","July","Aug","Sept","Oct","Nov","Dec"]),
                        xlims = (0.,275.),
                        ylabel = "Daily Incidence" )


# plt_health_usage = plot(0:364,median_Hocc_ts,
#                         lab = "Hospital beds",
#                         lw = 3,color = :blue,
#                         ribbon = (median_Hocc_ts.-lb_HOcc_ts,ub_HOcc_ts .- median_Hocc_ts),
#                         fillalpha = 0.25,
#                         xticks = (first_of_months,["Apr","May","June","July","Aug","Sept","Oct","Nov","Dec"]),
#                         xlims = (0.,275.),
#                         ylabel = "Daily Occupancy",
#                         title = "Health system usage")
# plot!(plt_health_usage,0:364,median_ICUocc_ts,
#                         lab = "ICU beds",
#                         lw = 3,color = :green,
#                         ribbon = (median_ICUocc_ts.-lb_ICUOcc_ts,ub_ICUOcc_ts .- median_ICUocc_ts),
#                         fillalpha = 0.25,
#                         xticks = (first_of_months,["Apr","May","June","July","Aug","Sept","Oct","Nov","Dec"]),
#                         xlims = (0.,275.),
#                         ylabel = "Daily Occupancy",
#                         title = "Health system usage")

plt_health_usage = plot(0:364,median_Hocc_ts,
                        lab = "",
                        lw = 3,color = :blue,ls=:dot,
                        ribbon = (median_Hocc_ts.-lb_HOcc_ts,ub_HOcc_ts .- median_Hocc_ts),
                        fillalpha = 0.25,
                        xticks = (first_of_months,["Apr","May","June","July","Aug","Sept","Oct","Nov","Dec"]),
                        xlims = (0.,275.),
                        ylabel = "Daily Occupancy",
                        title = "Health system usage")
plot!(plt_health_usage_both,0:364,median_ICUocc_ts,
                        lab = "",ls = :dot,
                        lw = 3,color = :green,
                        ribbon = (median_ICUocc_ts.-lb_ICUOcc_ts,ub_ICUOcc_ts .- median_ICUocc_ts),
                        fillalpha = 0.25,
                        xticks = (first_of_months,["Apr","May","June","July","Aug","Sept","Oct","Nov","Dec"]),
                        xlims = (0.,275.),
                        ylabel = "Daily Occupancy",
                        title = "Health system usage")

plot(0:364,[ median(cum_incidence_total[:,t+1]) for t in 1:365],lab = "Cum. incidence",
                        lw = 3,color = :blue,
                        fillalpha = 0.25,
                        xticks = (first_of_months,["Apr","May","June","July","Aug","Sept","Oct","Nov","Dec"]),
                        xlims = (0.,275.),
                        ylabel = "Total infections",title = "Controlled epidemic")
savefig("output/controlled_epidemic.png")

median([total_incidence[:,t] for t = 1:365])
savefig(plt_incidence,"output/plt_incidence.png")
savefig(plt_incidence_both,"output/plt_incidence_with_controls.png")

savefig(plt_health_usage,"output/plt_health_usage.png")
savefig(plt_health_usage_both,"output/plt_health_usage_with_controls.png")


# @save "output/plt_incidence.jld2" plt_incidence
# @save "output/plt_health_usage.jld2" plt_health_usage

@load "output/plt_incidence.jld2" plt_incidence
@load "output/plt_health_usage.jld2" plt_health_usage

plt_incidence_both = deepcopy(plt_incidence)
plt_health_usage_both = deepcopy(plt_health_usage)



# function cum_incidence_for_each_sim_by_type(sims)
#     n = length(sims.u)
#     T = length(sims.u[1])
#     cum_incidence_A = zeros(n,T)
#     cum_incidence_M = zeros(n,T)
#     cum_incidence_V = zeros(n,T)
#     cum_incidence_H = zeros(n,T)
#
#     for k = 1:n,t = 1:T
#         cum_incidence_A[k,t] = sum(sims.u[k][t][:,:,1])
#         cum_incidence_M[k,t] = sum(sims.u[k][t][:,:,2])
#         cum_incidence_V[k,t] = sum(sims.u[k][t][:,:,3])
#         cum_incidence_H[k,t] = sum(sims.u[k][t][:,:,4])
#     end
#     return cum_incidence_A,cum_incidence_M,cum_incidence_V,cum_incidence_H
# end


## County plots
using CSV
counties = CSV.read("data/2019_census_age_pyramids_counties.csv")
Nairobi_index = findfirst(counties.county .== "Nairobi")
Mombassa_index = findfirst(counties.county .== "Mombasa")
Kwale_index = findfirst(counties.county .== "Kwale")
Kilifi_index = findfirst(counties.county .== "Kilifi")
Mandera_index  = findfirst(counties.county .== "Mandera")

sims = sims_end_regional_lockdown



plt_incidence_nairobi,plt_usage_nairobi = give_plots_for_county(sims,30)
plot!(plt_incidence_nairobi,title = "Nairobi (SD + regional lockdown relaxed May 16th)")
plot!(plt_usage_nairobi,[0,365],[spare_capacity_H_by_county[30+1],spare_capacity_H_by_county[30+1]],
        title = "Nairobi (SD + regional lockdown relaxed May 16th)",lw = 2,ls = :dash,lab = "spare hosp. capacity",color = :blue)
plot!(plt_usage_nairobi,[0,365],[spare_capacity_ICU_by_county[30+1],spare_capacity_ICU_by_county[30+1]],
        title = "Nairobi (SD + regional lockdown relaxed May 16th)",lw = 2,ls = :dash,lab = "spare ICU capacity",color = :green)

savefig(plt_incidence_nairobi,"nairobi_incidence_regional_relaxation.png")
savefig(plt_usage_nairobi,"nairobi_health_usage_regional_relaxation.png")
plt_incidence_mombasa,plt_usage_mombasa = give_plots_for_county(sims,Mombassa_index)
plot!(plt_incidence_mombasa,title = "Mombasa (SD + regional lockdown relaxed May 16th)")
plot!(plt_usage_mombasa,title = "Mombasa (SD + regional lockdown relaxed May 16th)")
savefig(plt_incidence_mombasa,"mombasa_incidence_regional_relaxation.png")
savefig(plt_usage_mombasa,"mombasa_health_usage_regional_relaxation.png")
plt_incidence_rest,plt_usage_rest = give_plots_for_county(sims,setdiff(1:47,[28,30]))
plot!(plt_incidence_rest,title = "Rest of Kenya (SD + regional lockdown relaxed May 16th)")
plot!(plt_usage_rest,title = "Rest of Kenya (SD + regional lockdown relaxed May 16th")
savefig(plt_incidence_rest,"rest_of_country_incidence_regional_relaxation.png")
savefig(plt_usage_rest,"rest_of_country_health_usage_regional_relaxation.png")
