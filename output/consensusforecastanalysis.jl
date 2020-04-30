push!(LOAD_PATH, joinpath(homedir(),"GitHub/KenyaCoV/src"))
using Plots,Parameters,Distributions,DifferentialEquations,JLD2,DataFrames,StatsPlots,FileIO,MAT,RecursiveArrayTools
using Statistics: median, quantile
using LinearAlgebra: eigen


#Load data
@load joinpath(homedir(),"Github/KenyaCoVOutputs/sims_consensus_baseline.jld2") sims_baseline
@load joinpath(homedir(),"Github/KenyaCoVOutputs/sims_consensus_baseline_scaled.jld2") sims_baseline_scaled
@load joinpath(homedir(),"Github/KenyaCoVOutputs/sims_consensus_control.jld2") sims_controls
@load joinpath(homedir(),"Github/KenyaCoVOutputs/sims_consensus_control_scaled.jld2") sims_controls_scaled

## Functions --- to be used



include("hospitalisations.jl");


## Cumulative incidence
cum_incidence_A,cum_incidence_M,cum_incidence_V,cum_incidence_H = cum_incidence_for_each_sim_by_type(sims_controls)
cum_incidence_total = cum_incidence_A.+cum_incidence_M.+cum_incidence_V
total_incidence_peaks = find_peak_time_by_sim(sims_controls)
peak_days = [x[2] for x in total_incidence_peaks]
#Cumulative values for spreadsheet days 30,45,60,90,180,365(end)
#Total infections
median_total_infections = [ median(cum_incidence_total[:,t+1]) for t in [30,45,60,90,180,360,365] ]
lb_total_infections = [ quantile(cum_incidence_total[:,t+1],0.025) for t in [30,45,60,90,180,360,365] ]
ub_total_infections = [ quantile(cum_incidence_total[:,t+1],0.975) for t in [30,45,60,90,180,360,365] ]
median_total_infections_at_peak = median([ cum_incidence_total[k,t+1] for (k,t) in enumerate(peak_days) ])
lb_total_infections_at_peak = quantile([ cum_incidence_total[k,t+1] for (k,t) in enumerate(peak_days) ],0.025)
ub_total_infections_at_peak = quantile([ cum_incidence_total[k,t+1] for (k,t) in enumerate(peak_days) ],0.975)
#Mild infection
median_mild_infections = [ median(cum_incidence_M[:,t+1]) for t in [30,45,60,90,180,360,365] ]
lb_mild_infections = [ quantile(cum_incidence_M[:,t+1],0.025) for t in [30,45,60,90,180,360,365] ]
ub_mild_infections = [ quantile(cum_incidence_M[:,t+1],0.975) for t in [30,45,60,90,180,360,365] ]
median_mild_infections_at_peak = median([ cum_incidence_M[k,t+1] for (k,t) in enumerate(peak_days) ])
lb_mild_infections_at_peak = quantile([ cum_incidence_M[k,t+1] for (k,t) in enumerate(peak_days) ],0.025)
ub_mild_infections_at_peak = quantile([ cum_incidence_M[k,t+1] for (k,t) in enumerate(peak_days) ],0.975)

#Incidence
incidence_A = diff(cum_incidence_A,dims = 2)
incidence_M = diff(cum_incidence_M,dims = 2)
incidence_V = diff(cum_incidence_V,dims = 2)
incidence_H = diff(cum_incidence_H,dims = 2)
total_incidence = incidence_A.+incidence_M.+incidence_V

median_total_incidence = [ median(total_incidence[:,t]) for t in [30,45,60,90,180,360,365] ]
lb_total_incidence = [ quantile(total_incidence[:,t],0.025) for t in [30,45,60,90,180,360,365] ]
ub_total_incidence = [ quantile(total_incidence[:,t],0.975) for t in [30,45,60,90,180,360,365] ]
median_total_incidence_at_peak = median([ total_incidence[k,t] for (k,t) in enumerate(peak_days) ])
lb_total_incidence_at_peak = quantile([ total_incidence[k,t] for (k,t) in enumerate(peak_days) ],0.025)
ub_total_incidence_at_peak = quantile([ total_incidence[k,t] for (k,t) in enumerate(peak_days) ],0.975)
#Mild infection
median_mild_incidence = [ median(incidence_M[:,t]) for t in [30,45,60,90,180,360,365] ]
lb_mild_incidence = [ quantile(incidence_M[:,t],0.025) for t in [30,45,60,90,180,360,365] ]
ub_mild_incidence = [ quantile(incidence_M[:,t],0.975) for t in [30,45,60,90,180,360,365] ]
median_mild_incidence_at_peak = median([ incidence_M[k,t] for (k,t) in enumerate(peak_days) ])
lb_mild_incidence_at_peak = quantile([ incidence_M[k,t] for (k,t) in enumerate(peak_days) ],0.025)
ub_mild_incidence_at_peak = quantile([ incidence_M[k,t] for (k,t) in enumerate(peak_days) ],0.975)
#Severe infection
median_severe_incidence = [ median(incidence_H[:,t]) for t in [30,45,60,90,180,360,365] ]
lb_severe_incidence = [ quantile(incidence_H[:,t],0.025) for t in [30,45,60,90,180,360,365] ]
ub_severe_incidence = [ quantile(incidence_H[:,t],0.975) for t in [30,45,60,90,180,360,365] ]
median_severe_incidence_at_peak = median([ incidence_H[k,t] for (k,t) in enumerate(peak_days) ])
lb_severe_incidence_at_peak = quantile([ incidence_H[k,t] for (k,t) in enumerate(peak_days) ],0.025)
ub_severe_incidence_at_peak = quantile([ incidence_H[k,t] for (k,t) in enumerate(peak_days) ],0.975)


## Hospitalisation outcomes
hosp_occup_per_sim,ICU_occup_per_sim,new_ICU_per_sim,death_incidence_per_sim = total_hospital_outcomes_per_sim(sims_baseline)
cum_death_per_sim = similar(death_incidence_per_sim)
for k = 1:1000
    cum_death_per_sim[k,:] .= cumsum(death_incidence_per_sim[k,:])
end

median_hosp_occup = [ median(hosp_occup_per_sim[:,t]) for t in [30,45,60,90,180,360,365] ]
lb_hosp_occup  = [ quantile(hosp_occup_per_sim[:,t],0.025) for t in [30,45,60,90,180,360,365] ]
ub_hosp_occup  = [ quantile(hosp_occup_per_sim[:,t],0.975) for t in [30,45,60,90,180,360,365] ]


median_ICU_incidence = [ median(new_ICU_per_sim[:,t]) for t in [30,45,60,90,180,360,365] ]
lb_ICU_incidence  = [ quantile(new_ICU_per_sim[:,t],0.025) for t in [30,45,60,90,180,360,365] ]
ub_ICU_incidence  = [ quantile(new_ICU_per_sim[:,t],0.975) for t in [30,45,60,90,180,360,365] ]
median_ICU_incidence_ts = [ median(new_ICU_per_sim[:,t]) for t in 1:365 ]
(ICU_incidence_at,day_) = findmax(median_ICU_incidence_ts)
ICU_incidence_at_peak = [median_ICU_incidence_ts[day_],quantile(new_ICU_per_sim[:,day_],0.025),quantile(new_ICU_per_sim[:,day_],0.975)]


median_ICU_occup = [ median(ICU_occup_per_sim[:,t]) for t in [30,45,60,90,180,360,365] ]
lb_ICU_occup  = [ quantile(ICU_occup_per_sim[:,t],0.025) for t in [30,45,60,90,180,360,365] ]
ub_ICU_occup  = [ quantile(ICU_occup_per_sim[:,t],0.975) for t in [30,45,60,90,180,360,365] ]


median_death_incidence = [ median(death_incidence_per_sim[:,t]) for t in [30,45,60,90,180,360,365] ]
lb_death_incidence  = [ quantile(death_incidence_per_sim[:,t],0.025) for t in [30,45,60,90,180,360,365] ]
ub_death_incidence  = [ quantile(death_incidence_per_sim[:,t],0.975) for t in [30,45,60,90,180,360,365] ]

median_death_cum = [ median(cum_death_per_sim[:,t]) for t in [30,45,60,90,180,360,365] ]
lb_cum_death_incidence  = [ quantile(cum_death_per_sim[:,t],0.025) for t in [30,45,60,90,180,360,365] ]
ub_cum_death_incidence  = [ quantile(cum_death_per_sim[:,t],0.975) for t in [30,45,60,90,180,360,365] ]


## Plotting
using Dates
gr()
first_of_months = [dayofyear(2020,m,1) for m = 4:12] .- dayofyear(2020,3,13)
median_deaths_ts = [ median(death_incidence_per_sim[:,t]) for t in 1:365 ]
lb_deaths_ts = [ quantile(death_incidence_per_sim[:,t],0.025) for t in 1:365 ]
ub_deaths_ts = [ quantile(death_incidence_per_sim[:,t],0.975) for t in 1:365 ]
(max_D_inc,day_) = findmax(median_deaths_ts)
D_ts_at_peak = [median_deaths_ts[day_],lb_deaths_ts[day_],ub_deaths_ts[day_]]
cum_D_at_peak = [median(cum_death_per_sim[:,day_]),quantile(cum_death_per_sim[:,day_],0.025),quantile(cum_death_per_sim[:,day_],0.975)]

median_H_ts = [ median(incidence_H[:,t]) for t in 1:365 ]
lb_H_ts = [ quantile(incidence_H[:,t],0.025) for t in 1:365 ]
ub_H_ts = [ quantile(incidence_H[:,t],0.975) for t in 1:365 ]
(max_H_inc,day_) = findmax(median_H_ts)
H_ts_at_peak = [median_H_ts[day_],lb_H_ts[day_],ub_H_ts[day_]]



median_Hocc_ts = [ median(hosp_occup_per_sim[:,t]) for t in 1:365 ]
lb_HOcc_ts = [ quantile(hosp_occup_per_sim[:,t],0.025) for t in 1:365 ]
ub_HOcc_ts = [ quantile(hosp_occup_per_sim[:,t],0.975) for t in 1:365 ]
(max_H_occ,day_) = findmax(median_Hocc_ts)
H_occ_ts_at_peak = [median_Hocc_ts[day_],lb_HOcc_ts[day_],ub_HOcc_ts[day_]]


median_ICUocc_ts = [ median(ICU_occup_per_sim[:,t]) for t in 1:365 ]
lb_ICUOcc_ts = [ quantile(ICU_occup_per_sim[:,t],0.025) for t in 1:365 ]
ub_ICUOcc_ts = [ quantile(ICU_occup_per_sim[:,t],0.975) for t in 1:365 ]
(max_ICU_occ,day_) = findmax(median_ICUocc_ts)
ICU_occ_ts_at_peak = [median_ICUocc_ts[day_],lb_ICUOcc_ts[day_],ub_ICUOcc_ts[day_]]

plt_incidence = plot(0:364,median_H_ts,
                        lab = "Hospitalisations",
                        lw = 3,color = :red,
                        ribbon = (median_H_ts.-lb_H_ts,ub_H_ts .- median_H_ts),
                        fillalpha = 0.25,
                        xticks = (first_of_months,["Apr","May","June","July","Aug","Sept","Oct","Nov","Dec"]),
                        xlims = (0.,275.),
                        title = "Incidence of disease")

 plot!(plt_incidence,0:364,median_deaths_ts,
                        lab = "Deaths",
                        lw = 3,color = :black,
                        ribbon = (median_deaths_ts.-lb_deaths_ts,ub_deaths_ts .- median_deaths_ts),
                        fillalpha = 0.5,
                        xticks = (first_of_months,["Apr","May","June","July","Aug","Sept","Oct","Nov","Dec"]),
                        xlims = (0.,275.),
                        ylabel = "Daily Incidence" )


plt_health_usage = plot(0:364,median_Hocc_ts,
                        lab = "Hospital beds",
                        lw = 3,color = :blue,
                        ribbon = (median_Hocc_ts.-lb_HOcc_ts,ub_HOcc_ts .- median_Hocc_ts),
                        fillalpha = 0.25,
                        xticks = (first_of_months,["Apr","May","June","July","Aug","Sept","Oct","Nov","Dec"]),
                        xlims = (0.,275.),
                        ylabel = "Daily Occupancy",
                        title = "Health system usage")
plot!(plt_health_usage,0:364,median_ICUocc_ts,
                        lab = "ICU beds",
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
savefig(plt_health_usage,"output/plt_health_usage.png")

# @save "output/plt_incidence.jld2" plt_incidence
# @save "output/plt_health_usage.jld2" plt_health_usage




y_up = [median(incidence[:,t]) for t = 1:365]
plot(y_up)
