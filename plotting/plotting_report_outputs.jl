push!(LOAD_PATH, "/Users/Sam/GitHub/KenyaCoV/src")
using DataFrames,CSV,MAT,Statistics,LinearAlgebra,Optim,Plots,JLD2,RData,Distances,DelimitedFiles,Dates

county_populations = CSV.read("data/2019_census_age_pyramids_counties.csv")
N_pop = zeros(Int64,47,17)
for i = 1:47, a = 1:17
    N_pop[i,a] = county_populations[i,a+1]
end

## Plot Verity prediction of deaths in Mombasa after 54k infections (uniform attack rate)

age_cats = ["0-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39","40-44","45-49","50-54","55-59","60-64","65-69","70-74","75-79","80+"]
f = findfirst(county_populations.county .== "Mombasa")
N_mombasa = N_pop[f,:]
bar(N_mombasa,orientation=:horizonal,yticks = (1:17,age_cats),lab = "")
verity_IFR = [0.00161,0.00695,0.0309,
            0.0844,
            0.161,
            0.595,
            1.93,
            4.28,
            7.8]./100

verity_comparison = zeros(17)
for a = 1:8
    verity_comparison[2*(a-1)+1] = verity_IFR[a]
    verity_comparison[2*(a-1)+2] = verity_IFR[a]
end
verity_comparison[17] = verity_IFR[9]
plot(verity_comparison)
uniform_cases_mombasa = 54000*N_mombasa/sum(N_mombasa)
pred_deaths = verity_comparison.*uniform_cases_mombasa
plt_pred_deaths = bar(pred_deaths,orientation=:horizonal,
                        yticks = (1:17,age_cats),
                        xticks = 1:15,
                        lab = "",
                        title = "Predicted deaths in Mombasa after 54k cases. Total = $(round(sum(pred_deaths),digits = 1))",
                        xlabel = "Deaths per age group",
                        size = (1000,600))
sum(pred_deaths)
savefig(plt_pred_deaths,"predicted_deaths_mombasa.png")

## Plot the predicted times series for Mombasa under different scenarios

@load "reports/report_unmitigated/scenario_data_unmitigated.jld2"
scenariodata_unmitigated = deepcopy(scenariodata)
@load "reports/report_full_intervention/scenario_data_full_intervention.jld2"
scenariodata_fullmitigation = deepcopy(scenariodata)

t_today = (Date(2020,6,3) - Date(2020,3,13)).value

V_unmitigated = cumsum(scenariodata_unmitigated.incidence_V_ts.med,dims = 2)[f,:]
V_unmitigated_lpred = cumsum(scenariodata_unmitigated.incidence_V_ts.lpred,dims = 2)[f,:]
V_unmitigated_upred = cumsum(scenariodata_unmitigated.incidence_V_ts.upred,dims = 2)[f,:]
MV_unmitigated = cumsum(scenariodata_unmitigated.incidence_V_ts.med,dims = 2)[f,:] .+ cumsum(scenariodata_unmitigated.incidence_M_ts.med,dims = 2)[f,:]
MV_unmitigated_lpred = cumsum(scenariodata_unmitigated.incidence_V_ts.lpred,dims = 2)[f,:] .+ cumsum(scenariodata_unmitigated.incidence_M_ts.lpred,dims = 2)[f,:]
MV_unmitigated_upred = cumsum(scenariodata_unmitigated.incidence_V_ts.upred,dims = 2)[f,:] .+ cumsum(scenariodata_unmitigated.incidence_M_ts.upred,dims = 2)[f,:]

V_fullmitigation = cumsum(scenariodata_fullmitigation.incidence_V_ts.med,dims = 2)[f,:]
V_fullmitigation_lpred = cumsum(scenariodata_fullmitigation.incidence_V_ts.lpred,dims = 2)[f,:]
V_fullmitigation_upred = cumsum(scenariodata_fullmitigation.incidence_V_ts.upred,dims = 2)[f,:]
MV_fullmitigation = cumsum(scenariodata_fullmitigation.incidence_V_ts.med,dims = 2)[f,:] .+ cumsum(scenariodata_fullmitigation.incidence_M_ts.med,dims = 2)[f,:]
MV_fullmitigation_lpred = cumsum(scenariodata_fullmitigation.incidence_V_ts.lpred,dims = 2)[f,:] .+ cumsum(scenariodata_fullmitigation.incidence_M_ts.lpred,dims = 2)[f,:]
MV_fullmitigation_upred = cumsum(scenariodata_fullmitigation.incidence_V_ts.upred,dims = 2)[f,:] .+ cumsum(scenariodata_fullmitigation.incidence_M_ts.upred,dims = 2)[f,:]


plt_mombasa = plot(V_unmitigated.+1,
                    lab = "Cum. severe cases (unmitigated)",
                    title = "Mombasa: Mild and severe case with detected cases so far",
                    lw = 2,
                    yscale = :log10,
                    size = (1000,600),
                    legend = :bottomright,
                    ribbon = (V_unmitigated .- V_unmitigated_lpred,
                                    V_unmitigated_upred.- V_unmitigated ),
                        xlabel = "Days since 13th March",
                        ylabel = "Cum. infections + 1")
plot!(plt_mombasa,MV_unmitigated .+1,
        lw = 2,
        lab = "Cum. mild and severe cases (unmitigated)",
        ribbon = (MV_unmitigated .- MV_unmitigated_lpred,
                        MV_unmitigated_upred.- MV_unmitigated ))

plot!(plt_mombasa,V_fullmitigation.+1,
        lab = "Cum. severe cases (full mitigation)",
        lw = 2,
        ribbon = (V_fullmitigation .- V_fullmitigation_lpred,
                        V_fullmitigation_upred.- V_fullmitigation ))
plot!(plt_mombasa,MV_fullmitigation.+1,
        lab = "Cum. mild and severe cases (full mitigation)",
        lw = 2,
        ribbon = (MV_fullmitigation .- MV_fullmitigation_lpred,
                        MV_fullmitigation_upred.- MV_fullmitigation ))

scatter!([t_today],[571 + 1],ms = 10,lab = "Number of confirmed cases",color = :black)

pop_mombasa = sum(N_mombasa)
seropos_unmitigated = cumsum(scenariodata_unmitigated.incidence_A_ts.med,dims = 2)[f,:] .- cumsum(scenariodata_unmitigated.incidence_M_ts.med,dims = 2)[f,:] .- cumsum(scenariodata_unmitigated.incidence_V_ts.med,dims = 2)[f,:]
seropos_unmitigated_lpred = cumsum(scenariodata_unmitigated.incidence_A_ts.lpred,dims = 2)[f,:] .- cumsum(scenariodata_unmitigated.incidence_M_ts.lpred,dims = 2)[f,:] .- cumsum(scenariodata_unmitigated.incidence_V_ts.lpred,dims = 2)[f,:]
seropos_unmitigated_upred = cumsum(scenariodata_unmitigated.incidence_A_ts.upred,dims = 2)[f,:] .- cumsum(scenariodata_unmitigated.incidence_M_ts.upred,dims = 2)[f,:] .- cumsum(scenariodata_unmitigated.incidence_V_ts.upred,dims = 2)[f,:]
seropos_unmitigated = seropos_unmitigated./pop_mombasa
seropos_unmitigated_lpred = seropos_unmitigated_lpred./pop_mombasa
seropos_unmitigated_upred = seropos_unmitigated_upred./pop_mombasa
seropos_fullmitigation = cumsum(scenariodata_fullmitigation.incidence_A_ts.med,dims = 2)[f,:] .- cumsum(scenariodata_fullmitigation.incidence_M_ts.med,dims = 2)[f,:] .- cumsum(scenariodata_fullmitigation.incidence_V_ts.med,dims = 2)[f,:]
seropos_fullmitigation_lpred = cumsum(scenariodata_fullmitigation.incidence_A_ts.lpred,dims = 2)[f,:] .- cumsum(scenariodata_fullmitigation.incidence_M_ts.lpred,dims = 2)[f,:] .- cumsum(scenariodata_fullmitigation.incidence_V_ts.lpred,dims = 2)[f,:]
seropos_fullmitigation_upred = cumsum(scenariodata_fullmitigation.incidence_A_ts.upred,dims = 2)[f,:] .- cumsum(scenariodata_fullmitigation.incidence_M_ts.upred,dims = 2)[f,:] .- cumsum(scenariodata_fullmitigation.incidence_V_ts.upred,dims = 2)[f,:]
seropos_fullmitigation = seropos_fullmitigation./pop_mombasa
seropos_fullmitigation_lpred = seropos_fullmitigation_lpred./pop_mombasa
seropos_fullmitigation_upred = seropos_fullmitigation_upred./pop_mombasa

times = collect(11:(length(seropos_unmitigated)+10))

plt_seroprev = plot(times,seropos_unmitigated,
                ribbon = (seropos_unmitigated - seropos_unmitigated_lpred, seropos_unmitigated_upred - seropos_unmitigated),
                lab = "Unmitigated (assumed fixed 10 days delay to seroconversion)",
                legend = :topleft,
                xlabel = "Days since 13th March",
                ylabel = "Proportion of pop. seropositive",
                title = "Mombasa seroconversion")
plot!(plt_seroprev,times,seropos_fullmitigation,
                ribbon = (seropos_fullmitigation - seropos_fullmitigation_lpred, seropos_fullmitigation_upred - seropos_fullmitigation),
                lab = "Full mitigation (assumed fixed 10 days delay to seroconversion)")
scatter!(plt_seroprev,[t_today - 7.],[4.5/100],lab = "Seroprevalence estimate (assumed 7 days ago)",
        ms = 5)
savefig(plt_mombasa,"mombasa_incidence.png")

savefig(plt_seroprev,"mombasa_seroprev.png")
