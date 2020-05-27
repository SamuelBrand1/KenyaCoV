push!(LOAD_PATH, joinpath(homedir(),"GitHub/KenyaCoV/src"))
using JLD2,Plots,StatsPlots,CSV,DelimitedFiles,Dates
using Plots.PlotMeasures

using LinearAlgebra:eigen
using Statistics: median, quantile

@load "data/detection_rates_for_different_epsilons_model2.jld2" d_0 d_01 d_025 d_05 d_1
@load "data/posterior_distribution_R0.jld2"
χ_zhang = vcat(0.34*ones(3),ones(10),1.47*ones(4))

#Load age mixing matrices (these are all in to (row) from (col) format
@load "data/agemixingmatrix_Kenya_all_types.jld2" M_Kenya M_Kenya_ho M_Kenya_other M_Kenya_school M_Kenya_work
candidate_M_school = zeros(size(M_Kenya_school))
primary_M_school = zeros(size(M_Kenya_school))
secondary_M_school = zeros(size(M_Kenya_school)) #14-18
tertiary_M_school = zeros(size(M_Kenya_school)) #19-22


candidate_M_school[:,3] = 0.2*M_Kenya_school[:,3]#Standard 8
candidate_M_school[:,4] = 0.2*M_Kenya_school[:,4]#Form 4
primary_M_school[:,1:3] = M_Kenya_school[:,1:3]
secondary_M_school[:,3] = 0.2*M_Kenya_school[:,3]
secondary_M_school[:,4] = 0.8*M_Kenya_school[:,4]
tertiary_M_school[:,4] = 0.2*M_Kenya_school[:,4]
tertiary_M_school[:,5] = 0.6*M_Kenya_school[:,5]
heatmap(M_Kenya_ho)



twoweeks = 14.
regional_lkdown_start = (Date(2020,4,7) - Date(2020,3,13)).value
schools_open_june = (Date(2020,6,2) - Date(2020,3,13)).value
schools_close_august = (Date(2020,8,14) - Date(2020,3,13)).value
schools_open_august = (Date(2020,8,31) - Date(2020,3,13)).value
end_regional_lockdown = (Date(2020,5,16) - Date(2020,3,13)).value
schools_close_october = (Date(2020,10,30) - Date(2020,3,13)).value

schools_open_jan2021 = (Date(2021,1,4) - Date(2020,3,13)).value  # schools open 4th Jan 2021
schools_close_april2021 = (Date(2021,4,9) - Date(2020,3,13)).value  # Schools closed 9th April 2021
schools_open_may2021 = (Date(2021,5,3) - Date(2020,3,13)).value  # Schools open 3rd may 2021
schools_close_august2021 = (Date(2021,8,6) - Date(2020,3,13)).value  # Schools closed 6th August 2021
schools_open_august2021 = (Date(2021,8,30) - Date(2020,3,13)).value  # Schools open 30th August 2021
schools_close_october2021 = (Date(2021,10,22) - Date(2020,3,13)).value  # Schools closed 22nd October 2021


yearend = (Date(2021,12,31) - Date(2020,3,13)).value  # we now run intil Dec 2021

monthdates = [Date(2020,3,1) + Month(i) for i = 1:21 ]
monthnames = [monthname(d)[1:3]*"-$(year(d)-2000)" for d in monthdates]
tick_times = [(d - Date(2020,3,13)).value for d in monthdates]

"""
Set up parameters - ϵ = 1 scenario
"""

#Contact matrices at different points of epidemic
M_first_two_weeks = 0.8*M_Kenya_ho .+ 0.5*M_Kenya_school .+ M_Kenya_other .+ M_Kenya_work
M_schools_closed = 1.2*M_Kenya_ho  .+ 0.55*M_Kenya_other .+ 0.55*M_Kenya_work
M_schools_open_90perc = M_Kenya_ho  .+ 0.65*M_Kenya_other .+ 0.65*M_Kenya_work .+ 0.9*M_Kenya_school
M_schools_open_50perc = M_Kenya_ho  .+ 0.65*M_Kenya_other .+ 0.65*M_Kenya_work .+ 0.5*M_Kenya_school
M_schools_candidates_only_90perc = M_Kenya_ho  .+ 0.65*M_Kenya_other .+ 0.65*M_Kenya_work .+ 0.9*candidate_M_school
M_schools_candidates_only_50perc = M_Kenya_ho  .+ 0.65*M_Kenya_other .+ 0.65*M_Kenya_work .+ 0.5*candidate_M_school
M_schools_primary_only_90perc = M_Kenya_ho  .+ 0.65*M_Kenya_other .+ 0.65*M_Kenya_work .+ 0.9*primary_M_school
M_schools_primary_only_50perc = M_Kenya_ho  .+ 0.65*M_Kenya_other .+ 0.65*M_Kenya_work .+ 0.5*primary_M_school
M_schools_secondary_only_90perc = M_Kenya_ho  .+ 0.65*M_Kenya_other .+ 0.65*M_Kenya_work .+ 0.9*secondary_M_school
M_schools_secondary_only_50perc = M_Kenya_ho  .+ 0.65*M_Kenya_other .+ 0.65*M_Kenya_work .+ 0.5*secondary_M_school
M_schools_tertiary_only_90perc = M_Kenya_ho  .+ 0.65*M_Kenya_other .+ 0.65*M_Kenya_work .+ 0.9*tertiary_M_school
M_schools_tertiary_only_50perc = M_Kenya_ho  .+ 0.65*M_Kenya_other .+ 0.65*M_Kenya_work .+ 0.5*tertiary_M_school

M_array = [M_schools_closed,
            M_Kenya,
            M_first_two_weeks,
            M_schools_open_90perc,
            M_schools_open_50perc,
            M_schools_candidates_only_90perc,
            M_schools_candidates_only_50perc,
            M_schools_primary_only_90perc,
            M_schools_primary_only_50perc,
            M_schools_secondary_only_90perc,
            M_schools_secondary_only_50perc,
            M_schools_tertiary_only_90perc,
            M_schools_tertiary_only_50perc]
interventionnames = ["Schools closed",
                    "Unmitigated",
                    "Initial two weeks",
                    "Schools open (90% contacts)",
                    "Schools open (50% contacts)",
                    "Candidates only return (90% contacts)",
                    "Candidates only return (50% contacts)",
                    "Primary only return (90% contacts)",
                    "Primary only return (50% contacts)",
                    "Secondary only return (90% contacts)",
                    "Secondary only return (50% contacts)",
                    "Tertiary only return (90% contacts)",
                    "Tertiary only return (50% contacts)"]


#Set the susceptibility vector --- to specify the correct R₀ for closed scenario
sus_matrix = repeat(χ_zhang,1,17)
R_vector = [9. for a = 1:17]#This is only true for ϵ = 1, check otherwise
inf_matrix = repeat(R_vector',17,1)

eigs_kenya_closed, = eigen(sus_matrix.*M_schools_closed.*inf_matrix)
multiplier_for_kenya_closed = Real(eigs_kenya_closed[end])
χ = χ_zhang./multiplier_for_kenya_closed
fitted_sus_matrix = repeat(χ,1,17)

function get_multiplier(M)
    eigs, = eigen(fitted_sus_matrix.*M.*inf_matrix)
    return Real(eigs[end])
end

multipliers = map(get_multiplier,M_array)

fitted_median_R₀ = median(posterior_R₀)
fitted_lpred_R₀ = quantile(posterior_R₀,0.025)
fitted_upred_R₀ = quantile(posterior_R₀,0.975)

heatmap(M_Kenya_school)

plt_R₀ = bar(4:13,(multipliers.*fitted_median_R₀)[4:end],orientations = :horizonal,
        lab = "Includes SD",color = :green,
        yticks = (1:13,interventionnames),
        xlims = (0,3),
        title = "Rt under different interventions",
        xlabel = "Reproductive ratio")

bar!(plt_R₀,[2,3],[(multipliers.*fitted_median_R₀)[2:3]],orientations = :horizonal,
        color = :red,
        lab = "No SD")

bar!(plt_R₀,[(multipliers.*fitted_median_R₀)[1]],orientations = :horizonal,
            color = :blue,
            lab = "Fitted to cases")
plot!(plt_R₀,[(multipliers.*fitted_median_R₀)[1],(multipliers.*fitted_median_R₀)[1]],[0.,14],
        color = :black,lw =3 ,ls = :dot,lab = "")


scatter!(plt_R₀,multipliers.*fitted_median_R₀,1:13,
            xerror = (multipliers.*fitted_median_R₀ .- multipliers.*fitted_lpred_R₀,
                        multipliers.*fitted_upred_R₀ .- multipliers.*fitted_median_R₀),
            ms = 0.,color = :black,lab ="",lw=3)

savefig(plt_R₀,"plotting/Rt_under_different_interventions.png")

# baseline_lower_R₀_diff = baseline_median_R₀ - multiplier_for_kenya*quantile(KenyaCoV.d_R₀,0.025)
# baseline_higher_R₀_diff = multiplier_for_kenya*quantile(KenyaCoV.d_R₀,0.975) - baseline_median_R₀
#
# baseline_median_R₀_closed = multiplier_for_kenya_closed*quantile(KenyaCoV.d_R₀,0.5)
# baseline_lower_R₀_diff_closed = baseline_median_R₀_closed - multiplier_for_kenya_closed*quantile(KenyaCoV.d_R₀,0.025)
# baseline_higher_R₀_diff_closed = multiplier_for_kenya_closed*quantile(KenyaCoV.d_R₀,0.975) - baseline_median_R₀_closed
#
# baseline_median_R₀_open = multiplier_for_kenya_open*quantile(KenyaCoV.d_R₀,0.5)
# baseline_lower_R₀_diff_open = baseline_median_R₀_open - multiplier_for_kenya_open*quantile(KenyaCoV.d_R₀,0.025)
# baseline_higher_R₀_diff_open = multiplier_for_kenya_open*quantile(KenyaCoV.d_R₀,0.975) - baseline_median_R₀_open
#
#
#
#
# """
# Base line scenario
# """
#
#
#
# """
# SCENARIO 2 --- regional lockdown ending. Schools stay shut
# """
#
# times = 0:0.1:yearend
# median_Rt_sc2 = [baseline_median_R₀_closed*ramp_down(t) for t in times]
# lower_sc2 = [baseline_lower_R₀_diff_closed*ramp_down(t) for t in times]
# higher_sc2 = [baseline_higher_R₀_diff_closed*ramp_down(t) for t in times]
#
#
#
#
#
# """
# SCENARIO 3 --- Schools reopen in June
# """
# times_closed_until_june = 0:0.1:schools_open_june
# times_closed_until_august = 0:0.1:schools_open_august
# median_Rt_sc3 = [baseline_median_R₀_closed*ramp_down(t) for t in times_closed_until_june]
# lower_sc3 = [baseline_lower_R₀_diff_closed*ramp_down(t) for t in times_closed_until_june]
# higher_sc3 = [baseline_higher_R₀_diff_closed*ramp_down(t) for t in times_closed_until_june]
#
# median_Rt_sc4 = [baseline_median_R₀_closed*ramp_down(t) for t in times_closed_until_august]
# lower_sc4 = [baseline_lower_R₀_diff_closed*ramp_down(t) for t in times_closed_until_august]
# higher_sc4 = [baseline_higher_R₀_diff_closed*ramp_down(t) for t in times_closed_until_august]
#
# times_sc3 = [schools_open_june,schools_open_june,
#             schools_close_august,schools_close_august,
#             schools_open_august,schools_open_august,
#             schools_close_october,schools_close_october,
#             schools_open_jan2021,schools_open_jan2021,
#             schools_close_april2021,schools_close_april2021,
#             schools_open_may2021,schools_open_may2021,
#             schools_close_august2021,schools_close_august2021,
#             schools_open_august2021,schools_open_august2021,
#             schools_close_october2021,schools_close_october2021]
# median_Rt_sc3_close_and_open = repeat([baseline_median_R₀_closed*ramp_down(100.),baseline_median_R₀_open*ramp_down(100.),baseline_median_R₀_open*ramp_down(100.),baseline_median_R₀_closed*ramp_down(100.)],5,1)
# lower_Rt_sc3_open_and_close = repeat([baseline_lower_R₀_diff_closed*ramp_down(100.),baseline_lower_R₀_diff_open*ramp_down(100.),baseline_lower_R₀_diff_open*ramp_down(100.),baseline_lower_R₀_diff_closed*ramp_down(100.)],5,1)
# higher_Rt_sc3_open_and_close = repeat([baseline_higher_R₀_diff_closed*ramp_down(100.),baseline_higher_R₀_diff_open*ramp_down(100.),baseline_higher_R₀_diff_open*ramp_down(100.),baseline_higher_R₀_diff_closed*ramp_down(100.)],5,1)
#
#
# plt_Rt = plot([0.,yearend],[baseline_median_R₀,baseline_median_R₀],
#             lab = "No intervention",
#             ribbon = ([baseline_lower_R₀_diff,baseline_lower_R₀_diff],[baseline_higher_R₀_diff,baseline_higher_R₀_diff])
#             ,ylabel = "Rt",
#             xticks = (tick_times[1:2:end],monthnames[1:2:end]),
#             ylims = (0.,5.),lw = 3,
#             title = "Interventions and social distancing (30% contact reduction)")
#
#
# plot!(plt_Rt,times,median_Rt_sc2,
#         lab = "Schools closed, end regional lockdown in May",
#         ribbon = (lower_sc2,higher_sc2),lw=3)
# vcat(times_closed_until_june,times_sc3,[yearend])
# vcat(median_Rt_sc3,median_Rt_sc3_close_and_open,[baseline_median_R₀_closed*ramp_down(100.)])
#
# plot!(plt_Rt,vcat(times_closed_until_june,times_sc3,[yearend]),vcat(median_Rt_sc3,median_Rt_sc3_close_and_open,[baseline_median_R₀_closed*ramp_down(100.)]),
#         lab = "Schools open in June",lw = 3,
#         ribbon = (vcat(lower_sc3,lower_Rt_sc3_open_and_close,[baseline_lower_R₀_diff_closed*ramp_down(100.)])
#                     , vcat(higher_sc3,higher_Rt_sc3_open_and_close,[baseline_higher_R₀_diff_closed*ramp_down(100.)])))
#
#
# plot!(plt_Rt,vcat(times_closed_until_august,times_sc3[5:end],[yearend]),vcat(median_Rt_sc4,median_Rt_sc3_close_and_open[5:end],[baseline_median_R₀_closed*ramp_down(100.)]),
#     lab = "Schools open in August",lw = 3,
#     ribbon = (vcat(lower_sc4,lower_Rt_sc3_open_and_close[5:end],[baseline_lower_R₀_diff_closed*ramp_down(100.)])
#                 , vcat(higher_sc4,higher_Rt_sc3_open_and_close[5:end],[baseline_higher_R₀_diff_closed*ramp_down(100.)])))
#
# savefig(plt_Rt,"Rt_30_perc_reduction.png")
