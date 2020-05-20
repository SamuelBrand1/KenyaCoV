push!(LOAD_PATH, joinpath(homedir(),"GitHub/KenyaCoV/src"))
using Plots,Parameters,Distributions,DifferentialEquations,JLD2,DataFrames,StatsPlots,FileIO,MAT,RecursiveArrayTools,CSV
using Revise
import KenyaCoV
using LinearAlgebra:eigen
using Statistics: median, quantile

"""
Consensus modelling --- May week 1

1) estimate epidemic spread by county
2) estimate the effect of reopening schools on either 2 June 2020 or delaying opening to 31 August 2020
3) the effect of lifting travel restrictions between counties in [16th] May 2020. 
"""

"""
Load age structured data, define callback control measures, and effect of regional lockdown
"""




@load "data/detection_rates_for_different_epsilons_model2.jld2" d_0 d_01 d_025 d_05 d_1
χ_zhang = vcat(0.34*ones(3),ones(10),1.47*ones(4))

#Initial infecteds
"""
Load age mixing matrices (these are all in to (row) from (col) format)
"""

@load "data/agemixingmatrix_Kenya_all_types.jld2" M_Kenya M_Kenya_ho M_Kenya_other M_Kenya_school M_Kenya_work
@load "data/agemixingmatrix_china.jld2" M_China

#Function for changing contacts so as to have -45% over 14 days
function ramp_down(t)
    if t < 0.
        return 1.
    elseif t>= 0. && t <= 14.
        return (1-t/14.) + 0.70*t/14.
    elseif t > 14.
        return 0.70
    end
end
using Dates
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
Set up parameters
"""



u0,P,P_dest = KenyaCoV.model_ingredients_from_data("data/data_for_age_structuredmodel_with_counties.jld2",
                                            "data/flight_numbers.csv",
                                            "data/projected_global_prevelance.csv")
#Redistribute susceptibility JUST to rescale infectiousness so we get the correct R₀/r
P.χ = copy(χ_zhang)
P.rel_detection_rate = d_1
P.dt = 0.25;
P.ext_inf_rate = 0.;
P.ϵ = 1.


#Set the susceptibility vector --- just to specify the correct R₀
sus_matrix = repeat(χ_zhang,1,17)
R_A = P.ϵ*((1/P.σ₂) + (1/P.γ) ) #effective duration of asymptomatic
R_M = (P.ϵ/P.σ₂) + (P.ϵ_D/P.γ) #effective duration of mild
R_V = (P.ϵ/P.σ₂) + (P.ϵ_V/P.τ) #effective duration of severe
R_vector = [(1-P.rel_detection_rate[a])*R_A + P.rel_detection_rate[a]*(1-P.hₐ[a])*R_M + P.rel_detection_rate[a]*P.hₐ[a]*R_V for a = 1:17]
inf_matrix = repeat(R_vector',17,1)

eigs_china, = eigen(sus_matrix.*M_China.*inf_matrix)
max_eigval_china = Real(eigs_china[end])
eigs_kenya, = eigen(sus_matrix.*M_Kenya.*inf_matrix)
max_eigval_Kenya = Real(eigs_kenya[end])
multiplier_for_kenya = max_eigval_Kenya/max_eigval_china
P.χ .= χ_zhang ./max_eigval_china #This rescales everything so β is the same as R₀ for China

sus_matrix = repeat(P.χ,1,17)
M_closed = 1.2*M_Kenya_ho .+ M_Kenya_other .+ M_Kenya_work

eigs_kenya_closed, = eigen(sus_matrix.*M_closed.*inf_matrix)
multiplier_for_kenya_closed = Real(eigs_kenya_closed[end])
M_open = 1.1*M_Kenya_ho .+ M_Kenya_other .+ M_Kenya_work .+ 0.7*M_Kenya_school
eigs_kenya_open, = eigen(sus_matrix.*M_open.*inf_matrix)
multiplier_for_kenya_open = Real(eigs_kenya_open[end])

baseline_median_R₀ = multiplier_for_kenya*quantile(KenyaCoV.d_R₀,0.5)
baseline_lower_R₀_diff = baseline_median_R₀ - multiplier_for_kenya*quantile(KenyaCoV.d_R₀,0.025)
baseline_higher_R₀_diff = multiplier_for_kenya*quantile(KenyaCoV.d_R₀,0.975) - baseline_median_R₀

baseline_median_R₀_closed = multiplier_for_kenya_closed*quantile(KenyaCoV.d_R₀,0.5)
baseline_lower_R₀_diff_closed = baseline_median_R₀_closed - multiplier_for_kenya_closed*quantile(KenyaCoV.d_R₀,0.025)
baseline_higher_R₀_diff_closed = multiplier_for_kenya_closed*quantile(KenyaCoV.d_R₀,0.975) - baseline_median_R₀_closed

baseline_median_R₀_open = multiplier_for_kenya_open*quantile(KenyaCoV.d_R₀,0.5)
baseline_lower_R₀_diff_open = baseline_median_R₀_open - multiplier_for_kenya_open*quantile(KenyaCoV.d_R₀,0.025)
baseline_higher_R₀_diff_open = multiplier_for_kenya_open*quantile(KenyaCoV.d_R₀,0.975) - baseline_median_R₀_open




"""
Base line scenario
"""



"""
SCENARIO 2 --- regional lockdown ending. Schools stay shut
"""

times = 0:0.1:yearend
median_Rt_sc2 = [baseline_median_R₀_closed*ramp_down(t) for t in times]
lower_sc2 = [baseline_lower_R₀_diff_closed*ramp_down(t) for t in times]
higher_sc2 = [baseline_higher_R₀_diff_closed*ramp_down(t) for t in times]





"""
SCENARIO 3 --- Schools reopen in June
"""
times_closed_until_june = 0:0.1:schools_open_june
times_closed_until_august = 0:0.1:schools_open_august
median_Rt_sc3 = [baseline_median_R₀_closed*ramp_down(t) for t in times_closed_until_june]
lower_sc3 = [baseline_lower_R₀_diff_closed*ramp_down(t) for t in times_closed_until_june]
higher_sc3 = [baseline_higher_R₀_diff_closed*ramp_down(t) for t in times_closed_until_june]

median_Rt_sc4 = [baseline_median_R₀_closed*ramp_down(t) for t in times_closed_until_august]
lower_sc4 = [baseline_lower_R₀_diff_closed*ramp_down(t) for t in times_closed_until_august]
higher_sc4 = [baseline_higher_R₀_diff_closed*ramp_down(t) for t in times_closed_until_august]

times_sc3 = [schools_open_june,schools_open_june,
            schools_close_august,schools_close_august,
            schools_open_august,schools_open_august,
            schools_close_october,schools_close_october,
            schools_open_jan2021,schools_open_jan2021,
            schools_close_april2021,schools_close_april2021,
            schools_open_may2021,schools_open_may2021,
            schools_close_august2021,schools_close_august2021,
            schools_open_august2021,schools_open_august2021,
            schools_close_october2021,schools_close_october2021]
median_Rt_sc3_close_and_open = repeat([baseline_median_R₀_closed*ramp_down(100.),baseline_median_R₀_open*ramp_down(100.),baseline_median_R₀_open*ramp_down(100.),baseline_median_R₀_closed*ramp_down(100.)],5,1)
lower_Rt_sc3_open_and_close = repeat([baseline_lower_R₀_diff_closed*ramp_down(100.),baseline_lower_R₀_diff_open*ramp_down(100.),baseline_lower_R₀_diff_open*ramp_down(100.),baseline_lower_R₀_diff_closed*ramp_down(100.)],5,1)
higher_Rt_sc3_open_and_close = repeat([baseline_higher_R₀_diff_closed*ramp_down(100.),baseline_higher_R₀_diff_open*ramp_down(100.),baseline_higher_R₀_diff_open*ramp_down(100.),baseline_higher_R₀_diff_closed*ramp_down(100.)],5,1)


plt_Rt = plot([0.,yearend],[baseline_median_R₀,baseline_median_R₀],
            lab = "No intervention",
            ribbon = ([baseline_lower_R₀_diff,baseline_lower_R₀_diff],[baseline_higher_R₀_diff,baseline_higher_R₀_diff])
            ,ylabel = "Rt",
            xticks = (tick_times[1:2:end],monthnames[1:2:end]),
            ylims = (0.,5.),lw = 3,
            title = "Interventions and social distancing (30% contact reduction)")


plot!(plt_Rt,times,median_Rt_sc2,
        lab = "Schools closed, end regional lockdown in May",
        ribbon = (lower_sc2,higher_sc2),lw=3)
vcat(times_closed_until_june,times_sc3,[yearend])
vcat(median_Rt_sc3,median_Rt_sc3_close_and_open,[baseline_median_R₀_closed*ramp_down(100.)])

plot!(plt_Rt,vcat(times_closed_until_june,times_sc3,[yearend]),vcat(median_Rt_sc3,median_Rt_sc3_close_and_open,[baseline_median_R₀_closed*ramp_down(100.)]),
        lab = "Schools open in June",lw = 3,
        ribbon = (vcat(lower_sc3,lower_Rt_sc3_open_and_close,[baseline_lower_R₀_diff_closed*ramp_down(100.)])
                    , vcat(higher_sc3,higher_Rt_sc3_open_and_close,[baseline_higher_R₀_diff_closed*ramp_down(100.)])))


plot!(plt_Rt,vcat(times_closed_until_august,times_sc3[5:end],[yearend]),vcat(median_Rt_sc4,median_Rt_sc3_close_and_open[5:end],[baseline_median_R₀_closed*ramp_down(100.)]),
    lab = "Schools open in August",lw = 3,
    ribbon = (vcat(lower_sc4,lower_Rt_sc3_open_and_close[5:end],[baseline_lower_R₀_diff_closed*ramp_down(100.)])
                , vcat(higher_sc4,higher_Rt_sc3_open_and_close[5:end],[baseline_higher_R₀_diff_closed*ramp_down(100.)])))

savefig(plt_Rt,"Rt_30_perc_reduction.png")
