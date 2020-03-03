riskregionnames = ["Machakos/Muranga" "Mandera" "Baringo/Nakuru" "Nairobi" "Uasin Gishu" "Marsabit" "Garissa" "Nakuru/Narok" "Turkana" "Taita Taveta" "Kitui/Meru" "Kilifi/Mombasa" "Kericho/Kisumu" "Kilifi/Lamu" "Kakamega/Kisumu" "Wajir" "Kajiado/Kisumu" "Homa bay/Migori" "Samburu/Laikipia" "Kilifi/Kwale" "Total"]
age_cats = ["0-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39","40-44","45-49","50-54","55-59","60-64","65-69","70-74","75+"]

# function plot_observed_incidence_timeseries(results,treatment_rates)
#     inc_D = results[1][1]
#     inc_A = results[1][2]
#     plt = plot(1:365,inc_A[21,:,1],
#                 ribbon = (inc_A[21,:,2],inc_A[21,:,3]),
#                 lab = "",
#                 fillalpha = 0.25,ls=:dash,lw=2,
#                 color = :black);
#     plot!(plt,1:365,inc_D[21,:,1],
#             # ribbon = (inc_D[21,:,2],inc_D[21,:,3]),
#             lab = "No treatment/isolation",
#             color = :black,
#             fillalpha = 0.25,lw=3);
#     inc_A = results[2][1]
#     inc_D = results[2][2]
#     plot!(plt,1:365,inc_A[21,:,1],
#                 ribbon = (inc_A[21,:,2],inc_A[21,:,3]),
#                 lab = "",
#                 fillalpha = 0.25,
#                 color = :red,ls=:dash,lw=2 );
#     plot!(plt,1:365,inc_D[21,:,1],
#             # ribbon = (inc_D[21,:,2],inc_D[21,:,3]),
#             lab = "Av. $(round(1/(7*treatment_rates[2]),digits = 0)) weeks to treatment",
#             fillalpha = 0.25,
#             color = :red,lw=3);
#
#     inc_A = results[3][1]
#     inc_D = results[3][2]
#     plot!(plt,1:365,inc_A[21,:,1],
#                 ribbon = (inc_A[21,:,2],inc_A[21,:,3]),
#                 lab = "",
#                 fillalpha = 0.25,
#                 color = :green,ls=:dash,lw=2 );
#     plot!(plt,1:365,inc_D[21,:,1],
#             # ribbon = (inc_D[21,:,2],inc_D[21,:,3]),
#             lab = "Av. $(round(1/(7*treatment_rates[3]),digits = 0)) weeks to treatment",
#             fillalpha = 0.25,
#             color = :green,lw=3);
#     inc_A = results[4][1]
#     inc_D = results[4][2]
#     plot!(plt,1:365,inc_A[21,:,1],
#                 ribbon = (inc_A[21,:,2],inc_A[21,:,3]),
#                 lab = "",
#                 fillalpha = 0.25,
#                 color = :blue,ls=:dash,lw=2 );
#     plot!(plt,1:365,inc_D[21,:,1],
#             # ribbon = (inc_D[21,:,2],inc_D[21,:,3]),
#             lab = "Av. $(round(1/(7*treatment_rates[4]),digits = 0)) weeks to treatment",
#             fillalpha = 0.25,
#             color = :blue,lw=3);
#
#     xlabel!(plt,"days")
#     ylabel!(plt,"Incidence")
#     title!(plt,"Scenario A: Total daily incidence")
#
#     return plt
# end

function plot_incidence_spatial(results,treatment_rates,i)
    τ = treatment_rates[i]
    inc_D = results[i][1]
    plt = plot(1:365,inc_D[4,:,1].+1,
                ribbon = (inc_D[4,:,2],inc_D[4,:,3]),
                fillalpha = 0.1,
                lab = riskregionnames[4],
                # lab = "",
                lw=4,
                xlims = (50.,100.),
                # yscale = :log10,
                legend = :topleft,
                legendfontsize = 9.1);
    for i in [15,12]
        if i != 4
            plot!(plt,1:365,inc_D[i,:,1].+1,
                    # lab = "",
                    ribbon = (inc_D[i,:,2],inc_D[i,:,3]),
                    lab = riskregionnames[i],
                    lw=4);
        end
    end
    xlabel!(plt,"Days")
    ylabel!(plt,"Daily incidence")
    # ylims!(plt,(1,1e5))
    if τ == 0
        title!(plt,"No isolation")
    else
        title!(plt," Mean time to isolation $(round(1/τ,digits = 1)) days")
    end

    return plt
end

function plot_total_incidence_by_age(results,treatment_rates)
    cases_age = results[1][4][21,:,:]
    plt = scatter(collect(1:16)*2,cases_age[:,1],
            yerror = (cases_age[:,2],cases_age[:,3]),
            lab = "No treatment/isolation",
            xticks = (collect(1:16)*2,age_cats))
    for i=2:4
        τ = treatment_rates[i]
        cases_age = results[i][4][21,:,:]
        scatter!(plt,collect(1:16)*2,cases_age[:,1],
                yerror = (cases_age[:,2],cases_age[:,3]),
                lab = "Av. $(round(1/(7*τ),digits = 0)) weeks to treatment")
    end
    title!(plt,"Scenario A: Total cases by age and treatment")
    ylabel!(plt,"Total cases")
    return plt

end

function plot_total_incidence_by_treatment(results,treatment_rates)
    cases= [sum(results[i][4][21,:,1]) for i = 1:4]
    case_L = [sqrt(sum(results[i][4][21,:,2].^2)) for i = 1:4]
    case_U = [sqrt(sum(results[i][4][21,:,3].^2)) for i = 1:4]

    plt = scatter(collect(1:4)*2,cases,
            # yerror = (case_L,case_U),
            lab = "",
            xticks = (collect(1:4)*2,["n/t","3 wks","2 wks","1 wk"]),
            ms = 10)
    title!(plt,"Scenario A: Summed median total cases by treatment")
    ylabel!(plt,"Total cases")
    xlabel!(plt,"Av. time to treatment")
    ylims!(plt,(0.,2.5e7))
    return plt
end
