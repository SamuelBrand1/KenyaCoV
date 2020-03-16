riskregionnames = ["Machakos/Muranga" "Mandera" "Baringo/Nakuru" "Nairobi" "Uasin Gishu" "Marsabit" "Garissa" "Nakuru/Narok" "Turkana" "Taita Taveta" "Kitui/Meru" "Kilifi/Mombasa" "Kericho/Kisumu" "Kilifi/Lamu" "Kakamega/Kisumu" "Wajir" "Kajiado/Kisumu" "Homa bay/Migori" "Samburu/Laikipia" "Kilifi/Kwale" "Total"]
age_cats = ["0-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39","40-44","45-49","50-54","55-59","60-64","65-69","70-74","75-79","80+"]

function plot_total_incidence(results_group,treatments::Tuple{Float64,Real},i)
    τ,ϵ_D = treatments
    inc_D1 = results_group[1][i][1]
    inc_D2 = results_group[2][i][1]
    inc_A2 = results_group[2][i][2]
    plt = plot(1:365,inc_D1[21,:,1].+1,
                fontfamily="Helvetica",
                lw = 3,
                legend = :topright,
                ribbon =(inc_D1[21,:,2],inc_D1[21,:,3]),
                # fillalpha = 0.15,
                # xlims = (0.,100),
                yscale = :log10,
                xlabel = "Days",
                ylabel = "Daily incidence + 1")
        # plot!(plt,1:365,inc_D2[21,:,1].+1,
        #             lw = 3,
        #             lab="MERS-like scenario",
        #             ribbon =(inc_D2[21,:,2],inc_D2[21,:,3]),
        #             fillalpha = 0.15)
        #             # yscale = :log10,
        #             # xlabel = "Days",
        #             # ylabel = "Daily incidence")
        # plot!(plt,1:365,inc_A2[21,:,1].+1,
        #             lw = 3,
        #             lab="MERS-like scenario: undetected",
        #             ls = :dot,
        #             # ribbon =(inc_D2[21,:,2],inc_D2[21,:,3]),
        #             fillalpha = 0.15)
        #             # yscale = :log10,
        #             # xlabel = "Days",
        #             # ylabel = "Daily incidence")

        # plot!(plt,1:365,inc_D3[21,:,1].+inc_A3[21,:,1].+1,
        #             lw = 3,
        #             lab="MERS-like scenario, delta = 0.8",
        #             # ribbon =(inc_D3[21,:,2],inc_D3[21,:,3]),
        #             fillalpha = 0.15,
        #             yscale = :log10,
        #             xlabel = "Days",
        #             ylabel = "Daily incidence")
    return plt
end

function plot_total_incidence_group(scenario_group,treatment_group,treatment_num,rel_transmission_perc)
        plt = plot()
        tracing_rate,e_D = treatment_group[treatment_num]
        for (i,scenarioresults) in enumerate(scenario_group)
                inc_D = scenarioresults[treatment_num][1]
                plot!(plt,1:365,inc_D[21,:,1].+1,
                fontfamily="Helvetica",
                lw = 3,
                lab= "rel. infect. of undetecteds: $(rel_transmission_perc[i])%",
                legend = :topright,
                ribbon =(inc_D[21,:,2],inc_D[21,:,3]),
                fillalpha = 0.15,
                # xlims = (0.,100),
                yscale = :log10,
                xlabel = "Days",
                ylabel = "Daily incidence + 1")
        end
        return plt
end

function plot_total_incidence(results)
    inc_D1 = results[1][1]

    plt = plot(1:365,inc_D1[21,:,1].+1,
                fontfamily="Helvetica",
                lw = 3,
                legend = :topright,
                ribbon =(inc_D1[21,:,2],inc_D1[21,:,3]),
                # fillalpha = 0.15,
                # xlims = (0.,100),
                yscale = :log10,
                xlabel = "Days",
                ylabel = "Daily incidence + 1")
        # plot!(plt,1:365,inc_D2[21,:,1].+1,
        #             lw = 3,
        #             lab="MERS-like scenario",
        #             ribbon =(inc_D2[21,:,2],inc_D2[21,:,3]),
        #             fillalpha = 0.15)
        #             # yscale = :log10,
        #             # xlabel = "Days",
        #             # ylabel = "Daily incidence")
        # plot!(plt,1:365,inc_A2[21,:,1].+1,
        #             lw = 3,
        #             lab="MERS-like scenario: undetected",
        #             ls = :dot,
        #             # ribbon =(inc_D2[21,:,2],inc_D2[21,:,3]),
        #             fillalpha = 0.15)
        #             # yscale = :log10,
        #             # xlabel = "Days",
        #             # ylabel = "Daily incidence")

        # plot!(plt,1:365,inc_D3[21,:,1].+inc_A3[21,:,1].+1,
        #             lw = 3,
        #             lab="MERS-like scenario, delta = 0.8",
        #             # ribbon =(inc_D3[21,:,2],inc_D3[21,:,3]),
        #             fillalpha = 0.15,
        #             yscale = :log10,
        #             xlabel = "Days",
        #             ylabel = "Daily incidence")
    return plt
end

function plot_total_incidence(results,treatment_rates)
    inc_D = results[1][1]
    plt = plot(1:365,inc_D[21,:,1].+1,
                fontfamily="Helvetica",
                lw = 3,
                lab="No isolation",
                legend = :topleft,
                ribbon =(inc_D[21,:,2],inc_D[21,:,3]),
                fillalpha = 0.15,
                xlims = (30.,100),
                yscale = :log10,
                xlabel = "Days",
                ylabel = "Daily incidence")
    # for i = 2:3
    #     τ = treatment_rates[i]
    #     inc_D = results[i][1]
    #     plot!(plt,1:365,inc_D[21,:,1].+1,
    #                 lw = 3,
    #                 lab="Isolation in mean $(round(1/τ,digits = 1)) days",
    #                 ribbon =(inc_D[21,:,2],inc_D[21,:,3]),
    #                 fillalpha = 0.15,
    #                 xlims = (30.,100),
    #                 # yscale = :log10,
    #                 xlabel = "Days",
    #                 ylabel = "Daily incidence")
    # end
    return plt
end


function plot_incidence_spatial(results,treatment_rates,i,ordering)
    τ = treatment_rates[i]
    inc_D = results[i][1]
    plt = plot(1:365,inc_D[ordering[1],:,1].+1,
                # ribbon = (inc_D[4,:,2],inc_D[4,:,3]),
                fillalpha = 0.1,
                lab = riskregionnames[4],
                # lab = "",
                lw=2,
                xlims = (0.,100.),
                yscale = :log10,
                legend = :outertopright,
                legendfontsize = 8.9);
    for i in ordering[2:end]
        if i != 4
            plot!(plt,1:365,inc_D[i,:,1].+1,
                    # lab = "",
                    # ribbon = (inc_D[i,:,2],inc_D[i,:,3]),
                    lab = riskregionnames[i],
                    lw=2);
        end
    end
    xlabel!(plt,"Days")
    ylabel!(plt,"Daily incidence")
    # ylims!(plt,(1,1e5))
    # if τ == 0
    #     title!(plt,"No isolation")
    # else
    #     title!(plt," Mean time to isolation $(round(1/τ,digits = 1)) days")
    # end

    return plt
end

function plot_total_incidence_by_age(results,treatment_rates,i)
    d1,d2,d3 = size(results[1][4][:,:,:])
    cases_age = zeros(d3,d2)
    for k = 1:d3,a = 1:d2
        cases_age[k,a] = sum(results[i][4][:,a,k])
    end
    plt = boxplot(cases_age,lab = "",
                xticks = (collect(1:2:16),age_cats[1:2:16]),color = :blue,
                ylabel = "Detected cases",
                xlabel = "Age")
    τ = treatment_rates[i]
    if τ > 0
        title!(plt,"Mean time to isolation $(round(1/τ,digits = 1)) days")
    else
        title!(plt,"No isolation")
    end


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
