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


function plot_incidence_spatial(results,treatment_rates,i)
    ordering = sortperm(results[i][1][1:20,100,1],rev = true)
    τ = treatment_rates[i]
    inc_D = results[i][1]
    plt = plot()
    for i in ordering
        plot!(plt,1:365,inc_D[i,:,1].+1,
                    # ribbon = (inc_D[4,:,2],inc_D[4,:,3]),
                    fillalpha = 0.1,
                    lab = riskregionnames[i],
                    lw=2,
                    xlims = (0.,100.),
                    yscale = :log10,
                    legend = :outertopright,
                    legendfontsize = 8.9);
    end
    xlabel!(plt,"Days")
    ylabel!(plt,"Daily incidence")

    return plt
end

function plot_total_incidence_by_age(scenario_group,treatment_rates,rel_transmission_perc,i)
    d1,d2,d3 = size(scenario_group[1][i][4][:,:,:])
    n = length(scenario_group)
    cases_age = zeros(d2,n,3)
    label_scen = "$(rel_transmission_perc[1])% rel. infectiousness"
    for l in 2:length(rel_transmission_perc)
        label_scen = hcat(label_scen,"$(rel_transmission_perc[l])% rel. infectiousness")
    end
    for (n,results) in enumerate(scenario_group)
        for a = 1:d2
            cases_age[a,n,1] = median([sum(results[i][4][:,a,k]) for k =1:1000])
            cases_age[a,n,2] = quantile([sum(results[i][4][:,a,k]) for k =1:1000],0.025)
            cases_age[a,n,3] = quantile([sum(results[i][4][:,a,k]) for k =1:1000],0.975)
        end
    end
    plt = groupedbar(cases_age[:,:,1]./1e6,
                    fontfamily="Helvetica",
                    lab =label_scen,
                    xticks = (1:17,age_cats))

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
