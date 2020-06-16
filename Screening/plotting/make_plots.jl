push!(LOAD_PATH, "./src")
    push!(LOAD_PATH, "./Screening")
    using JLD2,Plots,Colors,ColorBrewer,Dates#DifferentialEquations,JLD2, CSV, DataFrames,Dates,Statistics
    #import KenyaCoV_screening

#######
fillalpha=.2
Plots.default(grid=false,legendfontsize=6,titlefontsize=10, xtickfontsize=10, ytickfontsize=9,fontfamily="Cambria",formatter=:scientific)
#StatsPlots.default(grid=false,legendfontsize=6, legendfontcolor=:gray30,titlefontsize=9, titlefontcolor=:gray30,xtickfontsize=9, ytickfontsize=9,tickfontcolor=:black,linewidth=false,xguidefontsize=9)
county_index_N=30;county_index_M=28
colors=[ColorBrewer.palette("Set2",7);ColorBrewer.palette("Set3",12)]
##
function make_plots(scenarios,output_folder)
    p_incidence=plot(title="Incidence in Kenya")
        #p_incidence_Nai=plot(title="Incidence in Nairobi")
        #p_incidence_Momb=plot(title="Incidence in Mombasa")
        #p_incidence_NaiMomb=plot(title="Incidence in Nairobi and Mombasa")
    p_incidence_V=plot(title="Incidence of severe cases in Kenya")

    p_H_ICU=plot(title="Hospital and ICU usage in Kenya")
    p_cases=bar(title="Cases in Kenya")
        cases_med=[];cases_lpred=[];cases_upred=[]
    p_cases_V=bar(title="Severe cases in Kenyaa")
        cases_V_med=[];cases_V_lpred=[];cases_V_upred=[]
    #p_cases_V_NaiMomb=bar(title="Severe cases in Nairobi and Mombasa")
        #cases_V_NM_med=[];cases_V_NM_lpred=[];cases_V_NM_upred=[]
    p_deaths=bar(title="Deaths in Kenya")
        deaths_med=[];deaths_lpred=[];deaths_upred=[]
    i=0
    for scenario in scenarios
        i+=1
        @load scenario[1]  results
        plot!(p_incidence,results.country_incidence_A_ts.med .+ results.country_incidence_M_ts.med .+ results.country_incidence_V_ts.med,
                ribbon=(results.country_incidence_A_ts.lpred .+ results.country_incidence_M_ts.lpred .+ results.country_incidence_V_ts.lpred,
                results.country_incidence_A_ts.upred .+ results.country_incidence_M_ts.upred .+ results.country_incidence_V_ts.upred),
                fillalpha=fillalpha,label=scenario[2],color=colors[i])
        #=plot!(p_incidence_Nai,results.incidence_A_ts.med[county_index_N,:] .+ results.incidence_M_ts.med[county_index_N,:] .+ results.incidence_V_ts.med[county_index_N,:],
                ribbon=(results.incidence_A_ts.lpred[county_index_N,:] .+ results.incidence_M_ts.lpred[county_index_N,:] .+ results.incidence_V_ts.lpred[county_index_N,:],
                results.incidence_A_ts.upred[county_index_N,:] .+ results.incidence_M_ts.upred[county_index_N,:] .+ results.incidence_V_ts.upred[county_index_N,:]),
                fillalpha=fillalpha,label=scenario[2])
        plot!(p_incidence_Momb,results.incidence_A_ts.med[county_index_M,:] .+ results.incidence_M_ts.med[county_index_M,:] .+ results.incidence_V_ts.med[county_index_M,:],
                ribbon=(results.incidence_A_ts.lpred[county_index_M,:] .+ results.incidence_M_ts.lpred[county_index_M,:] .+ results.incidence_V_ts.lpred[county_index_M,:],
                results.incidence_A_ts.upred[county_index_M,:] .+ results.incidence_M_ts.upred[county_index_M,:] .+ results.incidence_V_ts.upred[county_index_M,:]),
                fillalpha=fillalpha,label=scenario[2])
        plot!(p_incidence_NaiMomb,results.incidence_A_ts.med[county_index_N,:] .+ results.incidence_M_ts.med[county_index_N,:] .+ results.incidence_V_ts.med[county_index_N,:] .+results.incidence_A_ts.med[county_index_M,:] .+ results.incidence_M_ts.med[county_index_M,:] .+ results.incidence_V_ts.med[county_index_M,:],
                ribbon=(results.incidence_A_ts.lpred[county_index_N,:] .+ results.incidence_M_ts.lpred[county_index_N,:] .+ results.incidence_V_ts.lpred[county_index_N,:] .+ results.incidence_A_ts.lpred[county_index_M,:] .+ results.incidence_M_ts.lpred[county_index_M,:] .+ results.incidence_V_ts.lpred[county_index_M,:],
                results.incidence_A_ts.upred[county_index_N,:] .+ results.incidence_M_ts.upred[county_index_N,:] .+ results.incidence_V_ts.upred[county_index_N,:] .+ results.incidence_A_ts.upred[county_index_M,:] .+ results.incidence_M_ts.upred[county_index_M,:] .+ results.incidence_V_ts.upred[county_index_M,:]),
                fillalpha=fillalpha,label=scenario[2])=#
        plot!(p_incidence_V,results.country_incidence_V_ts.med,ribbon=(results.country_incidence_V_ts.lpred,results.country_incidence_V_ts.upred),fillalpha=fillalpha,label=scenario[2],color=colors[i])

        #plot!(p_H_ICU,results.prevalence_H_ts.med,ribbon=(results.prevalence_H_ts.lpred,results.prevalence_H_ts.upred),fillalpha=fillalpha,label=scenario[2])
        push!(cases_med,sum(results.total_cases_A_by_area_and_age.med)+sum(results.total_cases_M_by_area_and_age.med)+sum(results.total_cases_V_by_area_and_age.med))
            push!(cases_lpred,sum(results.total_cases_A_by_area_and_age.lpred)+sum(results.total_cases_M_by_area_and_age.lpred)+sum(results.total_cases_V_by_area_and_age.lpred))
            push!(cases_upred,sum(results.total_cases_A_by_area_and_age.upred)+sum(results.total_cases_M_by_area_and_age.upred)+sum(results.total_cases_V_by_area_and_age.upred))
        push!(cases_V_med,sum(results.total_cases_V_by_area_and_age.med))
            push!(cases_V_lpred,sum(results.total_cases_V_by_area_and_age.lpred))
            push!(cases_V_upred,sum(results.total_cases_V_by_area_and_age.upred))
        #=push!(cases_V_NM_med,sum(results.total_cases_V_by_area_and_age.med[county_index_N,:])+sum(results.total_cases_V_by_area_and_age.med[county_index_M,:]))
            push!(cases_V_NM_lpred,sum(results.total_cases_V_by_area_and_age.lpred[county_index_N,:])+sum(results.total_cases_V_by_area_and_age.lpred[county_index_M,:]))
            push!(cases_V_NM_upred,sum(results.total_cases_V_by_area_and_age.upred[county_index_N,:])+sum(results.total_cases_V_by_area_and_age.upred[county_index_M,:]))=#
        push!(deaths_med,results.total_deaths.med);push!(deaths_lpred,results.total_deaths.lpred);push!(deaths_upred,results.total_deaths.upred);
    end
    bar!(p_cases,[scenario[2] for scenario in scenarios],cases_med,errorbar=(cases_lpred,cases_upred),color=colors,legend=false)
        hline!(p_cases,[cases_med[1]])
    bar!(p_cases_V,[scenario[2] for scenario in scenarios],cases_V_med,errorbar=(cases_V_lpred,cases_V_upred),color=colors,legend=false)
        hline!(p_cases_V,[cases_V_med[1]])
    bar!(p_deaths,[scenario[2] for scenario in scenarios],deaths_med,errorbar=(deaths_lpred,deaths_upred),color=colors,legend=false)
        hline!(p_deaths,[deaths_med[1]])

    #saving plots
    if !isdir(output_folder)    mkdir(output_folder);   end
    savefig(p_incidence,output_folder*"incidence.png")
    savefig(p_incidence_V,output_folder*"incidence_V.png")
    savefig(p_cases,output_folder*"cases.png")
    savefig(p_cases_V,output_folder*"cases_V.png")
    savefig(p_deaths,output_folder*"deaths.png")
    #=display(p_incidence)#;  display(p_incidence_NaiMomb)#;  display(p_incidence_Momb);  display(p_incidence_Nai)
    display(p_incidence_V)
    display(p_cases);    display(p_cases_V)
    display(p_deaths)=#
end


######

scenarios=[("./Screening/output/0_NoInterv/session1_1000sims/noInterv_sc1_1000sims.jld2","No Intervention"),
         ("./Screening/output/1_CTH/session3_1000sims/CTH_sc1_1000sims.jld2","CTH (CTn=30)"),
            ("./Screening/output/1_CTH/session1_1000sims/CTH_sc1_1000sims.jld2","CTH (CTn=100)")#,
         #"./Screening/output/2_SympS/session1_1000sims/","./Screening/output/2_SympS/session2_1000sims/"#,
         #"./Screening/3_SympSCT/session3_1000sims/",#"./Screening/3_SympSCT/session4_1000sims/",
         #"./Screening/4_MS/session3_1000sims/","./Screening/4_MS/session4_1000sims/"#,
         #"./Screening/5_MSCT/session1_500sims/","./Screening/5_MSCT/session2_500sims/"
         ]
         scenarios=[("./Screening/output/0_NoInterv/session1_2sims/noInterv_sc1_2sims.jld2","No Intervention"),
                  ("./Screening/output/1_CTH/session3_2sims/CTH_sc1_2sims.jld2","CTH (CTn=30)"),
                     ("./Screening/output/1_CTH/session1_2sims/CTH_sc1_2sims.jld2","CTH (CTn=100)")#,
                  #"./Screening/output/2_SympS/session1_1000sims/","./Screening/output/2_SympS/session2_1000sims/"#,
                  #"./Screening/3_SympSCT/session3_1000sims/",#"./Screening/3_SympSCT/session4_1000sims/",
                  #"./Screening/4_MS/session3_1000sims/","./Screening/4_MS/session4_1000sims/"#,
                  #"./Screening/5_MSCT/session1_500sims/","./Screening/5_MSCT/session2_500sims/"
                  ]
    make_plots(scenarios,"Screening/plotting/plots__2020_"*Dates.format(now(), "mm_dd__HH_MM /"))
