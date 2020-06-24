push!(LOAD_PATH, "./src")
    push!(LOAD_PATH, "./Screening")
    using JLD2,Plots,Colors,ColorBrewer,ColorSchemes,Dates,CSV,DataFrames
    #import KenyaCoV_screening

#######
Plots.default(grid=false,legendfontsize=6,titlefontsize=10, xtickfontsize=9, ytickfontsize=9,fontfamily="Cambria",labelfontsize=9,formatter=:scientific)
    fillalpha=.2
    county_index_N=30;county_index_M=28
    colors=[:red;ColorBrewer.palette("Set2",7)[1:5];ColorBrewer.palette("Set3",12)]
    colors=[ColorBrewer.palette("Set2",8);ColorBrewer.palette("Accent",8)]
    #color_palette=:seaborn_bright
    colors=[ColorSchemes.seaborn_bright.colors[i]   for i=1:size(ColorSchemes.seaborn_bright.colors,1)]
        c=colors[1];colors[1]=colors[4];colors[4]=c;
    size_x=250;size_y=350;xrotation=60;
    S0= 4.7562088e7;
##
function make_plots(scenarios,output_folder)
    @load scenarios[end][1]  results
    sim_dur=size(results.country_incidence_A_ts.med,1)
    xticks=([i  for i=1:sim_dur    if i%30==0],[string(Int(i/30))  for i=1:sim_dur    if i%30==0])
    colors_cp=deepcopy(colors)
    if size(scenarios,1)>size(colors_cp,1)
        colors_cp=repeat(colors_cp,Int(ceil(size(scenarios,1)/size(colors_cp,1))))
    end
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
                fillalpha=fillalpha,label=scenario[2],color=colors_cp[i],xticks=xticks,xlabel="months")#,xlims=(2*30,6*))
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
        plot!(p_incidence_V,results.country_incidence_V_ts.med,ribbon=(results.country_incidence_V_ts.lpred,results.country_incidence_V_ts.upred),
                fillalpha=fillalpha,label=scenario[2],color=colors_cp[i],xticks=xticks,xlabel="months")

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
    #redoing the medians
    i=0
    for scenario in scenarios
        i+=1
        @load scenario[1]  results
        plot!(p_incidence,results.country_incidence_A_ts.med .+ results.country_incidence_M_ts.med .+ results.country_incidence_V_ts.med,
                label=false,color=colors_cp[i])#,xticks=xticks,xlabel="months")
    end

    bar!(p_cases,[scenario[2] for scenario in scenarios],cases_med,errorbar=(cases_lpred,cases_upred),
        color=colors_cp,legend=false,linecolor=:transparent,markerstrokecolor=:black,markerstrokewidth=.7,xrotation=xrotation,wide=:true)
        hline!(p_cases,[cases_med[1]],color=colors_cp[1])
    bar!(p_cases_V,[scenario[2] for scenario in scenarios],cases_V_med,errorbar=(cases_V_lpred,cases_V_upred),color=colors_cp,legend=false)
        hline!(p_cases_V,[cases_V_med[1]],color=colors_cp[1])
    bar!(p_deaths,[scenario[2] for scenario in scenarios],deaths_med,errorbar=(deaths_lpred,deaths_upred),
        color=colors_cp,legend=false,linecolor=:transparent,markerstrokecolor=:black,markerstrokewidth=.7,xrotation=xrotation,wide=:true)
        hline!(p_deaths,[deaths_med[1]],color=colors_cp[1])

    #saving plots
    if !isdir(output_folder)    mkdir(output_folder);   end
    savefig(p_incidence,output_folder*"incidence.png");savefig(plot!(p_incidence,xlims=(2*30+15,18*30),size=(310,320)),output_folder*"incidence_2.png")
    savefig(p_incidence_V,output_folder*"incidence_V.png")
    savefig(p_cases,output_folder*"cases.png");savefig(plot!(p_cases,size=(250,380)),output_folder*"cases_2.png")
    savefig(p_cases_V,output_folder*"cases_V.png")
    savefig(p_deaths,output_folder*"deaths.png");savefig(plot!(p_deaths,size=(310,380)),output_folder*"deaths_2.png")
    #display(p_incidence)#;  display(p_incidence_NaiMomb)#;  display(p_incidence_Momb);  display(p_incidence_Nai)
    #display(p_incidence_V)
    #display(p_cases);
    #display(p_cases_V)
    #display(p_deaths)

    df=DataFrame(files=[scenario[1] for scenario in scenarios],
                labels=[scenario[2] for scenario in scenarios],
                cases_fiff=[cases_med[1]-cases_med[i]    for i=1:size(cases_med,1)],cases_med=cases_med,cases_lpred=cases_lpred,cases_upred=cases_upred,
                deaths_fiff=[deaths_med[1]-deaths_med[i]    for i=1:size(deaths_med,1)],deaths_med=deaths_med,deaths_lpred=deaths_lpred,deaths_upred=deaths_upred)
    CSV.write(output_folder*"2020_"*Dates.format(now(), "mm_dd__HH_MM") *".csv",df)
    return cases_med,cases_lpred,cases_upred,deaths_med,deaths_lpred,deaths_upred
end

function make_plots_withTOTALCASES(scenarios,output_folder)
    @load scenarios[1][1]  results
    sim_dur=size(results.country_incidence_A_ts.med,1)
    xticks=([i  for i=1:sim_dur    if i%30==0],[string(Int(i/30))  for i=1:sim_dur    if i%30==0])
    colors_cp=deepcopy(colors)
    if size(scenarios,1)>size(colors_cp,1)
        colors_cp=repeat(colors_cp,Int(ceil(size(scenarios,1)/size(colors_cp,1))))
    end
    p_incidence=plot(title="Incidence in Kenya")
    p_incidence_V=plot(title="Incidence of severe cases in Kenya")

    p_H_ICU=plot(title="Hospital and ICU usage in Kenya")
    p_cases=bar(title="Cases in Kenya")
        cases_med=[];cases_lpred=[];cases_upred=[]
    p_cases_V=bar(title="Severe cases in Kenyaa")
        cases_V_med=[];cases_V_lpred=[];cases_V_upred=[]
    p_deaths=bar(title="Deaths in Kenya")
        deaths_med=[];deaths_lpred=[];deaths_upred=[]
    i=0
    for scenario in scenarios
        i+=1
        @load scenario[1]  results
        plot!(p_incidence,results.total_cases.med,
                ribbon=(results.total_cases.lpred, results.total_cases.upred),
                fillalpha=fillalpha,label=scenario[2],color=colors_cp[i],xticks=xticks,xlabel="months")#,xlims=(2*30,6*))
        plot!(p_incidence_V,results.country_incidence_V_ts.med,ribbon=(results.country_incidence_V_ts.lpred,results.country_incidence_V_ts.upred),
                fillalpha=fillalpha,label=scenario[2],color=colors_cp[i],xticks=xticks,xlabel="months")

        push!(cases_med,sum(results.total_cases.med))
            push!(cases_lpred,sum(results.total_cases.lpred))
            push!(cases_upred,sum(results.total_cases.upred))
        push!(cases_V_med,sum(results.total_cases_V_by_area_and_age.med))
            push!(cases_V_lpred,sum(results.total_cases_V_by_area_and_age.lpred))
            push!(cases_V_upred,sum(results.total_cases_V_by_area_and_age.upred))
        push!(deaths_med,results.total_deaths.med);push!(deaths_lpred,results.total_deaths.lpred);push!(deaths_upred,results.total_deaths.upred);
    end
    #redoing the medians
    i=0
    for scenario in scenarios
        i+=1
        @load scenario[1]  results
        plot!(p_incidence,results.total_cases.med,
                label=false,color=colors_cp[i])#,xticks=xticks,xlabel="months")
    end

    bar!(p_cases,[scenario[2] for scenario in scenarios],cases_med,errorbar=(cases_lpred,cases_upred),
        color=colors_cp,legend=false,linecolor=:transparent,markerstrokecolor=:black,markerstrokewidth=.7,xrotation=xrotation,wide=:true)
        hline!(p_cases,[cases_med[1]],color=colors_cp[1])
    bar!(p_cases_V,[scenario[2] for scenario in scenarios],cases_V_med,errorbar=(cases_V_lpred,cases_V_upred),color=colors_cp,legend=false)
        hline!(p_cases_V,[cases_V_med[1]],color=colors_cp[1])
    bar!(p_deaths,[scenario[2] for scenario in scenarios],deaths_med,errorbar=(deaths_lpred,deaths_upred),
        color=colors_cp,legend=false,linecolor=:transparent,markerstrokecolor=:black,markerstrokewidth=.7,xrotation=xrotation,wide=:true)
        hline!(p_deaths,[deaths_med[1]],color=colors_cp[1])

    #saving plots
    if !isdir(output_folder)    mkdir(output_folder);   end
    savefig(p_incidence,output_folder*"incidence.png");savefig(plot!(p_incidence,xlims=(2*30+15,11*30),size=(320,320)),output_folder*"incidence_2.png")
    savefig(p_incidence_V,output_folder*"incidence_V.png")
    savefig(p_cases,output_folder*"cases.png");savefig(plot!(p_cases,size=(250,380)),output_folder*"cases_2.png")
    savefig(p_cases_V,output_folder*"cases_V.png")
    savefig(p_deaths,output_folder*"deaths.png");savefig(plot!(p_deaths,size=(320,380)),output_folder*"deaths_2.png")

    df=DataFrame(files=[scenario[1] for scenario in scenarios],
                labels=[scenario[2] for scenario in scenarios],
                cases_fiff=[cases_med[1]-cases_med[i]    for i=1:size(cases_med,1)],cases_med=cases_med,cases_lpred=cases_lpred,cases_upred=cases_upred,
                deaths_fiff=[deaths_med[1]-deaths_med[i]    for i=1:size(deaths_med,1)],deaths_med=deaths_med,deaths_lpred=deaths_lpred,deaths_upred=deaths_upred)
    CSV.write(output_folder*"2020_"*Dates.format(now(), "mm_dd__HH_MM") *".csv",df)
    return cases_med,cases_lpred,cases_upred,deaths_med,deaths_lpred,deaths_upred
end

######
scenarios=[("./Screening/output/0_NoInterv/session1_500sims/noInterv_sc1_500sims.jld2","No Intervention"),
         #=("./Screening/output/1_CTH/session1_500sims/CTH_sc1_500sims.jld2","CTH (CTn=30)"),
         ("./Screening/output/1_CTH/session2_500sims/CTH_sc1_500sims.jld2","CTH (CTn=100)"),=#

         #=("./Screening/output/2_SympS/session1_500sims/SympS_sc1_500sims.jld2","SympS 1e3 m3 NM"),
             ("./Screening/output/2_SympS/session1_500sims/SympS_sc2_500sims.jld2","SympS 1e3 m4 NM"),
             ("./Screening/output/2_SympS/session1_500sims/SympS_sc3_500sims.jld2","SympS 1e3 m5 NM"),
             ("./Screening/output/2_SympS/session1_500sims/SympS_sc4_500sims.jld2","SympS 1e3 m6 NM"),
             ("./Screening/output/2_SympS/session1_500sims/SympS_sc5_500sims.jld2","SympS 1e3 m7 MN"),
             ("./Screening/output/2_SympS/session1_500sims/SympS_sc6_500sims.jld2","SympS 5e3 m3 NM"),
             ("./Screening/output/2_SympS/session1_500sims/SympS_sc7_500sims.jld2","SympS 5e3 m4 NM"),
             ("./Screening/output/2_SympS/session1_500sims/SympS_sc8_500sims.jld2","SympS 5e3 m5 NM"),
             ("./Screening/output/2_SympS/session1_500sims/SympS_sc9_500sims.jld2","SympS 5e3 m6 NM"),
             ("./Screening/output/2_SympS/session1_500sims/SympS_sc10_500sims.jld2","SympS 5e3 m7 NM"),

         ("./Screening/output/2_SympS/session2_500sims/SympS_sc1_500sims.jld2","SympS 1e3 m3 K"),
             ("./Screening/output/2_SympS/session2_500sims/SympS_sc2_500sims.jld2","SympS 1e3 m4 K"),
             ("./Screening/output/2_SympS/session2_500sims/SympS_sc3_500sims.jld2","SympS 1e3 m5 K"),
             ("./Screening/output/2_SympS/session2_500sims/SympS_sc4_500sims.jld2","SympS 1e3 m6 K"),
             ("./Screening/output/2_SympS/session2_500sims/SympS_sc5_500sims.jld2","SympS 1e3 m7 K"),
             ("./Screening/output/2_SympS/session2_500sims/SympS_sc6_500sims.jld2","SympS 5e3 m3 K"),
             ("./Screening/output/2_SympS/session2_500sims/SympS_sc7_500sims.jld2","SympS 5e3 m4 K"),
             ("./Screening/output/2_SympS/session2_500sims/SympS_sc8_500sims.jld2","SympS 5e3 m5 K"),
             ("./Screening/output/2_SympS/session2_500sims/SympS_sc9_500sims.jld2","SympS 5e3 m6 K"),
             ("./Screening/output/2_SympS/session2_500sims/SympS_sc10_500sims.jld2","SympS 5e3 m7 K"),=#

         #=("./Screening/output/3_SympSCT/session1_500sims/SympSCT_sc1_500sims.jld2","SympSCT 1e3 m3 NM CTn30"),
             ("./Screening/output/3_SympSCT/session1_500sims/SympSCT_sc2_500sims.jld2","SympSCT 1e3 m4 NM CTn30"),
             ("./Screening/output/3_SympSCT/session1_500sims/SympSCT_sc3_500sims.jld2","SympSCT 1e3 m5 NM CTn30"),
             ("./Screening/output/3_SympSCT/session1_500sims/SympSCT_sc4_500sims.jld2","SympSCT 1e3 m6 NM CTn30"),
             ("./Screening/output/3_SympSCT/session1_500sims/SympSCT_sc5_500sims.jld2","SympSCT 1e3 m7 MN CTn30"),
             ("./Screening/output/3_SympSCT/session1_500sims/SympSCT_sc6_500sims.jld2","SympSCT 5e3 m3 NM CTn30"),
             ("./Screening/output/3_SympSCT/session1_500sims/SympSCT_sc7_500sims.jld2","SympSCT 5e3 m4 NM CTn30"),
             ("./Screening/output/3_SympSCT/session1_500sims/SympSCT_sc8_500sims.jld2","SympSCT 5e3 m5 NM CTn30"),
             ("./Screening/output/3_SympSCT/session1_500sims/SympSCT_sc9_500sims.jld2","SympSCT 5e3 m6 NM CTn30"),
             ("./Screening/output/3_SympSCT/session1_500sims/SympSCT_sc10_500sims.jld2","SympSCT 5e3 m7 NM CTn30"),

         ("./Screening/output/3_SympSCT/session2_500sims/SympSCT_sc1_500sims.jld2","SympSCT 1e3 m3 K CTn30"),
             ("./Screening/output/3_SympSCT/session2_500sims/SympSCT_sc2_500sims.jld2","SympSCT 1e3 m4 K CTn30"),
             ("./Screening/output/3_SympSCT/session2_500sims/SympSCT_sc3_500sims.jld2","SympSCT 1e3 m5 K CTn30"),
             ("./Screening/output/3_SympSCT/session2_500sims/SympSCT_sc4_500sims.jld2","SympSCT 1e3 m6 K CTn30"),
             ("./Screening/output/3_SympSCT/session2_500sims/SympSCT_sc5_500sims.jld2","SympSCT 1e3 m7 K CTn30"),
             ("./Screening/output/3_SympSCT/session2_500sims/SympSCT_sc6_500sims.jld2","SympSCT 5e3 m3 K CTn30"),
             ("./Screening/output/3_SympSCT/session2_500sims/SympSCT_sc7_500sims.jld2","SympSCT 5e3 m4 K CTn30"),
             ("./Screening/output/3_SympSCT/session2_500sims/SympSCT_sc8_500sims.jld2","SympSCT 5e3 m5 K CTn30"),
             ("./Screening/output/3_SympSCT/session2_500sims/SympSCT_sc9_500sims.jld2","SympSCT 5e3 m6 K CTn30"),
             ("./Screening/output/3_SympSCT/session2_500sims/SympSCT_sc10_500sims.jld2","SympSCT 5e3 m7 K CTn30"),

         ("./Screening/output/3_SympSCT/session3_500sims/SympSCT_sc1_500sims.jld2","SympSCT 1e3 m3 NM CTn30"),
             ("./Screening/output/3_SympSCT/session3_500sims/SympSCT_sc2_500sims.jld2","SympSCT 1e3 m4 NM CTn100"),
             ("./Screening/output/3_SympSCT/session3_500sims/SympSCT_sc3_500sims.jld2","SympSCT 1e3 m5 NM CTn100"),
             ("./Screening/output/3_SympSCT/session3_500sims/SympSCT_sc4_500sims.jld2","SympSCT 1e3 m6 NM CTn100"),
             ("./Screening/output/3_SympSCT/session3_500sims/SympSCT_sc5_500sims.jld2","SympSCT 1e3 m7 MN CTn100"),
             ("./Screening/output/3_SympSCT/session3_500sims/SympSCT_sc6_500sims.jld2","SympSCT 5e3 m3 NM CTn100"),
             ("./Screening/output/3_SympSCT/session3_500sims/SympSCT_sc7_500sims.jld2","SympSCT 5e3 m4 NM CTn100"),
             ("./Screening/output/3_SympSCT/session3_500sims/SympSCT_sc8_500sims.jld2","SympSCT 5e3 m5 NM CTn100"),
             ("./Screening/output/3_SympSCT/session3_500sims/SympSCT_sc9_500sims.jld2","SympSCT 5e3 m6 NM CTn100"),
             ("./Screening/output/3_SympSCT/session3_500sims/SympSCT_sc10_500sims.jld2","SympSCT 5e3 m7 NM CTn100"),

         ("./Screening/output/3_SympSCT/session4_500sims/SympSCT_sc1_500sims.jld2","SympSCT 1e3 m3 K CTn100"),
             ("./Screening/output/3_SympSCT/session4_500sims/SympSCT_sc2_500sims.jld2","SympSCT 1e3 m4 K CTn100"),
             ("./Screening/output/3_SympSCT/session4_500sims/SympSCT_sc3_500sims.jld2","SympSCT 1e3 m5 K CTn100"),
             ("./Screening/output/3_SympSCT/session4_500sims/SympSCT_sc4_500sims.jld2","SympSCT 1e3 m6 K CTn100"),
             ("./Screening/output/3_SympSCT/session4_500sims/SympSCT_sc5_500sims.jld2","SympSCT 1e3 m7 K CTn100"),
             ("./Screening/output/3_SympSCT/session4_500sims/SympSCT_sc6_500sims.jld2","SympSCT 5e3 m3 K CTn100"),
             ("./Screening/output/3_SympSCT/session4_500sims/SympSCT_sc7_500sims.jld2","SympSCT 5e3 m4 K CTn100"),
             ("./Screening/output/3_SympSCT/session4_500sims/SympSCT_sc8_500sims.jld2","SympSCT 5e3 m5 K CTn100"),
             ("./Screening/output/3_SympSCT/session4_500sims/SympSCT_sc9_500sims.jld2","SympSCT 5e3 m6 K CTn100"),
             ("./Screening/output/3_SympSCT/session4_500sims/SympSCT_sc10_500sims.jld2","SympSCT 5e3 m7 K CTn100")#,=#

         ("./Screening/output/4_MS/session1_500sims/MS_sc1_500sims.jld2","MS 1e3 m3 NM"),
             ("./Screening/output/4_MS/session1_500sims/MS_sc2_500sims.jld2","MS 1e3 m4 NM"),
             ("./Screening/output/4_MS/session1_500sims/MS_sc3_500sims.jld2","MS 1e3 m5 NM"),
             ("./Screening/output/4_MS/session1_500sims/MS_sc4_500sims.jld2","MS 1e3 m6 NM"),
             ("./Screening/output/4_MS/session1_500sims/MS_sc5_500sims.jld2","MS 1e3 m7 MN"),
             ("./Screening/output/4_MS/session1_500sims/MS_sc6_500sims.jld2","MS 5e3 m3 NM"),
             ("./Screening/output/4_MS/session1_500sims/MS_sc7_500sims.jld2","MS 5e3 m4 NM"),
             ("./Screening/output/4_MS/session1_500sims/MS_sc8_500sims.jld2","MS 5e3 m5 NM"),
             ("./Screening/output/4_MS/session1_500sims/MS_sc9_500sims.jld2","MS 5e3 m6 NM"),
             ("./Screening/output/4_MS/session1_500sims/MS_sc10_500sims.jld2","MS 5e3 m7 NM"),

         ("./Screening/output/4_MS/session2_500sims/MS_sc1_500sims.jld2","MS 1e3 m3 K"),
             ("./Screening/output/4_MS/session2_500sims/MS_sc2_500sims.jld2","MS 1e3 m4 K"),
             ("./Screening/output/4_MS/session2_500sims/MS_sc3_500sims.jld2","MS 1e3 m5 K"),
             ("./Screening/output/4_MS/session2_500sims/MS_sc4_500sims.jld2","MS 1e3 m6 K"),
             ("./Screening/output/4_MS/session2_500sims/MS_sc5_500sims.jld2","MS 1e3 m7 K"),
             ("./Screening/output/4_MS/session2_500sims/MS_sc6_500sims.jld2","MS 5e3 m3 K"),
             ("./Screening/output/4_MS/session2_500sims/MS_sc7_500sims.jld2","MS 5e3 m4 K"),
             ("./Screening/output/4_MS/session2_500sims/MS_sc8_500sims.jld2","MS 5e3 m5 K"),
             ("./Screening/output/4_MS/session2_500sims/MS_sc9_500sims.jld2","MS 5e3 m6 K"),
             #("./Screening/output/4_MS/session2_500sims/MS_sc10_500sims.jld2","MS 5e3 m7 K")#,

         #=("./Screening/output/5_MSCT/session1_500sims/MSCT_sc1_500sims.jld2","MSCT 1e3 m3 NM CTn30"),
             ("./Screening/output/5_MSCT/session1_500sims/MSCT_sc2_500sims.jld2","MSCT 1e3 m4 NM CTn30"),
             ("./Screening/output/5_MSCT/session1_500sims/MSCT_sc3_500sims.jld2","MSCT 1e3 m5 NM CTn30"),
             ("./Screening/output/5_MSCT/session1_500sims/MSCT_sc4_500sims.jld2","MSCT 1e3 m6 NM CTn30"),
             ("./Screening/output/5_MSCT/session1_500sims/MSCT_sc5_500sims.jld2","MSCT 1e3 m7 MN CTn30"),
             ("./Screening/output/5_MSCT/session1_500sims/MSCT_sc6_500sims.jld2","MSCT 5e3 m3 NM CTn30"),
             ("./Screening/output/5_MSCT/session1_500sims/MSCT_sc7_500sims.jld2","MSCT 5e3 m4 NM CTn30"),
             ("./Screening/output/5_MSCT/session1_500sims/MSCT_sc8_500sims.jld2","MSCT 5e3 m5 NM CTn30"),
             ("./Screening/output/5_MSCT/session1_500sims/MSCT_sc9_500sims.jld2","MSCT 5e3 m6 NM CTn30"),
             ("./Screening/output/5_MSCT/session1_500sims/MSCT_sc10_500sims.jld2","MSCT 5e3 m7 NM CTn30"),

         ("./Screening/output/5_MSCT/session2_500sims/MSCT_sc1_500sims.jld2","MSCT 1e3 m3 K CTn30"),
             ("./Screening/output/5_MSCT/session2_500sims/MSCT_sc2_500sims.jld2","MSCT 1e3 m4 K CTn30"),
             ("./Screening/output/5_MSCT/session2_500sims/MSCT_sc3_500sims.jld2","MSCT 1e3 m5 K CTn30"),
             ("./Screening/output/5_MSCT/session2_500sims/MSCT_sc4_500sims.jld2","MSCT 1e3 m6 K CTn30"),
             ("./Screening/output/5_MSCT/session2_500sims/MSCT_sc5_500sims.jld2","MSCT 1e3 m7 K CTn30"),
             ("./Screening/output/5_MSCT/session2_500sims/MSCT_sc6_500sims.jld2","MSCT 5e3 m3 K CTn30"),
             ("./Screening/output/5_MSCT/session2_500sims/MSCT_sc7_500sims.jld2","MSCT 5e3 m4 K CTn30"),
             ("./Screening/output/5_MSCT/session2_500sims/MSCT_sc8_500sims.jld2","MSCT 5e3 m5 K CTn30"),
             ("./Screening/output/5_MSCT/session2_500sims/MSCT_sc9_500sims.jld2","MSCT 5e3 m6 K CTn30"),
             ("./Screening/output/5_MSCT/session2_500sims/MSCT_sc10_500sims.jld2","MSCT 5e3 m7 K CTn30"),

         ("./Screening/output/5_MSCT/session3_500sims/MSCT_sc1_500sims.jld2","MSCT 1e3 m3 NM CTn100"),
             ("./Screening/output/5_MSCT/session3_500sims/MSCT_sc2_500sims.jld2","MSCT 1e3 m4 NM CTn100"),
             ("./Screening/output/5_MSCT/session3_500sims/MSCT_sc3_500sims.jld2","MSCT 1e3 m5 NM CTn100"),
             ("./Screening/output/5_MSCT/session3_500sims/MSCT_sc4_500sims.jld2","MSCT 1e3 m6 NM CTn100"),
             ("./Screening/output/5_MSCT/session3_500sims/MSCT_sc5_500sims.jld2","MSCT 1e3 m7 MN CTn100"),
             ("./Screening/output/5_MSCT/session3_500sims/MSCT_sc6_500sims.jld2","MSCT 5e3 m3 NM CTn100"),
             ("./Screening/output/5_MSCT/session3_500sims/MSCT_sc7_500sims.jld2","MSCT 5e3 m4 NM CTn100"),
             ("./Screening/output/5_MSCT/session3_500sims/MSCT_sc8_500sims.jld2","MSCT 5e3 m5 NM CTn100"),
             ("./Screening/output/5_MSCT/session3_500sims/MSCT_sc9_500sims.jld2","MSCT 5e3 m6 NM CTn100"),
             ("./Screening/output/5_MSCT/session3_500sims/MSCT_sc10_500sims.jld2","MSCT 5e3 m7 NM CTn100"),

         ("./Screening/output/5_MSCT/session4_500sims/MSCT_sc1_500sims.jld2","MSCT 1e3 m3 K CTn100"),
             ("./Screening/output/5_MSCT/session4_500sims/MSCT_sc2_500sims.jld2","MSCT 1e3 m4 K CTn100"),
             ("./Screening/output/5_MSCT/session4_500sims/MSCT_sc3_500sims.jld2","MSCT 1e3 m5 K CTn100"),
             ("./Screening/output/5_MSCT/session4_500sims/MSCT_sc4_500sims.jld2","MSCT 1e3 m6 K CTn100"),
             ("./Screening/output/5_MSCT/session4_500sims/MSCT_sc5_500sims.jld2","MSCT 1e3 m7 K CTn100"),
             ("./Screening/output/5_MSCT/session4_500sims/MSCT_sc6_500sims.jld2","MSCT 5e3 m3 K CTn100"),
             ("./Screening/output/5_MSCT/session4_500sims/MSCT_sc7_500sims.jld2","MSCT 5e3 m4 K CTn100"),
             ("./Screening/output/5_MSCT/session4_500sims/MSCT_sc8_500sims.jld2","MSCT 5e3 m5 K CTn100"),
             ("./Screening/output/5_MSCT/session4_500sims/MSCT_sc9_500sims.jld2","MSCT 5e3 m6 K CTn100"),
             ("./Screening/output/5_MSCT/session4_500sims/MSCT_sc10_500sims.jld2","MSCT 5e3 m7 K CTn100")=#

         ]

#=scenarios=[("./Screening/output/0_NoInterv/session1_500sims/noInterv_sc1_500sims.jld2","No Intervention"),
            ("./Screening/output/1_CTH/session1_500sims/CTH_sc1_500sims.jld2","CTH CTn=30"),
            ("./Screening/output/1_CTH/session2_500sims/CTH_sc1_500sims.jld2","CTH CTn=100")]=#
#=scenarios=[("./Screening/output/0_NoInterv/session1_500sims/noInterv_sc1_500sims.jld2","No Intervention"),#BEST OF MSCT
            ("./Screening/output/5_MSCT/session1_500sims/MSCT_sc5_500sims.jld2","MS 1e3 NM CTn30"),
            ("./Screening/output/5_MSCT/session1_500sims/MSCT_sc10_500sims.jld2","MS 5e3 NM CTn30"),
            ("./Screening/output/5_MSCT/session3_500sims/MSCT_sc5_500sims.jld2","MS 1e3 NM CTn100"),
            ("./Screening/output/5_MSCT/session3_500sims/MSCT_sc9_500sims.jld2","MS 5e3 NM CTn100"),
            ("./Screening/output/5_MSCT/session2_500sims/MSCT_sc4_500sims.jld2","MS 1e3 K CTn30"),
            ("./Screening/output/5_MSCT/session2_500sims/MSCT_sc9_500sims.jld2","MS 5e3 K CTn30"),
            ("./Screening/output/5_MSCT/session4_500sims/MSCT_sc4_500sims.jld2","MS 1e3 K CTn100"),
            ("./Screening/output/5_MSCT/session4_500sims/MSCT_sc9_500sims.jld2","MS 5e3 K CTn100")]=#
#=scenarios=[("./Screening/output/3_SympSCT/session1_500sims/SympSCT_sc4_500sims.jld2","SympS 1e3 NM CTn30"),#BEST OF SympSCT
            ("./Screening/output/3_SympSCT/session1_500sims/SympSCT_sc8_500sims.jld2","SympS 5e3 NM CTn30"),
            ("./Screening/output/3_SympSCT/session3_500sims/SympSCT_sc4_500sims.jld2","SympS 1e3 NM CTn100"),
            ("./Screening/output/3_SympSCT/session3_500sims/SympSCT_sc6_500sims.jld2","SympS 5e3 NM CTn100"),
            ("./Screening/output/3_SympSCT/session2_500sims/SympSCT_sc3_500sims.jld2","SympS 1e3 K CTn30"),
            ("./Screening/output/3_SympSCT/session2_500sims/SympSCT_sc9_500sims.jld2","SympS 5e3 K CTn30"),
            ("./Screening/output/3_SympSCT/session4_500sims/SympSCT_sc4_500sims.jld2","SympS 1e3 K CTn100"),
            ("./Screening/output/3_SympSCT/session4_500sims/SympSCT_sc9_500sims.jld2","SympS 5e3 K CTn100")]=#
scenarios=[("./Screening/output/0_NoInterv/session1_500sims/noInterv_sc1_500sims.jld2","No Intervention"),
            ("./Screening/output/6_SD/session1_500sims/SD_sc1_500sims.jld2","SD 20%"),
            ("./Screening/output/6_SD/session2_500sims/SD_sc1_500sims.jld2","SD 40%"),
            ("./Screening/output/7_CTHSD/session1_500sims/CTHSD_sc1_500sims.jld2","CTH CTn=30 SD20%"),
            ("./Screening/output/7_CTHSD/session2_500sims/CTHSD_sc1_500sims.jld2","CTH CTn=30 SD40%"),
            ("./Screening/output/7_CTHSD/session3_500sims/CTHSD_sc1_500sims.jld2","CTH CTn=100 SD20%"),
            ("./Screening/output/7_CTHSD/session4_500sims/CTHSD_sc1_500sims.jld2","CTH CTn=100 SD40%")]
scenarios=[("./Screening/output/0_NoInterv/session1_500sims/noInterv_sc1_500sims.jld2","No Intervention"),
            ("./Screening/output/1_CTH/session1_500sims/CTH_sc1_500sims.jld2","CTH CTn=30"),
            ("./Screening/output/1_CTH/session2_500sims/CTH_sc1_500sims.jld2","CTH CTn=100"),
            ("./Screening/output/6_SD/session1_500sims/SD_sc1_500sims.jld2","SD 20%"),
            ("./Screening/output/7_CTHSD/session1_500sims/CTHSD_sc1_500sims.jld2","CTH CTn=30 SD20%"),
            ("./Screening/output/7_CTHSD/session3_500sims/CTHSD_sc1_500sims.jld2","CTH CTn=100 SD20%"),
            ("./Screening/output/6_SD/session2_500sims/SD_sc1_500sims.jld2","SD 40%"),
            ("./Screening/output/7_CTHSD/session2_500sims/CTHSD_sc1_500sims.jld2","CTH CTn=30 SD40%"),
            ("./Screening/output/7_CTHSD/session4_500sims/CTHSD_sc1_500sims.jld2","CTH CTn=100 SD40%")]
         rez=make_plots(scenarios,"Screening/plotting/plots__2020_"*Dates.format(now(), "mm_dd__HH_MM")*"__CTH_SD_all/")

## Social distancing
#=scenarios=[#("./Screening/output/0_NoInterv/session1_500sims/noInterv_sc1_500sims.jld2","No Intervention"),
            ("./Screening/output/6_SD/session1_500sims/SD_sc1_500sims.jld2","SD 20%"),
            ("./Screening/output/6_SD/session2_500sims/SD_sc1_500sims.jld2","SD 40%")]
            rez=make_plots_withTOTALCASES(scenarios,"Screening/plotting/plots__2020_"*Dates.format(now(), "mm_dd__HH_MM")*"__SD/")=#



#######
[rez[1][1]-rez[1][i]    for i=2:size(rez[1],1)]#diff in cases_med

[rez[4][1]-rez[4][i]    for i=2:size(rez[4],1)]#diff in deaths_med

#@load "./Screening/output/3_SympSCT/session1_2sims/SympSCT_sc7_2sims.jld2" results
#    sum(results.total_q_by_area_and_age.med)

##get main numbers from a specific scenario
@load "./Screening/output/0_NoInterv/session1_500sims/noInterv_sc1_500sims.jld2" results
@load "./Screening/output/6_SD/session1_2sims/SD_sc1_2sims.jld2" results
    R1=deepcopy(results)
@load "./Screening/output/1_CTH/session2_500sims/CTH_sc1_500sims.jld2" results
@load "./Screening/output/7_CTHSD/session1_2sims/CTHSD_sc1_2sims.jld2" results
    R2=deepcopy(results)

#CumI:
"CumI=",sum(R1.total_cases_A_by_area_and_age.med)+sum(R1.total_cases_M_by_area_and_age.med)+sum(R1.total_cases_V_by_area_and_age.med),
            sum(R1.total_cases_A_by_area_and_age.lpred)+sum(R1.total_cases_M_by_area_and_age.lpred)+sum(R1.total_cases_V_by_area_and_age.lpred),
            sum(R1.total_cases_A_by_area_and_age.upred)+sum(R1.total_cases_M_by_area_and_age.upred)+sum(R1.total_cases_V_by_area_and_age.upred)
"cumV=",sum(R1.total_cases_V_by_area_and_age.med),
            sum(R1.total_cases_V_by_area_and_age.lpred),
            sum(R1.total_cases_V_by_area_and_age.upred)
"Deaths=",R1.total_deaths.med,R1.total_deaths.lpred,R1.total_deaths.upred

"cumQ=",sum(R2.total_q_by_area_and_age.med),
            sum(R2.total_q_by_area_and_age.lpred),
            sum(R2.total_q_by_area_and_age.upred)

"cumI difference",sum(R1.total_cases_A_by_area_and_age.med)+sum(R1.total_cases_M_by_area_and_age.med)+sum(R1.total_cases_V_by_area_and_age.med) - sum(R2.total_cases_A_by_area_and_age.med)-sum(R2.total_cases_M_by_area_and_age.med)-sum(R2.total_cases_V_by_area_and_age.med)

"Deaths difference=",R1.total_deaths.med - R2.total_deaths.med
