push!(LOAD_PATH, "./src")
    push!(LOAD_PATH, "./Screening")
    using JLD2,Plots,Colors,ColorBrewer,ColorSchemes,Dates,CSV,DataFrames

#######
Plots.default(grid=false,legendfontsize=6,titlefontsize=8, xtickfontsize=8, ytickfontsize=8,fontfamily="Cambria",labelfontsize=9)#,formatter=:scientific)
    fillalpha=.2
    county_index_N=30;county_index_M=28
    colors=[:red;ColorBrewer.palette("Set2",7)[1:5];ColorBrewer.palette("Set3",12)]
    colors=[ColorBrewer.palette("Set2",8);ColorBrewer.palette("Accent",8)]
    #color_palette=:seaborn_bright
    colors=[ColorSchemes.seaborn_bright.colors[i]   for i=1:size(ColorSchemes.seaborn_bright.colors,1)]
        c=colors[1];colors[1]=colors[4];colors[4]=c;
    size_x=250;size_y=350;xrotation=50;
    S0= 4.7562088e7;

###
function make_plots_bestOfs(scenarios,output_folder,intervention_label)
    @load scenarios[end][1]  results
    sim_dur=size(results.country_incidence_I_ts.med,1)
    xticks=([i  for i=1:sim_dur    if i%30==0],[string(Int(i/30))  for i=1:sim_dur    if i%30==0])
    xlabels=["No interv"]
    for scenario in scenarios[2:end]
        if occursin("hospitalized",intervention_label)
            push!(xlabels, #=string(scenario[4])*" tests for "*=#string(scenario[2])*" m");
        else
            push!(xlabels, string(scenario[4])*" tests for "*string(scenario[2])*" m");
        end
    end
    colors_cp=deepcopy(colors)
    if size(scenarios,1)>size(colors_cp,1)
        colors_cp=repeat(colors_cp,Int(ceil(size(scenarios,1)/size(colors_cp,1))))
    end
    p_incidence=plot(title=intervention_label*"\nIncidence in Kenya")
    p_cases=plot(title="Total infecteds in Kenya-"*intervention_label)
        cases_med=[];cases_lpred=[];cases_upred=[]
    p_MV_cases=plot(title="Total symptomatic and severe in Kenya-"*intervention_label)
        cases_MV_med=[];cases_MV_lpred=[];cases_MV_upred=[]
    p_deaths=plot(title="Deaths in Kenya-"*intervention_label)
        deaths_med=[];deaths_lpred=[];deaths_upred=[]
    p_q=plot(title="Quarantined infecteds in Kenya-"*intervention_label)
        q_med=[];q_lpred=[];q_upred=[]
    p_qs=plot(title="Quarantined S in Kenya-"*intervention_label)
        qs_med=[];qs_lpred=[];qs_upred=[]
    #starting with the vspans of the incidence plot
    i=1
    for scenario in scenarios[2:end]
        i+=1
        if occursin("hospitalized",intervention_label)
            plot!(p_incidence,[(scenario[3]-1)*30, (scenario[3]-1)*30+30*scenario[2]],
                [(i-1)*10^floor(log10(maximum(results.country_incidence_I_ts.med)))/2,(i-1)*10^floor(log10(maximum(results.country_incidence_I_ts.med)))/2],
                color=colors_cp[i],label=string(scenario[2])*" m duration",markershape=:vline,ls=(scenario[4]==500 ? :dash : :solid))
        else
            plot!(p_incidence,[(scenario[3]-1)*30, (scenario[3]-1)*30+30*scenario[2]],
                [(i-1)*10^floor(log10(maximum(results.country_incidence_I_ts.med)))/2,(i-1)*10^floor(log10(maximum(results.country_incidence_I_ts.med)))/2],
                color=colors_cp[i],label=string(scenario[4])*" tests",markershape=:vline,ls=(scenario[4]==500 ? :dash : :solid))
        end
    end
    i=0
    for scenario in scenarios
        i+=1
        @load scenario[1]  results;
        L=i!=1 ? string(scenario[2])*"m" : "No interv"
        plot!(p_incidence,results.country_incidence_I_ts.med,
                ribbon=(results.country_incidence_I_ts.lpred,
                results.country_incidence_I_ts.upred),ls=(scenario[4]==500 ? :dash : :solid),
                fillalpha=fillalpha,label=L,color=colors_cp[i],xticks=xticks,xlabel="months",
                yformatter = yi -> yi==0 ? "0" : string(#=Int=#(yi/10^(floor(log10(yi))))),ylabel="x10^"*string(Int(floor(log10(maximum(results.country_incidence_I_ts.med))))))

        push!(cases_med,results.total_cases.med)#sum(results.total_cases_A_by_area_and_age.med)+sum(results.total_cases_M_by_area_and_age.med)+sum(results.total_cases_V_by_area_and_age.med))
            push!(cases_lpred,results.total_cases.lpred)#sum(results.total_cases_A_by_area_and_age.lpred)+sum(results.total_cases_M_by_area_and_age.lpred)+sum(results.total_cases_V_by_area_and_age.lpred))
            push!(cases_upred,results.total_cases.upred)#sum(results.total_cases_A_by_area_and_age.upred)+sum(results.total_cases_M_by_area_and_age.upred)+sum(results.total_cases_V_by_area_and_age.upred))
        push!(cases_MV_med,results.total_cases_M.med+results.total_severe_cases.med) #sum(results.total_cases_M_by_area_and_age.med)+sum(results.total_cases_V_by_area_and_age.med))
            push!(cases_MV_lpred,results.total_cases_M.lpred+results.total_severe_cases.lpred) #sum(results.total_cases_M_by_area_and_age.lpred)+sum(results.total_cases_V_by_area_and_age.lpred))
            push!(cases_MV_upred,results.total_cases_M.upred+results.total_severe_cases.upred) #sum(results.total_cases_M_by_area_and_age.upred)+sum(results.total_cases_V_by_area_and_age.upred))
        push!(deaths_med,results.total_deaths.med);push!(deaths_lpred,results.total_deaths.lpred);push!(deaths_upred,results.total_deaths.upred);
        push!(q_med,results.total_q.med)
            push!(q_lpred,results.total_q.lpred)
            push!(q_upred,results.total_q.upred)
        push!(qs_med,results.total_qs.med)
            push!(qs_lpred,results.total_qs.lpred)
            push!(qs_upred,results.total_qs.upred)
    end
    #redoing the medians
    i=0
    for scenario in scenarios
        i+=1
        @load scenario[1]  results
        plot!(p_incidence,results.country_incidence_I_ts.med,label=false,color=colors_cp[i],ls=(scenario[4]==500 ? :dash : :solid))#,xticks=xticks,xlabel="months")
    end

    scatter!(p_cases,xlabels,cases_med,errorbar=(abs.(cases_med.-cases_lpred),abs.(cases_med.-cases_upred)),markershape=:square,
        yformatter = yi -> yi==0 ? "0" : string(#=Int=#(yi/10^(floor(log10(yi))))),ylabel="x10^"*string(Int(floor(log10(maximum(cases_med))))),
        color=colors_cp,legend=false,linecolor=:transparent,markerstrokecolor=:black,markerstrokewidth=.7,xrotation=xrotation,wide=:true)
        hline!(p_cases,[cases_med[1]],color=colors_cp[1])

    scatter!(p_MV_cases,xlabels,cases_MV_med,errorbar=(abs.(cases_MV_med.-cases_MV_lpred),abs.(cases_MV_med.-cases_MV_upred)),markershape=:square,
        yformatter = yi -> yi==0 ? "0" : string(#=Int=#(yi/10^(floor(log10(yi))))),ylabel="x10^"*string(Int(floor(log10(maximum(cases_MV_med))))),
        color=colors_cp,legend=false,linecolor=:transparent,markerstrokecolor=:black,markerstrokewidth=.7,xrotation=xrotation,wide=:true)
        hline!(p_MV_cases,[cases_MV_med[1]],color=colors_cp[1])
    scatter!(p_deaths,xlabels,deaths_med,errorbar=(abs.(deaths_med.-deaths_lpred),abs.(deaths_med.-deaths_upred)),markershape=:square,
        yformatter = yi -> yi==0 ? "0" : string(yi/10^(floor(log10(yi)))),ylabel="x10^"*string(Int(floor(log10(maximum(deaths_med))))),
        color=colors_cp,legend=false,linecolor=:transparent,markerstrokecolor=:black,markerstrokewidth=.7,xrotation=xrotation,wide=:true)
        hline!(p_deaths,[deaths_med[1]],color=colors_cp[1])
    if maximum(qs_med)!=0
        scatter!(p_q,xlabels[1:end],q_med[1:end],errorbar=(abs.(q_med.-q_lpred),abs.(q_med.-q_upred)),markershape=:square,
            yformatter = yi -> yi==0 ? "0" : string(yi/10^(floor(log10(yi)))),ylabel="x10^"*string(Int(floor(log10(maximum(q_med))))),
            color=colors_cp[1:end],legend=false,linecolor=:transparent,markerstrokecolor=:black,markerstrokewidth=.7,xrotation=xrotation,wide=:true)
            #hline!(p_q,[q_med[1]],color=colors_cp[1])
        scatter!(p_qs,xlabels[1:end],qs_med[1:end],errorbar=(abs.(qs_med.-qs_lpred),abs.(qs_med.-qs_upred)),markershape=:square,
            yformatter = yi -> yi==0 ? "0" : string(yi/10^(floor(log10(yi)))),ylabel="x10^"*string(Int(floor(log10(maximum(qs_med))))),
            color=colors_cp[2:end],legend=false,linecolor=:transparent,markerstrokecolor=:black,markerstrokewidth=.7,xrotation=xrotation,wide=:true)
    end

    #saving plots
    if !isdir(output_folder)    mkdir(output_folder);   end
    savefig(plot!(p_incidence,xlims=(2*30,15*30),size=(380,320)),output_folder*"incidence_2.png")
    #savefig(plot!(p_cases,ylims=(2*10^Int(floor(log10(maximum(cases_med)))),Inf),size=(220,320)),output_folder*"cases_2.png")
    #savefig(plot!(p_MV_cases,ylims=(2*10^Int(floor(log10(maximum(cases_MV_med)))),Inf),size=(220,320)),output_folder*"cases_MV_2.png")
    #savefig(plot!(p_deaths,ylims=(2*10^Int(floor(log10(maximum(deaths_med)))),Inf),size=(220,320)),output_folder*"deaths_2.png")
    #savefig(plot!(p_q,size=(200,320)),output_folder*"Q.png")
    #savefig(plot!(p_qs,size=(200,320)),output_folder*"Qs.png")
    #=savefig(plot(layout=(1,2),size=(450,230),
                plot!(p_cases,#=ylims=(2*10^Int(floor(log10(maximum(cases_med)))),Inf),=#title="Total infecteds",xticks=false),
                plot!(p_MV_cases,#=ylims=(2*10^Int(floor(log10(maximum(cases_MV_med)))),Inf),=#title="Total symptomatic and severe",xticks=false)),
                output_folder*"cases_AMV_MV.png")
    savefig(plot(layout=(1,2),size=(450,300),
                plot!(p_deaths,#=ylims=(2*10^Int(floor(log10(maximum(deaths_med)))),Inf),=#title="Deaths"),
                plot!(p_q,title="Quarantined EPAMV")),
                output_folder*"cases_D_QI.png")=#

    savefig(plot(layout=(1,4),size=(600,300),
                plot!(p_cases,title="Total infecteds"),
                plot!(p_MV_cases,title="Symptomatic and severe"),
                plot!(p_deaths,title="Deaths"),
                plot!(p_q,title="Quarantined infecteds")),
                output_folder*"cases_4plots.png")

    df=DataFrame(files=[scenario[1] for scenario in scenarios],
                labels=xlabels,
                cases_diff=[cases_med[1]-cases_med[i]    for i=1:size(cases_med,1)],cases_med=cases_med,cases_lpred=cases_lpred,cases_upred=cases_upred,
                deaths_diff=[deaths_med[1]-deaths_med[i]    for i=1:size(deaths_med,1)],deaths_med=deaths_med,deaths_lpred=deaths_lpred,deaths_upred=deaths_upred,
                q=[q_med[i]    for i=1:size(q_med,1)],q_med=q_med,q_lpred=q_lpred,q_upred=q_upred,
                qs=[qs_med[i]    for i=1:size(qs_med,1)],qs_med=qs_med,qs_lpred=qs_lpred,qs_upred=qs_upred)
    CSV.write(output_folder*"2020_"*Dates.format(now(), "mm_dd__HH_MM") *".csv",df)
end
