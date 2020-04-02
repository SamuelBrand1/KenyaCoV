using Plots,MAT, StatsPlots, Statistics,JLD2, Images, ImageView,ImageDraw,CSV,DataFrames,Printf,ColorBrewer
Plots.default(grid=false,legendfontsize=6, legendfontcolor=:gray30,titlefontsize=8, titlefontcolor=:gray30,xtickfontsize=6, ytickfontsize=6,tickfontcolor=:gray30,titlefontfamily="Cambria")
StatsPlots.default(grid=false,legendfontsize=6, legendfontcolor=:gray30,#:midnightblue#legendfont="Cambria"
                    titlefontsize=8, titlefontcolor=:gray30,#titlefontfamily="Cambria",
                    xtickfontsize=6, ytickfontsize=6,tickfontcolor=:gray30,linewidth=false)#tickfontfamily="Cambria",
colors=[:orange,:purple3,:maroon,:gold,:orangered,:grey,:purple,:ivory3,:chocolate1,:tan1,:rosybrown,:rosybrown2,:brown2,:brown3,:brown4,:deeppink3,:deeppink4]
colors2=ColorBrewer.palette("Pastel1",8);colors2=repeat(colors2,outer=[10])
markershapes=[:cross,:star4, :vcross, :star6, :hline]
riskregionnames = ["Machakos/Muranga" "Mandera" "Baringo/Nakuru" "Nairobi" "Uasin Gishu" "Marsabit" "Garissa" "Nakuru/Narok" "Turkana" "Taita Taveta" "Kitui/Meru" "Kilifi/Mombasa" "Kericho/Kisumu" "Kilifi/Lamu" "Kakamega/Kisumu" "Wajir" "Kajiado/Kisumu" "Homa bay/Migori" "Samburu/Laikipia" "Kilifi/Kwale" "Total"]
wa_coords=[[300,450], [515,85], [165,360], [235,465], [115,300], [300,140], [510,380], [180, 430], [120, 130], [355, 615], [340, 375], [465, 630], [100, 380], [495, 530], [40, 355], [490, 255], [155, 495], [30, 440], [250, 290], [400, 670]]
ages=[string(i)*"-"*string(i+4) for i=0:5:70];push!(ages,"75+")#"0-4" "5-9" "10-14" "15-19" "20-24" "25-29" "30-34" "35-39" "40-44" "45-49" ]

########
function save_sessionParams(folder,cumIs,cumI_diffs)
    data_files=readdir(folder)
    data_files=[folder*s for s in data_files if endswith(s,".jld2")&&s!="stats.jld2"]
    params=[]#("sc_nb","n_traj","taup","stop_Q","Κ_max_capacity12","k_per_event4","IDs_cfirst","dt","ext_inf_rate","E","delta","gamma","sigma","beta","tau","k","km","Dt","Κ_max_capacity","Κ_max_capacity4")]
    for i=1:size(data_files,1)
        @load data_files[i]  sessionParams
        sessionParams.cumI=cumIs[i];    sessionParams.cumI_diff=cumI_diffs[i];
        params=[params;sessionParams]
    end
    sort!(params, by = x -> x.sc_nb)
    #params=[params cumIs]
    CSV.write(folder*"sessionParams.csv", params)
    return params
end
function mysort(data_files)
    A=[]#["",0] for i=1:size(data_files,1)]
    for i=1:size(data_files,1)
        @load data_files[i]  sessionParams
        #push!(A,[data_files[i],sessionParams.sc_nb])
        A=[A; [[data_files[i] sessionParams.sc_nb]]]
    end
    sort!(A, by = x -> x[2]);
    return [A[i][1] for i=1:size(A,1)]
end
#=function make_plots_oneδ(folder,S0)
    data_files=readdir(folder);
    data_files=[folder*s for s in data_files if endswith(s,".jld2")&&s!="stats.jld2"]
    data_files=mysort(data_files)
    cumI_diffs=[];peakss=[];peak_diffss=[];    peaks=[];    peak_diffs=[];    median_peak0=0;   Cpeaks=[]
    wa=12;wa_name="Kilifi"

    @load data_files[1]  sessionParams
    b=boxplot(title="Total number of cases in "*wa_name*" - R0="*string(sessionParams.R₀)*", d="*string(sessionParams.δ))
    b4=boxplot(title="Gain in number of cases in "*wa_name)
    b7=boxplot(title="Peak delays in "*wa_name)
    final_cum0detection=0;        cumI_diffs=[];        cumIs=[];    medians_cumI=[]
    xticks_labels=["p,k_max"];
    for i=1:size(data_files,1)
        @load data_files[i]  sims_vector
        @load data_files[i]  sessionParams
        push!(xticks_labels,string(Int(sessionParams.τₚ[wa]*100))*"%,")#*@sprintf("%.0E", sessionParams.Κ_max_capacity12))
        ##final_cumI_wa:#peaks code
        final_cum=[];           peak=[];    Cpeak=[]
        for sim=1:size(sims_vector,1) #for each sim
            push!(final_cum,sum(sims_vector[sim][9][wa,:,1:3])) # final cumulative I=IA+ID+IQ, summed for all ages

            I=[sims_vector[sim][4][t][wa][1] for t=1:size(sims_vector[sim][4],1)]
            tpeak=findall(x->x==maximum(I),I)[1]
            push!(peak,tpeak)

            C=[sims_vector[sim][2][t][wa][1] for t=1:size(sims_vector[sim][2],1)]
            push!(Cpeak,findall(x->x==maximum(C),C)[1])
        end
        push!(cumIs,median(final_cum))
        boxplot!(b,final_cum,color=repeat([colors2[i]],5),legend=false)

        if i==1     final_cum0detection=median(final_cum)       end
        push!(cumI_diffs,final_cum0detection-median(final_cum))
        boxplot!(b4,final_cum0detection .- final_cum,color=repeat([colors2[i]],5),legend=false)

        ##final_cumI_barplot:
        cumI=[[sum(sims_vector[sim][1][end][wa,1:3]) for wa=1:20]   for sim=1:size(sims_vector,1)]
        median_cumI=[median([cumI[sim][wa] for sim=1:size(sims_vector,1)])      for wa=1:20]
        push!(medians_cumI,median_cumI)

        push!(peaks,median(peak))
        if i==1     median_peak0=median(peaks)       end
        push!(peak_diffs,median(peaks)-median_peak0)
        boxplot!(b7,peaks .- median_peak0,legend=false,color=repeat([colors2[i]],5))

        push!(Cpeaks,median(Cpeak))
    end
    ##final_cumI_wa:
    #println("cumI_diff=",cumI_diffs,"   cumIs=",cumIs)
    @save folder*"stats.jld2"   cumI_diffs,cumIs,peak_diffs,peaks
    CSV.write(folder*"stats.csv", DataFrame(xticks_labels=xticks_labels[2:end],cumI_diffs=cumI_diffs, cumIs=cumIs,peak_diffs=peak_diffs,peaks=peaks,Cpeaks=Cpeaks))

    b8=scatter(Cpeaks,color=colors2,title="CPeak times in "*wa_name,
            legend=false,xticks=([0:1:size(xticks_labels,1);], xticks_labels),size=(width2,length1),xrotation=Rotation)#,markershape)

    boxplot!(b, xticks=([0:1:size(xticks_labels,1);], xticks_labels),size=(width1,length1),xrotation=Rotation)#,ylims=(1.04e6,Inf))
    boxplot!(b4,xticks=([0:1:size(xticks_labels,1);], xticks_labels),size=(width2,length1),xrotation=Rotation)
    boxplot!(b7,xticks=([0:1:size(xticks_labels,1);], xticks_labels),size=(width1,length1),xrotation=Rotation)#,color=colors2)
    b2=bar(cumI_diffs,color=colors2,title="Gain in number of cases in "*wa_name*" - R0="*string(sessionParams.R₀)*", d="*string(sessionParams.δ),linewidth=false,
            legend=false,xticks=([0:1:size(xticks_labels,1);], xticks_labels),size=(width1,length1),xrotation=Rotation)
    b3=bar(cumIs,color=colors2,title="Total number of cases in "*wa_name,
            legend=false,xticks=([0:1:size(xticks_labels,1);], xticks_labels),size=(width2,length1),xrotation=Rotation)
    b5=bar(peak_diffs,color=colors2,title="Peak time difference in "*wa_name,
            legend=false,xticks=([0:1:size(xticks_labels,1);], xticks_labels),size=(width2,length1),xrotation=Rotation)
    b6=bar(peaks,color=colors2,title="Peak times in "*wa_name,
            legend=false,xticks=([0:1:size(xticks_labels,1);], xticks_labels),size=(width2,length1),xrotation=Rotation)
    savefig(b,folder*"jl_cumI_"*wa_name*"_box.png");        savefig(b2,folder*"jl_cumI_gain_"*wa_name*"_bar.png")
    savefig(b3,folder*"jl_cumI_"*wa_name*"_bar.png");       savefig(b4,folder*"jl_cumI_gain_"*wa_name*"_box.png")
    savefig(b5,folder*"jl_peak_diffs_"*wa_name*"_bar.png"); savefig(b6,folder*"jl_peaks_"*wa_name*"_bar.png")
    savefig(b7,folder*"jl_peaks_"*wa_name*"_box.png");      savefig(b8,folder*"jl_cpeaks_"*wa_name*"_bar.png")

    savefig(plot(b,b2,layout=(2,1),size=(width1,length1*2)),folder*"jl_2plots_"*wa_name*".png")
    savefig(plot(b,b2,layout=(1,2),size=(width1*2,length1)),folder*"jl_2plots_"*wa_name*"_2.png")

    ##final_cumI_barplot:
    S02=[]; for wa=1:20,i=1:size(data_files,1) push!(S02,S0[wa]);   end
    barplotALL=bar(S02,color=:green,xticks=([1:size(data_files,1):20*size(data_files,1);],riskregionnames),size=(width2,length1),xrotation=Rotation,legend=false)

    medians_cumI2=[];
    for wa=1:20,i=1:size(data_files,1)    push!(medians_cumI2,medians_cumI[i][wa]);   end
    bar!(barplotALL,medians_cumI2,color=colors2)
    savefig(barplotALL,folder*"jl_cumI_barplot.png")
end=#
function make_plots_twoδs(folders,plot_folder,S0,δs)
    if !isdir(plot_folder)   mkdir(plot_folder)   end
    wa=12;wa_name="Kilifi"
    data_files=[]
    for folder in folders
        data_files2=readdir(folder);
        data_files2=[folder*s for s in data_files2 if endswith(s,".jld2")&&s!="stats.jld2"]
        data_files=[data_files ;data_files2]
    end
    global data_files=mysort(data_files)
    cumI_diffs=[];peakss=[];peak_diffss=[];    peaks=[];    peak_diffs=[];    median_peak0=0;   Cpeaks=[]

    @load data_files[1]  sessionParams
    b=boxplot(title="Total number of cases in "*wa_name*" - R0="*string(sessionParams.R₀))#*", d="*string(sessionParams.δ))
    b_v2=boxplot(title="Total number of cases in "*wa_name*" - R0="*string(sessionParams.R₀))
    b4=boxplot(title="Gain in number of cases in "*wa_name)
    b7=boxplot(title="Peak delays in "*wa_name)
    final_cum0detection=0;        cumI_diffs=[];        cumIs=[];    medians_cumI=[]
    cumCs=[];cumICs=[]
    xticks_labels=["det, dur"];
    for i=1:size(data_files,1)
        @load data_files[i]  sims_vector
        @load data_files[i]  sessionParams
        push!(xticks_labels,string(Int(sessionParams.τₚ[wa]*100))*"% "*string(Int(sessionParams.CT_dur[12]))*"days")#*@sprintf("%.0E", sessionParams.Κ_max_capacity12))
        ##final_cumI_wa:#peaks code
        final_cum=[];           peak=[];    #Cpeak=[]
        cumC=[];cumIC=[];
        for sim=1:size(sims_vector,1) #for each sim
            push!(final_cum,sum(sims_vector[sim][9][wa,:,1:3])) # final cumulative I=IA+ID+IQ, summed for all ages

            I=[sims_vector[sim][4][t][wa][1] for t=1:size(sims_vector[sim][4],1)]
            tpeak=findall(x->x==maximum(I),I)[1]
            push!(peak,tpeak)

            #C=[sims_vector[sim][2][t][wa][1] for t=1:size(sims_vector[sim][2],1)]
            #push!(Cpeak,findall(x->x==maximum(C),C)[1])
            push!(cumC,sims_vector[sim][2][end][wa][1])
            push!(cumIC,sims_vector[sim][3][end][wa][1])
        end
        push!(cumIs,median(final_cum))
        push!(cumCs,median(cumC))
        push!(cumICs,median(cumIC))

        boxplot!(b,final_cum,color=repeat([colors2[i]],5),legend=false)
        if i<=9
            boxplot!(b_v2,final_cum,color=repeat([colors2[i]],5),legend=false,ls=:dash)
        else
            boxplot!(b_v2,final_cum,color=repeat([colors2[i]],5),legend=false,ls=:solid)
        end

        if i==1     final_cum0detection=median(final_cum)       end
        push!(cumI_diffs,final_cum0detection-median(final_cum))
        boxplot!(b4,final_cum0detection .- final_cum,color=repeat([colors2[i]],5),legend=false)#,line=:dot)

        ##final_cumI_barplot:
        cumI=[[sum(sims_vector[sim][1][end][wa,1:3]) for wa=1:20]   for sim=1:size(sims_vector,1)]
        median_cumI=[median([cumI[sim][wa] for sim=1:size(sims_vector,1)])      for wa=1:20]
        push!(medians_cumI,median_cumI)

        push!(peaks,median(peak))
        if i==1     median_peak0=median(peaks)       end
        push!(peak_diffs,median(peaks)-median_peak0)
        boxplot!(b7,peaks .- median_peak0,legend=false,color=repeat([colors2[i]],5))

        #push!(Cpeaks,median(Cpeak))
    end
    ##final_cumI_wa:
    #@save plot_folder*"stats.jld2"   cumI_diffs,cumIs,peak_diffs,peaks
    CSV.write(plot_folder*"stats.csv", DataFrame(xticks_labels=xticks_labels[2:end],cumI_diffs=cumI_diffs, cumIs=cumIs,peak_diffs=peak_diffs,peaks=peaks#=,Cpeaks=Cpeaks=#))

    #b8=scatter(Cpeaks,color=colors2,title="CPeak times in "*wa_name,legend=false,
    #        xticks=([0:1:size(xticks_labels,1);], xticks_labels),size=(width2,length1),xrotation=Rotation)#,markershape)

    boxplot!(b, xticks=([0:1:size(xticks_labels,1);], xticks_labels),size=(width1,length1),xrotation=Rotation)#,ylims=(1.04e6,Inf))
        boxplot!(b_v2, xticks=([0:1:size(xticks_labels,1);], xticks_labels),size=(width1,length1),xrotation=Rotation,ylims=(6.5e5,Inf))
        #display(b_v2)

    boxplot!(b4,xticks=([0:1:size(xticks_labels,1);], xticks_labels),size=(width2,length1),xrotation=Rotation)
    boxplot!(b7,xticks=([0:1:size(xticks_labels,1);], xticks_labels),size=(width1,length1),xrotation=Rotation)

    b2=bar(cumI_diffs[2:end],color=colors2,title="Gain in number of cases in "*wa_name*" - R0="*string(sessionParams.R₀),linestyles=[[:dot for i=1:4] ;[:solid for i=1:5]],lw=0,#linewidth=false,
            legend=false,xticks=([0:1:size(xticks_labels,1);], xticks_labels),size=(width1,length1),xrotation=Rotation)
            #display(b2)
            b2_v2=bar(cumI_diffs[9:end],color=colors2,title="Gain in number of cases in "*wa_name*" - R0="*string(sessionParams.R₀),lc=false,fillalpha=.8,#linewidth=false,
                    legend=false,xticks=([0:1:size(cumI_diffs[2:7],1);], [xticks_labels[1];xticks_labels[3:end]]),size=(width1,length1),xrotation=Rotation)
            bar!(b2_v2,cumI_diffs[2:7],color=colors2,ls=:dash,lw=0,bar_widths=.5,
                    legend=false,size=(width1,length1),xrotation=Rotation)
                    #display(b2_v2)
    b3=bar(cumIs,color=colors2,title="Total number of cases in "*wa_name,
            legend=false,xticks=([0:1:size(xticks_labels,1);], xticks_labels),size=(width2,length1),xrotation=Rotation)
    b5=bar(peak_diffs,color=colors2,title="Peak time difference in "*wa_name,
            legend=false,xticks=([0:1:size(xticks_labels,1);], xticks_labels),size=(width2,length1),xrotation=Rotation)
    b6=bar(peaks,color=colors2,title="Peak times in "*wa_name,
            legend=false,xticks=([0:1:size(xticks_labels,1);], xticks_labels),size=(width2,length1),xrotation=Rotation)

    b8=bar(cumCs[9:end],color=colors2,title="Total number of contacted in "*wa_name*" - R0="*string(sessionParams.R₀),lc=false,fillalpha=.8,#linewidth=false,
            legend=false,xticks=([0:1:size(cumI_diffs[2:7],1);], [xticks_labels[1];xticks_labels[3:end]]),size=(width1,length1),xrotation=Rotation)
        bar!(b8,cumCs[2:7],color=colors2,ls=:dash,lw=0,bar_widths=.5,
                legend=false,size=(width1,length1),xrotation=Rotation);#display(b8)
    b9=bar(cumICs[9:end],color=colors2,title="Contacted infecteds (E or Ix) in "*wa_name*" - R0="*string(sessionParams.R₀),lc=false,fillalpha=.8,#linewidth=false,
            legend=false,xticks=([0:1:size(cumI_diffs[2:7],1);], [xticks_labels[1];xticks_labels[3:end]]),size=(width1,length1),xrotation=Rotation)
        bar!(b9,cumICs[2:7],color=colors2,ls=:dash,lw=0,bar_widths=.5,
                legend=false,size=(width1,length1),xrotation=Rotation);#display(b9)


    savefig(b,plot_folder*"jl_cumI_"*wa_name*"_box.png");        savefig(b2,plot_folder*"jl_cumI_gain_"*wa_name*"_bar.png");
                                                                 savefig(b2_v2,plot_folder*"jl_cumI_gain_"*wa_name*"_bar_v2.png")
    savefig(b3,plot_folder*"jl_cumI_"*wa_name*"_bar.png");       savefig(b4,plot_folder*"jl_cumI_gain_"*wa_name*"_box.png")
    savefig(b5,plot_folder*"jl_peak_diffs_"*wa_name*"_bar.png"); savefig(b6,plot_folder*"jl_peaks_"*wa_name*"_bar.png")
    savefig(b7,plot_folder*"jl_peaks_"*wa_name*"_box.png");      savefig(b8,plot_folder*"jl_cpeaks_"*wa_name*"_bar.png")

    savefig(b8,plot_folder*"jl_cumC_"*wa_name*"_bar.png");      savefig(b9,plot_folder*"jl_cumIC_"*wa_name*"_bar.png");
end
function make_plots_oneExample(sessions,τₚ_list,S0,Κ_max_capacity_Kilifi)
    wa=12;wa_name="Kilifi"
    p1=plot(title="1sim example I=A+D+IQ "*wa_name,legend=:topright);
    p2=plot(title="1sim example Cumulative I=A+D+IQ "*wa_name,legend=:topright);
    p2Q=plot(title="1sim example Cumulative IQ "*wa_name,legend=:topright);
    p3=plot(title="1sim example Q "*wa_name,legend=:topright);
    p4=plot(title="1sim example S, R Kilifi",legend=:topleft);
    p5=plot(title="1sim example infecteds quarantined after contact Kilifi",legend=:topleft);
    p6=plot();p7=plot();
    j=0
    for session in sessions
        j+=1
        println("Working plots for session ",session)
        results_folder=".\\contacts\\results_session"*string(Int(floor(session/10)))*"0s\\results_session"*string(session)*"\\"
        data_files=readdir(results_folder)
        data_files=[s for s in data_files if endswith(s,".jld2")]
        @load results_folder*data_files[1]  sims_vector
        n_traj=size(sims_vector,1)

        #one_sim_wa(τₚ_list,results_folder,data_files,n_traj,wa_coords,12,"Kilifi");#one_sim_wa(τₚ_list,results_folder,data_files,n_traj,wa_coords,4,"Nairobi")
        for i=2#1:size(τₚ_list,1)
            @load results_folder*data_files[i]  sims_vector
            I=[sims_vector[2][4][t][wa][1] for t=1:size(sims_vector[1][2],1)]
            plot!(p1,I)#,label="session"*string(session)*"-"*string(τₚ_list[i]*100)*"%-I", color=colors[j])
            cumIA=[sims_vector[2][1][t][wa,1] for t=1:size(sims_vector[1][1],1)]
            cumID=[sims_vector[2][1][t][wa,2] for t=1:size(sims_vector[1][1],1)]
            cumIQ=[sims_vector[2][1][t][wa,3] for t=1:size(sims_vector[1][1],1)]
            cumI=cumIA .+ cumID .+cumIQ #OR I=[sum(sims_vector[1][1][t][wa,1:3]) for t=1:size(sims_vector[1][1],1)]
            plot!(p2,cumI)#,label="session"*string(session)*"-"*string(τₚ_list[i]*100)*"%-cumI", color=colors[j])
            plot!(p2,cumIA)#,label="session"*string(session)*"-"*string(τₚ_list[i]*100)*"%-IA", color=colors[j], line=:dash)
            plot!(p2,cumID)#,label="session"*string(session)*"-"*string(τₚ_list[i]*100)*"%-ID", color=colors[j], line=:dot)
            plot!(p2,cumIQ)#,label="session"*string(session)*"-"*string(τₚ_list[i]*100)*"%-IQ", color=colors[j], line=(2,:dashdotdot))
            plot!(p2Q,cumIQ)#,label="session"*string(session)*"-"*string(τₚ_list[i]*100)*"%-IQ", color=colors[j], line=(2,:dashdotdot))
            Q=[sims_vector[2][6][t][wa][1] for t=1:size(sims_vector[1][2],1)]
            plot!(p3,Q)#,label="detection"*string(τₚ_list[i]*100)*"%-Q", color=colors[j])

            ##TO DELETE WHEN RUNNING LARGE SIMS:
            S=[sims_vector[1][5][t][wa][1] for t=1:size(sims_vector[1][1],1)]
            plot!(p4,S)#,label="detection"*string(τₚ_list[i]*100)*"%-S", color=colors[j])
            R=[sims_vector[1][6][t][wa][1] for t=1:size(sims_vector[1][1],1)]
            plot!(p4,R)#,label="detection"*string(τₚ_list[i]*100)*"%-R", color=colors[j])
            IQ=[sims_vector[2][5][t][wa][1] for t=1:size(sims_vector[1][1],1)];plot(IQ)
            cumIC=[sims_vector[2][3][t][wa][1] for t=1:size(sims_vector[1][1],1)]
            plot!(p5,cumIC)#,label="detection"*string(τₚ_list[i]*100)*"%-cumIC", color=colors[j])
            cumC=[sims_vector[2][2][t][wa][1] for t=1:size(sims_vector[1][1],1)]
            plot!(p6,cumC)
        end
    end
    savefig(p1,".\\contacts\\results_session70s\\"*"jl_ONE_SIM_I_"*wa_name*".png")
    savefig(p2,".\\contacts\\results_session70s\\"*"jl_ONE_SIM_CumI_"*wa_name*".png")
    savefig(p2Q,".\\contacts\\results_session70s\\"*"jl_ONE_SIM_CumIQ_"*wa_name*".png")
    savefig(p3,".\\contacts\\results_session70s\\"*"jl_ONE_SIM_Q_"*wa_name*".png")
    savefig(p4,".\\contacts\\results_session70s\\"*"jl_ONE_SIM_SR_"*wa_name*".png")
    savefig(p5,".\\contacts\\results_session70s\\"*"jl_ONE_SIM_cumIC_"*wa_name*".png")
end
########
S0=[4.138758e6, 867417.0, 2.326182e6, 8.084069e6, 3.229145e6, 459761.0, 999280.0, 1.979082e6, 926952.0, 340661.0, 2.381706e6, 2.126254e6, 2.960717e6, 786461.0, 7.478259e6, 781212.0, 2.114588e6, 4.094022e6, 569586.0, 917976.0]
S012A=[279992, 264585, 246992, 206735, 220761, 211260, 181931, 130091, 111761, 82036, 55221, 42469, 34381, 23781, 16214, 18044]
folder="./contacts/results_session10-11_10sims/"
width1=300;length1=300;
width2=500;length2=length1*2;
Rotation=40
#make_plots_oneδ(folder,S0)

folders=["./contacts/results_session10_20sims/","./contacts/results_session11_20sims/"];plot_folder="./contacts/results_session10_20sims/plots_sessions10-11/"
#folders=["./contacts/results_session12_20sims/","./contacts/results_session13_20sims/"];plot_folder="./contacts/results_session10_20sims/plots_sessions12-13/"
δs=[.3,.7]
make_plots_twoδs(folders,plot_folder,S0,δs)
#make_plots_ages(folder,S0)
#make_plots_oneExample(sessions,τₚ_list,S0,Κ_max_capacity_Kilifi)
#save_sessionParams(folder,cumIs,cumI_diffs)

#=folders=[]
data_files=readdir(folder)
data_files=[folder*s for s in data_files if endswith(s,".jld2")&&s!="stats.jld2"]
for i=1:size(data_files,1)#for i=1:size(τₚ_list,1)
    #@load data_files[i]  sims_vector
    @load data_files[i]  sessionParams
    if sessionParams.stop_Q==true
        print("sc_nb=",sessionParams.sc_nb," / ")
    end
end=#
