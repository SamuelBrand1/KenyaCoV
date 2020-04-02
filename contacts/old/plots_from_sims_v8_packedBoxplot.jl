using Plots,MAT, StatsPlots, Statistics,JLD2, Images, ImageView,ImageDraw,CSV,DataFrames,Printf,ColorBrewer
Plots.default(grid=false,legendfontsize=6, legendfontcolor=:gray30,titlefontsize=8, titlefontcolor=:gray30,xtickfontsize=7, ytickfontsize=6,tickfontcolor=:gray30,titlefontfamily="Cambria")
StatsPlots.default(grid=false,legendfontsize=6, legendfontcolor=:gray30,#:midnightblue#legendfont="Cambria"
                    titlefontsize=8, titlefontcolor=:gray30,#titlefontfamily="Cambria",
                    xtickfontsize=7, ytickfontsize=6,tickfontcolor=:gray30,linewidth=false)#tickfontfamily="Cambria",
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
    params=[]
    for i=1:size(data_files,1)
        @load data_files[i]  sessionParams
        sessionParams.cumI=cumIs[i];    sessionParams.cumI_diff=cumI_diffs[i];
        params=[params;sessionParams]
    end
    sort!(params, by = x -> x.sc_nb)
    CSV.write(folder*"sessionParams.csv", params)
    return params
end
function mysort(data_files)
    A=[]
    for i=1:size(data_files,1)
        @load data_files[i]  sessionParams
        A=[A; [[data_files[i] sessionParams.sc_nb]]]
    end
    sort!(A, by = x -> x[2]);
    return [A[i][1] for i=1:size(A,1)]
end
function make_plots_v7(folder,S0)
    data_files=readdir(folder);
    data_files=[folder*s for s in data_files if endswith(s,".jld2")&&s!="stats.jld2"]
    data_files=mysort(data_files)
    cumI_diffs=[];peakss=[];peak_diffss=[];    peaks=[];    peak_diffs=[];    median_peak0=0
    wa=12;wa_name="Kilifi"

    @load data_files[1]  sessionParams
    b=boxplot(title="Total number of cases in "*wa_name*" - R0="*string(sessionParams.R₀))#*", d="*string(sessionParams.δ))
    b4=boxplot(title="Gain in number of cases in "*wa_name)
    b7=boxplot(title="Peak delays in "*wa_name)
    final_cum0detection=0;        cumI_diffs=[];        cumIs=[];    medians_cumI=[];   Cpeak=[]
    xticks_labels=["p,k_max"];

    folder="./contacts/results_session1_50sims_NTD/session1_50sims_true/"
    data_files=readdir(folder);
    data_files=[folder*s for s in data_files if endswith(s,".jld2")&&s!="stats.jld2"]
    data_files=mysort(data_files)
    for i=1:size(data_files,1)
        @load data_files[i]  sims_vector
        @load data_files[i]  sessionParams
        push!(xticks_labels,string(Int(sessionParams.τₚ*100))*"%,"*@sprintf("%.0E", sessionParams.Κ_max_capacity12))
        ##final_cumI_wa:#peaks code
        final_cum=[];           peak=[]
        for sim=1:size(sims_vector,1) #for each sim
            push!(final_cum,sum(sims_vector[sim][4][wa,:,1:3])) # final cumulative I=IA+ID+IQ, summed for all ages

            #contacts:
            C=[sims_vector[sim][end][t][wa][1] for t=1:size(sims_vector[sim][3],1)]
            push!(Cpeak,findall(x->x==maximum(C),C)[1])
        end
        push!(cumIs,median(final_cum))
        if i==1 boxplot!(b,final_cum,color=repeat([colors2[1]],5),legend=false)
        else    boxplot!(b,final_cum,color=repeat([colors2[2]],5),legend=false,shape=:circle,bar_width=.5) end

        push!(Cpeaks,median(Cpeak))

        #=if i==1     final_cum0detection=median(final_cum)       end
        push!(cumI_diffs,final_cum0detection-median(final_cum))
        boxplot!(b4,final_cum0detection .- final_cum,color=repeat([colors2[i]],5),legend=false)

        ##final_cumI_barplot:
        cumI=[[sum(sims_vector[sim][1][end][wa,1:3]) for wa=1:20]   for sim=1:size(sims_vector,1)]
        median_cumI=[median([cumI[sim][wa] for sim=1:size(sims_vector,1)])      for wa=1:20]
        push!(medians_cumI,median_cumI)

        push!(peaks,median(peak))
        if i==1     median_peak0=median(peaks)       end
        push!(peak_diffs,median(peaks)-median_peak0)
        boxplot!(b7,peaks .- median_peak0,legend=false,color=repeat([colors2[i]],5))=#
    end


    final_cum0detection=0;        cumI_diffs=[];        cumIs=[];    medians_cumI=[]
    folder="./contacts/results_session2_50sims_NTD/session2_50sims_true/"
    data_files=readdir(folder);
    data_files=[folder*s for s in data_files if endswith(s,".jld2")&&s!="stats.jld2"]
    data_files=mysort(data_files)
    for i=2:size(data_files,1)
        @load data_files[i]  sims_vector
        final_cum=[];
        for sim=1:size(sims_vector,1) #for each sim
            push!(final_cum,sum(sims_vector[sim][4][wa,:,1:3])) # final cumulative I=IA+ID+IQ, summed for all ages
        end
        push!(cumIs,median(final_cum))
        boxplot!(b,[i+.05],final_cum,color=repeat([colors2[5]],5),legend=false,bar_width=.5)
    end

    final_cum0detection=0;        cumI_diffs=[];        cumIs=[];    medians_cumI=[]
    folder="./contacts/results_session3_50sims_NTD/session3_50sims_true/"
    data_files=readdir(folder);
    data_files=[folder*s for s in data_files if endswith(s,".jld2")&&s!="stats.jld2"]
    data_files=mysort(data_files)
    for i=2:size(data_files,1)
        @load data_files[i]  sims_vector
        final_cum=[];
        for sim=1:size(sims_vector,1) #for each sim
            push!(final_cum,sum(sims_vector[sim][4][wa,:,1:3])) # final cumulative I=IA+ID+IQ, summed for all ages
        end
        push!(cumIs,median(final_cum))
        boxplot!(b,[i+.2],final_cum,color=repeat([colors2[4]],5),legend=false,bar_width=.5)
    end




    ##final_cumI_wa:
    #@save folder*"stats.jld2"   cumI_diffs,cumIs,peak_diffs,peaks
    #CSV.write(folder*"stats.csv", DataFrame(xticks_labels=xticks_labels[2:end],cumI_diffs=cumI_diffs, cumIs=cumIs,peak_diffs=peak_diffs,peaks=peaks))

    boxplot!(b, xticks=([0:1:size(xticks_labels,1);], xticks_labels),size=(width1,length1),xrotation=Rotation,ylims=(1.04e6,Inf))
    display(b)
    #=boxplot!(b4,xticks=([0:1:size(xticks_labels,1);], xticks_labels),size=(width2,length1),xrotation=Rotation)
    boxplot!(b7,xticks=([0:1:size(xticks_labels,1);], xticks_labels),size=(width1,length1),xrotation=Rotation)#,color=colors2)
    b2=bar(cumI_diffs,color=colors2,title="Gain in number of cases in "*wa_name*" - R0="*string(sessionParams.R₀)*", d="*string(sessionParams.δ),linewidth=false,
            legend=false,xticks=([0:1:size(xticks_labels,1);], xticks_labels),size=(width1,length1),xrotation=Rotation)
    b3=bar(cumIs,color=colors2,title="Total number of cases in "*wa_name,
            legend=false,xticks=([0:1:size(xticks_labels,1);], xticks_labels),size=(width2,length1),xrotation=Rotation)
    b5=bar(peak_diffs,color=colors2,title="Peak time difference in "*wa_name,
            legend=false,xticks=([0:1:size(xticks_labels,1);], xticks_labels),size=(width2,length1),xrotation=Rotation)
    b6=bar(peaks,color=colors2,title="Peak times in "*wa_name,
            legend=false,xticks=([0:1:size(xticks_labels,1);], xticks_labels),size=(width2,length1),xrotation=Rotation)=#
    savefig(b,folder*"jl_cumI_"*wa_name*"_box2.png")
    #=savefig(b4,folder*"jl_cumI_gain_"*wa_name*"_box.png")
    savefig(b2,folder*"jl_cumI_gain_"*wa_name*"_bar.png")
    savefig(b3,folder*"jl_cumI_"*wa_name*"_bar.png")
    savefig(plot(b,b2,layout=(2,1),size=(width1,length1*2)),folder*"jl_2plots_"*wa_name*".png")
    savefig(plot(b,b2,layout=(1,2),size=(width1*2,length1)),folder*"jl_2plots_"*wa_name*"_2.png")
    savefig(b5,folder*"jl_peak_diffs_"*wa_name*"_bar.png")
    savefig(b6,folder*"jl_peaks_"*wa_name*"_bar.png")
    savefig(b7,folder*"jl_peaks_"*wa_name*"_box.png")
    ##final_cumI_barplot:
    S02=[]; for wa=1:20,i=1:size(data_files,1) push!(S02,S0[wa]);   end
    barplotALL=bar(S02,color=:green,xticks=([1:size(data_files,1):20*size(data_files,1);],riskregionnames),size=(width2,length1),xrotation=Rotation,legend=false)

    medians_cumI2=[];
    for wa=1:20,i=1:size(data_files,1)    push!(medians_cumI2,medians_cumI[i][wa]);   end
    bar!(barplotALL,medians_cumI2,color=colors2)
    savefig(barplotALL,folder*"jl_cumI_barplot.png")# =#

        ##FOR cumIgain_many_sessions and peaks_many_sessions
        #push!(cumI_diffs,get_cumIgain_wa(τₚ_list,results_folder,data_files,n_traj,wa_coords,12,"Kilifi"))
        #peaks,peak_diffs=get_peaks_wa(τₚ_list,results_folder,data_files,n_traj,wa_coords,12,"Kilifi")
        #push!(peakss,peaks);push!(peak_diffss,peak_diffs)

    #cumIgain_many_sessions(cumI_diffs,["500","1e3","5e3","1e4"],sessions,τₚ_list,12,"Kilifi")
    #peaks_many_sessions(peakss,peak_diffss,["500","1e3","5e3","1e4"],sessions,τₚ_list,12,"Kilifi")
    #return peakss,peak_diffss
end

function make_plots_ages(folder,S0)
    data_files=readdir(folder);
    data_files=[folder*s for s in data_files if endswith(s,".jld2")&&s!="stats.jld2"]
    data_files=mysort(data_files)
    wa=12;wa_name="Kilifi"

    @load data_files[1]  sims_vector
    cumIsA=[[] for i=1:16]
    for sim=1:size(sims_vector,1),a=1:16    push!(cumIsA[a],sum(sims_vector[sim][4][wa,a,1:3]));end

    @load data_files[5]  sims_vector# 9 is sc53
    cumIsA_sc53=[[] for i=1:16]
    for sim=1:size(sims_vector,1),a=1:16    push!(cumIsA_sc53[a],sum(sims_vector[sim][4][wa,a,1:3])); end

    @load data_files[9]  sims_vector# 9 is sc58
    cumIsA_sc58=[[] for i=1:16]
    for sim=1:size(sims_vector,1),a=1:16    push!(cumIsA_sc58[a],sum(sims_vector[sim][4][wa,a,1:3])); end

    S012A_2=[]
    for i=1:16  push!(S012A_2,S012A[i]);push!(S012A_2,S012A[i]) ;end
    b=bar(S012A_2,bar_width=1.01,color=:green,xticks=([1:2:32;],["A"*string(i) for i=1:16]),fillalpha=.5,linecolor=false,xrotation=Rotation,size=(width2,length1))#bar_edges=false)

    cumIfinal=[];
    for a=1:16  push!(cumIfinal,median(cumIsA[a]));push!(cumIfinal,median(cumIsA_sc58[a])); end
    bar!(b,cumIfinal,color=repeat(colors2[1:2],16),xrotation=Rotation,contour_labels=true,size=(width2,length1))

    ##Grouped barplot number infecteds
    cumIs2 = [median(cumIsA[a]) for a=1:16]#[20, 35, 30, 35, 27,25, 32, 34, 20, 25]
    append!(cumIs2,[median(cumIsA_sc58[a]) for a=1:16])
    append!(cumIs2,[median(cumIsA_sc53[a]) for a=1:16])
    sx = repeat(["No detection", "sc58(50%,1e5,false)", "sc53(50%,1e5,true)"], inner = 16)
    xlabels = repeat(string.(1:16), outer = 3)

    gb=bar(S012A,color=:green,fillalpha=.5,linecolor=false,label="Age group size")
    groupedbar!(gb, cumIs2, group = sx, ylabel = "Final infected",linecolor=false,color=[colors2[2] colors2[5] colors2[1]],
    title = "Final cumulative infecteds in Kilifi by age group",xticks=false#=([1:1:16;],ages)=#,xrotation=Rotation,size=(width2,length1));

    ##Grouped barplot %infecteds
    cumIs3 = [median(cumIsA[a])*100/S012A[a] for a=1:16]#[20, 35, 30, 35, 27,25, 32, 34, 20, 25]
    append!(cumIs3,[median(cumIsA_sc58[a])*100/S012A[a] for a=1:16])
    append!(cumIs3,[median(cumIsA_sc53[a])*100/S012A[a] for a=1:16])
    sx = repeat(["No detection", "sc58(50%,1e5,false)", "sc53(50%,1e5,true)"], inner = 16)
    xlabels = repeat(string.(1:16), outer = 3)

    #gb=bar(S012A,color=:green,fillalpha=.5,linecolor=false,label="Age group size")
    gb2=groupedbar(cumIs3, group = sx, ylabel = "Final infected",linecolor=false,color=[colors2[2] colors2[5] colors2[1]],
    title = "% infecteds in Kilifi by age group",xticks=([1:1:16;],ages),xrotation=Rotation,size=(width2,length1),legend=false);

    savefig(plot(gb,gb2,layout=(2,1),size=(width2,length2)),folder*"jl_cumI_"*wa_name*"_groupedbar.png")
end
########
S0=[4.138758e6, 867417.0, 2.326182e6, 8.084069e6, 3.229145e6, 459761.0, 999280.0, 1.979082e6, 926952.0, 340661.0, 2.381706e6, 2.126254e6, 2.960717e6, 786461.0, 7.478259e6, 781212.0, 2.114588e6, 4.094022e6, 569586.0, 917976.0]
S012A=[279992, 264585, 246992, 206735, 220761, 211260, 181931, 130091, 111761, 82036, 55221, 42469, 34381, 23781, 16214, 18044]
folder="./contacts/results_session1_50sims_NTD/session1_50sims_false/"
width1=400;length1=300;
width2=500;length2=length1*2;
Rotation=40
make_plots_v7(folder,S0)
#make_plots_ages(folder,S0)
#make_plots_oneExample(sessions,τₚ_list,S0,Κ_max_capacity_Kilifi)
#save_sessionParams(folder,cumIs,cumI_diffs)

folders=[]
#=
data_files=readdir(folder)
data_files=[folder*s for s in data_files if endswith(s,".jld2")&&s!="stats.jld2"]
for i=1:size(data_files,1)#for i=1:size(τₚ_list,1)
    #@load data_files[i]  sims_vector
    @load data_files[i]  sessionParams
    if sessionParams.stop_Q==true
        print("sc_nb=",sessionParams.sc_nb," / ")
    end
end=#
