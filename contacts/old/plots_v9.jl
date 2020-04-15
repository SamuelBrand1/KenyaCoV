using Plots,MAT, StatsPlots, Statistics,JLD2,CSV,DataFrames,Printf,ColorBrewer#, Images, ImageView,ImageDraw
Plots.default(grid=false,legendfontsize=6, legendfontcolor=:gray30,titlefontsize=8, titlefontcolor=:gray30,xtickfontsize=6, ytickfontsize=6,tickfontcolor=:gray30,titlefontfamily="Cambria")
StatsPlots.default(grid=false,legendfontsize=5, legendfontcolor=:gray30,#:midnightblue#legendfont="Cambria"
                    titlefontsize=8, titlefontcolor=:gray30,#titlefontfamily="Cambria",
                    xtickfontsize=6, ytickfontsize=6,tickfontcolor=:gray30,linewidth=false)#tickfontfamily="Cambria",
colors=[:orange,:purple3,:maroon,:gold,:orangered,:grey,:purple,:ivory3,:chocolate1,:tan1,:rosybrown,:rosybrown2,:brown2,:brown3,:brown4,:deeppink3,:deeppink4]
colors2=ColorBrewer.palette("Pastel1",8);colors2=repeat(colors2,outer=[10])
markershapes=[:cross,:star4, :vcross, :star6, :hline]
riskregionnames = ["Machakos/Muranga" "Mandera" "Baringo/Nakuru" "Nairobi" "Uasin Gishu" "Marsabit" "Garissa" "Nakuru/Narok" "Turkana" "Taita Taveta" "Kitui/Meru" "Kilifi/Mombasa" "Kericho/Kisumu" "Kilifi/Lamu" "Kakamega/Kisumu" "Wajir" "Kajiado/Kisumu" "Homa bay/Migori" "Samburu/Laikipia" "Kilifi/Kwale" "Total"]
wa_coords=[[300,450], [515,85], [165,360], [235,465], [115,300], [300,140], [510,380], [180, 430], [120, 130], [355, 615], [340, 375], [465, 630], [100, 380], [495, 530], [40, 355], [490, 255], [155, 495], [30, 440], [250, 290], [400, 670]]
ages=[string(i)*"-"*string(i+4) for i=0:5:70];push!(ages,"75+yrs")#"0-4" "5-9" "10-14" "15-19" "20-24" "25-29" "30-34" "35-39" "40-44" "45-49" ]
######## No intervention plots
function kenya_barplot(data_files,plot_folder,S0)
    if !isdir(plot_folder)   mkdir(plot_folder)   end
    wa=12;wa_name="Kilifi"

    final_cum0detection=0;    #=cumIs=[];=#   cumI_medians=[];    cumI_medians_kenya=[];
    cumDeaths_medians=[];cumDeaths_medians_kenya=[];cumDeaths_stds_kenya=[]
    for i=1:size(data_files,1)
        @load data_files[i]  sims_vector
        @load data_files[i]  sessionParams
        #=final_cum=[];
        for sim=1:size(sims_vector,1) #for each sim
            push!(final_cum,sum(sims_vector[sim][10][wa,:,1:3])) # final cumulative I=IA+ID+IQ, summed for all ages
        end
        push!(cumIs,median(final_cum));=#

        ##final_cumI_barplot:
        cumI=[[sum(sims_vector[sim][1][end][wa,1:3]) for wa=1:20]   for sim=1:size(sims_vector,1)]
        median_cumI=[median([cumI[sim][wa] for sim=1:size(sims_vector,1)])      for wa=1:20]
        push!(cumI_medians,median_cumI)

        ##deaths for barplot:
        cumDeaths=[[sum(sims_vector[sim][4][end][wa]) for wa=1:20]   for sim=1:size(sims_vector,1)]
        median_cumDeaths=[median([cumDeaths[sim][wa] for sim=1:size(sims_vector,1)])      for wa=1:20]
        push!(cumDeaths_medians,median_cumDeaths)

        ##All kenya cumI
        median_cumI_kenya=median([sum([cumI[sim][wa]  for wa=1:20])  for sim=1:size(sims_vector,1)])
        push!(cumI_medians_kenya,median_cumI_kenya)
        std_cumI_kenya=std([sum([cumI[sim][wa]  for wa=1:20])  for sim=1:size(sims_vector,1)])
        push!(cumDeaths_stds_kenya,std_cumI_kenya)

        ## All Kenya cumDeaths
        median_cumDeaths_kenya=median([sum([cumDeaths[sim][wa]  for wa=1:20])  for sim=1:size(sims_vector,1)])
        push!(cumDeaths_medians_kenya,median_cumDeaths_kenya)

    end
    S02=[]; for wa=1:20 push!(S02,S0[wa]);   end
    barplotALL=bar([1.5:2:40.5;],S02,color=[:green,:chartreuse4],xticks=([1:2:40;],riskregionnames),fillalpha=.5,linewidth=false,bar_width=2,linecolor=false,
                    size=(width1,length1),rotation=40,label=false,title="Estimated cases in Kenya without intervention")

    cumI_medians2=[];cumDeaths_medians2=[]
    for wa=1:20,i=1:1size(data_files,1)    push!(cumI_medians2,cumI_medians[i][wa]);   push!(cumDeaths_medians2,cumDeaths_medians[i][wa]);end

    bar!(barplotALL,cumI_medians2,color=[colors2[1],colors2[6]],linecolor=false,label=false,bar_width=1.005)#,yerrorbar=cumI_allfiles_sd,markercolor=:black)
    bar!(barplotALL,cumDeaths_medians2,color=:red,linecolor=false,label=false,bar_width=1.005);#display(barplotALL)
    savefig(barplotALL,plot_folder*"jl_cumI_barplot.png")

    bdeaths=bar(cumDeaths_medians2,color=:red,yflip=:true,xaxis=:false,size=(width1*.95,length1*.4),label=false);#display(bdeaths)
    savefig(bdeaths,plot_folder*"jl_cumDeaths_Kenya_barplot.png")
    return cumI_medians_kenya,cumDeaths_stds_kenya,cumI_medians,cumDeaths_medians
end

function kilifi_ages(data_files,plot_folder,S0)
    if !isdir(plot_folder)   mkdir(plot_folder)   end
    wa=12;wa_name="Kilifi"

    @load data_files[1]  sims_vector
    cumIsA=[[] for i=1:16];     cumDeaths_1=[[] for i=1:16];
    for sim=1:size(sims_vector,1),a=1:16    push!(cumIsA[a],sum(sims_vector[sim][10][wa,a,1:3]));end
    for sim=1:size(sims_vector,1),a=1:16    push!(cumDeaths_1[a],sims_vector[sim][11][wa,a]); end

    @load data_files[2]  sims_vector#
    cumIsA_2=[[] for i=1:16];   cumDeaths_2=[[] for i=1:16];
    for sim=1:size(sims_vector,1),a=1:16    push!(cumIsA_2[a],sum(sims_vector[sim][10][wa,a,1:3])); end
    for sim=1:size(sims_vector,1),a=1:16    push!(cumDeaths_2[a],sims_vector[sim][11][wa,a]); end

    S012A_2=[]
    for i=1:16  push!(S012A_2,S012A[i]);push!(S012A_2,S012A[i]) ;end
    b=bar(S012A_2,bar_width=1.01,color=[:green,:chartreuse4],xticks=([1.5:2:32;],ages),fillalpha=.5,linecolor=false,rotation=30,label=false)#bar_edges=false)

    cumIfinal=[];
    for a=1:16  push!(cumIfinal,median(cumIsA[a]));push!(cumIfinal,median(cumIsA_2[a])); end
    bar!(b,cumIfinal,color=[colors2[1],colors2[6]],rotation=30,contour_labels=true,size=(width1*1.2,length1*.8),
            linecolor=false,label=false,title="Estimated cases and deaths per age in Kilifi without intervention")

    cumDeathsfinal=[];
    for a=1:16  push!(cumDeathsfinal,median(cumDeaths_1[a]));push!(cumDeathsfinal,median(cumDeaths_2[a])); end
    fillalphas=[cumDeathsfinal[a]==0 ? 1 : 0  for a=1:16]
    bar!(b,cumDeathsfinal,color=:red,rotation=30,size=(width1*1.2,length1*.8),#fillalpha=fillalphas,contour_labels=true,
            linecolor=false,label=false)#,title="Estimated cases per age in Kilifi without intervention")
    #        scatter!(b,cumDeathsfinal,color=:red,markershapes=:hline,fillalpha=fillalphas,#rotation=30,contour_labels=true,#size=(width1*1.2,length1*.8),
    #                #=linecolor=false,=#label=false,text=[[string(cumDeathsfinal[a]) for a=1:16] [string(cumDeathsfinal[a]) for a=1:16] [string(cumDeathsfinal[a]) for a=1:16]])
    savefig(b,plot_folder*"jl_cumI_ages_kilifi.png")

    #plotting deaths
    bdeaths=bar(cumDeathsfinal,color=:red,linecolor=[colors2[1],colors2[6]],yflip=:true,xaxis=:false,size=(width1*1.15,length1*.3),label=false);#display(bdeaths)
    savefig(bdeaths,plot_folder*"jl_cumDeaths_ages_kilifi.png")
    #=l=@layout [ a{0.75h}
                b{0.25h}  ]
    savefig(plot(b,bdeaths,layout=l,size=(width1,length1)),plot_folder*"jl_cumIandcumDeaths_ages_kilifi.png")=#
    return [cumIsA,cumIsA_2]
    #=##Grouped barplot number infecteds
    cumIs2 = [median(cumIsA[a]) for a=1:16]#[20, 35, 30, 35, 27,25, 32, 34, 20, 25]
    append!(cumIs2,[median(cumIsA_2[a]) for a=1:16])
    sx = repeat(["No detection", "R0=2.5", "R0=3"], inner = 16)
    xlabels = repeat(string.(1:16), outer = 3)

    gb=bar(S012A,color=:green,fillalpha=.5,linecolor=false,label="Age group size")
    groupedbar!(gb, cumIs2, group = sx, ylabel = "Final infected",linecolor=false,color=[colors2[2] colors2[5] colors2[1]],
    title = "Final cumulative infecteds in Kilifi by age group",xticks=false#=([1:1:16;],ages)=#,rotation=30,size=(width2,length1));

    ##Grouped barplot %infecteds
    cumIs3 = [median(cumIsA[a])*100/S012A[a] for a=1:16]#[20, 35, 30, 35, 27,25, 32, 34, 20, 25]
    append!(cumIs3,[median(cumIsA_sc58[a])*100/S012A[a] for a=1:16])
    append!(cumIs3,[median(cumIsA_sc53[a])*100/S012A[a] for a=1:16])
    sx = repeat(["No detection", "R0=2.5", "R0=3"], inner = 16)
    xlabels = repeat(string.(1:16), outer = 3)

    #gb=bar(S012A,color=:green,fillalpha=.5,linecolor=false,label="Age group size")
    gb2=groupedbar(cumIs3, group = sx, ylabel = "Final infected",linecolor=false,color=[colors2[2] colors2[5] colors2[1]],
    title = "% infecteds in Kilifi by age group",xticks=([1:1:16;],ages),rotation=30,size=(width2,length1),legend=false);
    display(gb2)
    savefig(plot(gb,gb2,layout=(2,1),size=(width2,length2)),plot_folder*"jl_cumI_"*wa_name*"_groupedbar.png")=#
end

######## ???
function kilifi_ages_α5_90days(data_files,plot_folder,S0)
    if !isdir(plot_folder)   mkdir(plot_folder)   end
    wa=12;wa_name="Kilifi"

    @load data_files[1]  sims_vector
    cumIsA=[[] for i=1:16];    for sim=1:size(sims_vector,1),a=1:16    push!(cumIsA[a],sum(sims_vector[sim][10][wa,a,1:3]));end

    @load data_files[2]  sims_vector#
    cumIsA_2=[[] for i=1:16];    for sim=1:size(sims_vector,1),a=1:16    push!(cumIsA_2[a],sum(sims_vector[sim][10][wa,a,1:3])); end

    @load data_files[3]  sims_vector#
    cumIsA_3=[[] for i=1:16];    for sim=1:size(sims_vector,1),a=1:16    push!(cumIsA_3[a],sum(sims_vector[sim][10][wa,a,1:3])); end

    S012A_2=[]
    for i=1:16  push!(S012A_2,S012A[i]);push!(S012A_2,S012A[i]) ;end
    b=bar(S012A_2,bar_width=1.01,color=[:green,:chartreuse4],xticks=([1:2:32;],ages),fillalpha=.5,linecolor=false,rotation=30,label=false)#bar_edges=false)

    cumIfinal=[];
    for a=1:16  push!(cumIfinal,median(cumIsA[a]));push!(cumIfinal,median(cumIsA_2[a])); end
    bar!(b,cumIfinal,color=color=[colors2[1],colors2[6]],rotation=30,contour_labels=true,size=(width1*1.2,length1*.8),
            linecolor=false,label=false,title="Estimated cases per age in Kilifi R0=2.5,a=50%,90days")

            cumIfinal3=[];
            for a=1:16  push!(cumIfinal3,median(cumIsA_3[a])); end
    bar!(b,[2:2:32;],cumIfinal3,color=color=colors2[6],rotation=30,bar_width=.5,lw=.5,#contour_labels=true,size=(width1*1.2,length1*.8),
            ls=:dash,label=false)#,title="Estimated cases per age in Kilifi without intervention")
    #display(b)
    savefig(b,plot_folder*"jl_cumI_ages_kilifi.png")
    return cumIsA
end

function kilifi_ages_α5_9_90days(data_files,plot_folder,S0)
    if !isdir(plot_folder)   mkdir(plot_folder)   end
    wa=12;wa_name="Kilifi"

    @load data_files[1]  sims_vector
    cumIsA=[[] for i=1:16];    for sim=1:size(sims_vector,1),a=1:16    push!(cumIsA[a],sum(sims_vector[sim][10][wa,a,1:3]));end

    @load data_files[2]  sims_vector#
    cumIsA_2=[[] for i=1:16];    for sim=1:size(sims_vector,1),a=1:16    push!(cumIsA_2[a],sum(sims_vector[sim][10][wa,a,1:3])); end

    @load data_files[3]  sims_vector#
    cumIsA_3=[[] for i=1:16];    for sim=1:size(sims_vector,1),a=1:16    push!(cumIsA_3[a],sum(sims_vector[sim][10][wa,a,1:3])); end

    @load data_files[4]  sims_vector#
    cumIsA_4=[[] for i=1:16];    for sim=1:size(sims_vector,1),a=1:16    push!(cumIsA_4[a],sum(sims_vector[sim][10][wa,a,1:3])); end

    @load data_files[5]  sims_vector#
    cumIsA_5=[[] for i=1:16];    for sim=1:size(sims_vector,1),a=1:16    push!(cumIsA_5[a],sum(sims_vector[sim][10][wa,a,1:3])); end


    S012A_2=[]
    for i=1:16  push!(S012A_2,S012A[i]);push!(S012A_2,S012A[i]) ;push!(S012A_2,S012A[i]) ;end
    b=bar(S012A_2,bar_width=1.01,color=[:green,:green,:green,:chartreuse4,:chartreuse4,:chartreuse4],xticks=([1.5:3:48;],ages),fillalpha=.5,linecolor=false,rotation=30,label=false)#bar_edges=false)

    cumIfinal=[];
    for a=1:16  push!(cumIfinal,median(cumIsA[a]));push!(cumIfinal,median(cumIsA_2[a]));push!(cumIfinal,median(cumIsA_3[a])); end
    bar!(b,cumIfinal,color=[colors2[1],colors2[6],colors2[4]],rotation=30,contour_labels=true,size=(width1*1.5,length1*.8),
            ls=:solid,bar_width=1.01,linecolor=[false,:black,:black],label=false,title="Estimated cases per age in Kilifi R0=2.5,a=50%,90days")

            cumIfinal4=[];
            for a=1:16  push!(cumIfinal4,median(cumIsA_4[a])); end
    bar!(b,[2:3:48;],cumIfinal4,color=colors2[6],rotation=30,bar_width=.5,lw=.5,ls=:dash,label=false)#,title="Estimated cases per age in Kilifi without intervention")

            cumIfinal5=[];
            for a=1:16  push!(cumIfinal5,median(cumIsA_5[a])); end
    bar!(b,[3:3:48;],cumIfinal5,color=colors2[4],rotation=30,bar_width=.5,lw=.5,ls=:dash,label=false)
    #display(b)
    savefig(b,plot_folder*"jl_cumI_ages_kilifi_2.png")

    ##### Over 60s
    S012A_2=[]
    for i=13:16  push!(S012A_2,S012A[i]);push!(S012A_2,S012A[i]) ;push!(S012A_2,S012A[i]) ;end
    b=bar(S012A_2,bar_width=1.01,color=[:green,:green,:green,:chartreuse4,:chartreuse4,:chartreuse4],xticks=([1.5:3:12;],ages[13:end]),fillalpha=.5,linecolor=false,rotation=30,label=false)#bar_edges=false)

    cumIfinal=[];
    for a=13:16  push!(cumIfinal,median(cumIsA[a]));push!(cumIfinal,median(cumIsA_2[a]));push!(cumIfinal,median(cumIsA_3[a])); end
    bar!(b,cumIfinal,color=[colors2[1],colors2[6],colors2[4]],rotation=30,contour_labels=true,size=(width1*1.5,length1*.8),
            ls=:solid,bar_width=1.01,linecolor=[false,:black,:black],label=false,title="Death numbers per age in Kilifi R0=2.5,90days")

            cumIfinal4=[];
            for a=13:16  push!(cumIfinal4,median(cumIsA_4[a])); end
    bar!(b,[2:3:12;],cumIfinal4,color=colors2[6],rotation=30,bar_width=.5,lw=.5,ls=:dash,label=false)#,title="Estimated cases per age in Kilifi without intervention")

            cumIfinal5=[];
            for a=13:16  push!(cumIfinal5,median(cumIsA_5[a])); end
    bar!(b,[3:3:12;],cumIfinal5,color=colors2[4],rotation=30,bar_width=.5,lw=.5,ls=:dash,label=false);#display(b)
    savefig(b,plot_folder*"jl_cumI_ages_kilifi_over60s.png")

    ##### Over 60s GAIN
    #=S012A_2=[]
    for i=13:16  push!(S012A_2,S012A[i]);push!(S012A_2,S012A[i]) ;push!(S012A_2,S012A[i]) ;end
    b=bar(S012A_2,bar_width=1.01,color=[:green,:green,:green,:chartreuse4,:chartreuse4,:chartreuse4],xticks=([1.5:3:12;],ages[13:end]),fillalpha=.5,linecolor=false,rotation=30,label=false)#bar_edges=false)
    =#
    cumIfinalgain=[];
    for a=13:16  push!(cumIfinalgain,median(cumIsA[a])-median(cumIsA_2[a]));push!(cumIfinalgain,median(cumIsA[a])-median(cumIsA_3[a])); end
    b=bar(cumIfinalgain,color=[colors2[6],colors2[4]],rotation=30,contour_labels=true,size=(width1*1.5,length1*.8),
            ls=:solid,bar_width=1.01,linecolor=[false,:black,:black],label=false,title="Gain in death numbers per age in Kilifi R0=2.5,90days")

            cumIfinal4gain=[];
            for a=13:16  push!(cumIfinal4gain,median(cumIsA[a])-median(cumIsA_4[a])); end
    bar!(b,[1:2:8;],cumIfinal4gain,color=colors2[6],rotation=30,bar_width=.5,lw=.5,ls=:dash,label=false)#,title="Estimated cases per age in Kilifi without intervention")

            cumIfinal5gain=[];
            for a=13:16  push!(cumIfinal5gain,median(cumIsA[a])-median(cumIsA_5[a])); end
    bar!(b,[2:2:8;],cumIfinal5gain,color=colors2[4],rotation=30,bar_width=.5,lw=.5,ls=:dash,label=false);#display(b)
    savefig(b,plot_folder*"jl_cumI_ages_kilifi_over60sGAIN.png")
    return cumIsA
end

######## INTERVENTION
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
function make_plots_twoδs(folders,plot_folder,S0,δs)
    if !isdir(plot_folder)   mkdir(plot_folder)   end
    wa=12;wa_name="Kilifi"
    data_files=[]
    for folder in folders
        data_files2=readdir(folder);
        data_files2=[folder*s for s in data_files2 if endswith(s,".jld2")&&s!="stats.jld2"]
        data_files=[data_files ;data_files2]
    end
    data_files=mysort(data_files)
    @load data_files[1]  sessionParams
    b=boxplot(title="Estimated cases in "*wa_name*" - R0="*string(sessionParams.R₀))
    b_v2=boxplot(title="Estimated cases in "*wa_name*" - R0="*string(sessionParams.R₀))
    b4=boxplot(title="Gain in number of cases in "*wa_name)
    b7=boxplot(title="Peak delays in "*wa_name)

    cumI_diffs=[];peakss=[];peak_diffss=[];    peaks=[];    peak_diffs=[];    median_peak0=0;   Cpeaks=[]
    final_cum0detection=0;        cumI_diffs=[];        cumIs=[];    medians_cumI=[]
    cumCs=[];cumICs=[];           cumDeaths=[]

    xticks_labels=["a, dur"];
    for i=1:size(data_files,1)
        @load data_files[i]  sims_vector;        @load data_files[i]  sessionParams
        push!(xticks_labels,string(Int(sessionParams.τₚ[wa]*100))*"% "*string(Int(sessionParams.CT_dur[12]))*"days")
        ##final_cumI_wa:#peaks code
        final_cum=[];           peak=[];    #Cpeak=[]
        cumC=[];cumIC=[];cumDeath=[]
        for sim=1:size(sims_vector,1) #for each sim
            push!(final_cum,sum(sims_vector[sim][10][wa,:,1:3])) # final cumulative I=IA+ID+IQ, summed for all ages

            I=[sims_vector[sim][5][t][wa][1] for t=1:size(sims_vector[sim][5],1)]
            tpeak=findall(x->x==maximum(I),I)[1];   push!(peak,tpeak)

            #C=[sims_vector[sim][2][t][wa][1] for t=1:size(sims_vector[sim][2],1)]
            #push!(Cpeak,findall(x->x==maximum(C),C)[1])
            push!(cumC,sims_vector[sim][2][end][wa][1]);    push!(cumIC,sims_vector[sim][3][end][wa][1])
            push!(cumDeath,sims_vector[sim][4][end][wa][1])
        end
        push!(cumIs,median(final_cum)); push!(cumCs,median(cumC));  push!(cumICs,median(cumIC));    push!(cumDeaths,median(cumDeath))

        boxplot!(b,final_cum,color=repeat([colors2[i]],5),legend=false)
        if i<=9 boxplot!(b_v2,final_cum,color=repeat([colors2[i]],5),legend=false,ls=:dash)
        else    boxplot!(b_v2,final_cum,color=repeat([colors2[i]],5),legend=false,ls=:solid);   end

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
    CSV.write(plot_folder*"stats.csv", DataFrame(xticks_labels=xticks_labels[2:end],cumI_diffs=cumI_diffs, cumIs=cumIs,
                peak_diffs=peak_diffs,peaks=peaks#=,Cpeaks=Cpeaks=#,
                cumDeaths=cumDeaths))

    #b8=scatter(Cpeaks,color=colors2,title="CPeak times in "*wa_name,legend=false,
    #        xticks=([0:1:size(xticks_labels,1);], xticks_labels),size=(width2,length1),xrotation=Rotation)#,markershape)

    boxplot!(b, xticks=([0:1:size(xticks_labels,1);], xticks_labels),size=(width1,length1),xrotation=Rotation)#,ylims=(1.04e6,Inf))
        boxplot!(b_v2, xticks=([0:1:size(xticks_labels,1);], xticks_labels),size=(width1,length1),xrotation=Rotation,ylims=(6.5e5,Inf))

    boxplot!(b4,xticks=([0:1:size(xticks_labels,1);], xticks_labels),size=(width1,length1),xrotation=Rotation)
    boxplot!(b7,xticks=([0:1:size(xticks_labels,1);], xticks_labels),size=(width1,length1),xrotation=Rotation)

    b2=bar(cumI_diffs[2:end],color=colors2,title="Gain in number of cases in "*wa_name*" - R0="*string(sessionParams.R₀),linestyles=[[:dot for i=1:4] ;[:solid for i=1:5]],lw=0,#linewidth=false,
            legend=false,xticks=([0:1:size(xticks_labels,1);], xticks_labels),size=(width1,length1),xrotation=Rotation)
        b2_v2=bar(cumI_diffs[9:end],color=colors2,ls=:solid,title="Gain in number of cases in "*wa_name*" - R0="*string(sessionParams.R₀),#=lc=false,=#fillalpha=.8,#linewidth=false,
                legend=false,xticks=([0:1:size(cumI_diffs[2:7],1);], [xticks_labels[1];xticks_labels[3:end]]),size=(width1*.4,length1*.4),xrotation=Rotation)
        bar!(b2_v2,cumI_diffs[2:7],color=colors2,ls=:dash,lw=0,bar_widths=.5,
                legend=false,size=(width1,length1),xrotation=Rotation)
    b3=bar(cumIs,color=colors2,title="Total number of cases in "*wa_name,
            legend=false,xticks=([0:1:size(xticks_labels,1);], xticks_labels),size=(width2,length1),xrotation=Rotation)
    b5=bar(peak_diffs,color=colors2,title="Peak time difference in "*wa_name,
            legend=false,xticks=([0:1:size(xticks_labels,1);], xticks_labels),size=(width2,length1),xrotation=Rotation)
    b6=bar(peaks,color=colors2,title="Peak times in "*wa_name,
            legend=false,xticks=([0:1:size(xticks_labels,1);], xticks_labels),size=(width2,length1),xrotation=Rotation)
    b8=bar(cumCs[2:7],color=colors2,ls=:solid,title="Total number of contacted in "*wa_name*" - R0="*string(sessionParams.R₀),fillalpha=.8,
            legend=false,xticks=([0:1:size(cumI_diffs[2:7],1);], [xticks_labels[1];xticks_labels[3:end]]),size=(width1,length1),xrotation=Rotation)
        bar!(b8,cumCs[9:end],color=colors2,ls=:dash,lw=0,bar_widths=.5,
                legend=false,size=(width1,length1),xrotation=Rotation);
    b9=bar(cumICs[2:7],color=colors2,ls=:solid,title="Contacted infecteds (E or Ix) in "*wa_name*" - R0="*string(sessionParams.R₀),fillalpha=.8,
            legend=false,xticks=([0:1:size(cumI_diffs[2:7],1);], [xticks_labels[1];xticks_labels[3:end]]),size=(width1,length1),xrotation=Rotation)
        bar!(b9,cumICs[9:end],color=colors2,ls=:dash,lw=0,bar_widths=.5,
                legend=false,size=(width1,length1),xrotation=Rotation);
    b10=bar(cumDeaths[8:end],color=colors2,ls=:dash,title="Total deaths in "*wa_name*" - R0="*string(sessionParams.R₀),fillalpha=.8,
            legend=false,xticks=([0:1:size(cumDeaths[1:7],1);], [xticks_labels[1];xticks_labels[2:end]]),size=(width1,length1),xrotation=Rotation)
        bar!(b10,cumDeaths[1:7],color=colors2,ls=:solid,lw=0,bar_widths=.5,
                legend=false,size=(width1,length1),xrotation=Rotation);

    savefig(b,plot_folder*"jl_cumI_"*wa_name*"_box.png");        savefig(b2,plot_folder*"jl_cumI_gain_"*wa_name*"_bar.png");
                                                                 savefig(b2_v2,plot_folder*"jl_cumI_gain_"*wa_name*"_bar_v2.png")
    savefig(b3,plot_folder*"jl_cumI_"*wa_name*"_bar.png");       savefig(b4,plot_folder*"jl_cumI_gain_"*wa_name*"_box.png")
    savefig(b5,plot_folder*"jl_peak_diffs_"*wa_name*"_bar.png"); savefig(b6,plot_folder*"jl_peaks_"*wa_name*"_bar.png")
    savefig(b7,plot_folder*"jl_peaks_"*wa_name*"_box.png");      savefig(b8,plot_folder*"jl_cpeaks_"*wa_name*"_bar.png")

    savefig(b8,plot_folder*"jl_cumC_"*wa_name*"_bar.png");      savefig(b9,plot_folder*"jl_cumIC_"*wa_name*"_bar.png");
    savefig(b10,plot_folder*"jl_cumDeaths_"*wa_name*"_bar.png");
end

######### INFECTEDS AND CT AREA
function oneCurve_CT_comparedtoNoIntervention(no_int_data_file,data_file,plot_folder)
    if !isdir(plot_folder)   mkdir(plot_folder)   end
    wa=12;wa_name="Kilifi"
    @load data_file  sessionParams
    @load data_file  sims_vector
    sim=10
    if no_int_data_file==""    ##do not add no intervention red curve to the plots
        IKilifi=[sims_vector[sim][5][t][12][1] for t=1:366]
        p=plot(IKilifi,fillalpha=.7,label="Infected",size=(width1*.8,length1*.7),xlims=(0,300),
                title="Example of an epidemic curve in Kilifi with 3 months \ncontact tracing, R0="*string(sessionParams.R₀)*",det="*string(Int(sessionParams.δ*100))*"%,a="*string(Int(sum(sessionParams.τₚ)*100))*"%")
        CT_start=findfirst(x->x!=0,[sims_vector[sim][2][t][12][1] for t=1:366])
        vspan!(p,[CT_start, CT_start+30], linecolor = false, fillcolor = colors2[5], fillalpha=.4,label="1st month CT")
        vspan!(p,[CT_start+31, CT_start+60], linecolor = false, fillcolor = colors2[4], fillalpha=.4,label="2nd month CT")
        vspan!(p,[CT_start+61, CT_start+90], linecolor = false, fillcolor = colors2[1], fillalpha=.4,label="3rd month CT")
        savefig(p,plot_folder*"jl_exampleI_"*wa_name*"_R0="*string(sessionParams.R₀)*",det="*string(Int(sessionParams.δ*100))*",a="*string(Int(sum(sessionParams.τₚ)*100))*".png");
    else                         ##Add no intervention red curve to the plots
        @load no_int_data_file  sessionParams
        @load no_int_data_file  sims_vector
        IKilifi_nointervention=[sims_vector[sim][5][t][12][1] for t=1:366]
        p=plot(IKilifi_nointervention,fillalpha=.7,label="No interv",size=(width1*.8,length1*.7),xlims=(0,300),color=:red,
                title="Example of an epidemic curve in Kilifi with 3 months \ncontact tracing, R0="*string(sessionParams.R₀)*",det="*string(Int(sessionParams.δ*100))*"%,a="*string(Int(sum(sessionParams.τₚ)*100))*"%")
        @load data_file  sessionParams
        @load data_file  sims_vector
        IKilifi=[sims_vector[sim][5][t][12][1] for t=1:366]
        plot!(p,IKilifi,fillalpha=.7,label="Infected",size=(width1*.8,length1*.7),xlims=(0,300),color=:blue)
        CT_start=findfirst(x->x!=0,[sims_vector[sim][2][t][12][1] for t=1:366])
        vspan!(p,[CT_start, CT_start+30], linecolor = false, fillcolor = colors2[5], fillalpha=.4,label="1st month CT")
        vspan!(p,[CT_start+31, CT_start+60], linecolor = false, fillcolor = colors2[4], fillalpha=.4,label="2nd month CT")
        vspan!(p,[CT_start+61, CT_start+90], linecolor = false, fillcolor = colors2[1], fillalpha=.4,label="3rd month CT")
        savefig(p,plot_folder*"jl_exampleI_"*wa_name*"_R0="*string(sessionParams.R₀)*",det="*string(Int(sessionParams.δ*100))*",a="*string(Int(sum(sessionParams.τₚ)*100))*"_v2.png");
    end
end
function CT_comparedtoNoIntervention(no_int_data_file,data_file,plot_folder)
    if !isdir(plot_folder)   mkdir(plot_folder)   end
    wa=12;wa_name="Kilifi"

    @load no_int_data_file  sessionParams
    @load no_int_data_file  sims_vector
    p=plot(size=(width1*.8,length1*.75),xlims=(0,300))

    for sim=1:size(sims_vector,1)
        IKilifi_nointervention=[sims_vector[sim][5][t][12][1] for t=1:366]
        if sim==1       plot!(p,IKilifi_nointervention,linealpha=.3,color=:red,label="No interv")
        else            plot!(p,IKilifi_nointervention,linealpha=.3,xlims=(0,300),color=:red,label=false)   end
    end
    @load data_file  sessionParams
    @load data_file  sims_vector
    plot!(p,title="Stochastic epidemic curves in Kilifi with 3 months \ncontact tracing, R0="*string(sessionParams.R₀)*",sympt="*string(Int(sessionParams.δ*100))*"%,det="*string(Int(sum(sessionParams.τₚ)*100))*"%")
    for sim=1:size(sims_vector,1)
        IKilifi=[sims_vector[sim][5][t][12][1] for t=1:366]
        if sim==1        plot!(p,IKilifi,linealpha=.3,size=(width1*.8,length1*.7),xlims=(0,300),color=:blue,label="Intervention")
        else            plot!(p,IKilifi,linealpha=.3,size=(width1*.8,length1*.7),xlims=(0,300),color=:blue,label=false); end
        if sim==1
            CT_start=findfirst(x->x!=0,[sims_vector[sim][2][t][12][1] for t=1:366])
            vspan!(p,[CT_start, CT_start+30], linecolor = false, fillcolor = colors2[5], fillalpha=.4,label="1st month CT")
            vspan!(p,[CT_start+31, CT_start+60], linecolor = false, fillcolor = colors2[4], fillalpha=.4,label="2nd month CT")
            vspan!(p,[CT_start+61, CT_start+90], linecolor = false, fillcolor = colors2[1], fillalpha=.4,label="3rd month CT")
        end
    end
    #plot!(legend=false)
    savefig(p,plot_folder*"jl_I_"*wa_name*"_R0="*string(sessionParams.R₀)*",det="*string(Int(sessionParams.δ*100))*",a="*string(Int(sum(sessionParams.τₚ)*100))*"_v2.png");
end

########
S0=[4.138758e6, 867417.0, 2.326182e6, 8.084069e6, 3.229145e6, 459761.0, 999280.0, 1.979082e6, 926952.0, 340661.0, 2.381706e6, 2.126254e6, 2.960717e6, 786461.0, 7.478259e6, 781212.0, 2.114588e6, 4.094022e6, 569586.0, 917976.0]
S012A=[279992, 264585, 246992, 206735, 220761, 211260, 181931, 130091, 111761, 82036, 55221, 42469, 34381, 23781, 16214, 18044]
#folder="./contacts/results_session1-2-3-4_50sims_0_COPY/"
width1=400;length1=400;
width2=1000;length2=800;
#cumI_medians_kenya,cumI_stds_kenya,cumI_medians,cumI_stds,peaks,peak_stds=make_plots_v7(folder,S0)

####### Kenya and Kilifi-Ages with No Intervention
data_files=["./contacts/results_session20_50sims/sims_sc2000.jld2","./contacts/results_session22_50sims/sims_sc2200.jld2"]
plot_folder="./contacts/results_sessions20-22_det0_50sims/"
#cumI_medians_kenya,cumDeaths_stds_kenya,cumI_medians,cumDeaths_medians=kenya_barplot(data_files,plot_folder,S0)
#=cumI_medians_kenya[2]/sum(S0)
cumDeaths_stds_kenya[1]
cumI_medians[2][12]/S0[12]=#

#cumIsA=kilifi_ages(data_files,plot_folder,S0)
#medians=[median(cumIsA[1][a])   for a=1:16];sum(medians[13:end])/sum(S012A[13:end]);ages[13]

###### Ages with intervention
data_files=["./contacts/results_session20_50sims/sims_sc2000.jld2","./contacts/results_session20_50sims/sims_sc2053.jld2","./contacts/results_session21_50sims/sims_sc2153.jld2"]
plot_folder="./contacts/results_sessions20to21_det5_90days_50sims/"
#cumIfinalKilifiA=kilifi_ages_α5_90days(data_files,plot_folder,S0)

###### CumIs 2δs
Rotation=40
folders=["./contacts/results_session22_50sims/","./contacts/results_session23_50sims/"];
plot_folder="./contacts/results_session22_50sims/plots_sessions22-23/"
δs=[.3,.7]
#make_plots_twoδs(folders,plot_folder,S0,δs)

#######
data_files=["./contacts/results_session20_50sims/sims_sc2000.jld2","./contacts/results_session20_50sims/sims_sc2053.jld2","./contacts/results_session20_50sims/sims_sc2093.jld2","./contacts/results_session21_50sims/sims_sc2153.jld2","./contacts/results_session21_50sims/sims_sc2193.jld2"]
plot_folder="./contacts/results_sessions20to21_det5_90days_50sims/"

#kilifi_ages_α5_9_90days(data_files,plot_folder,S0)
#kilifi_ages_α5_9_90days_over60s(data_files,plot_folder,S0)


########
plot_folder="./contacts/results_session20_50sims/oneExample20-23/"
#=oneCurve_CT_comparedtoNoIntervention("./contacts/results_session20_50sims/sims_sc2000.jld2","./contacts/results_session20_50sims/sims_sc2093.jld2",plot_folder)
oneCurve_CT_comparedtoNoIntervention("./contacts/results_session21_50sims/sims_sc2100.jld2","./contacts/results_session21_50sims/sims_sc2193.jld2",plot_folder)
oneCurve_CT_comparedtoNoIntervention("./contacts/results_session22_50sims/sims_sc2200.jld2","./contacts/results_session22_50sims/sims_sc2293.jld2",plot_folder)
oneCurve_CT_comparedtoNoIntervention("./contacts/results_session23_50sims/sims_sc2300.jld2","./contacts/results_session23_50sims/sims_sc2393.jld2",plot_folder)
=#
plot_folder="./contacts/results_session20_50sims/Icurves20-23/"
CT_comparedtoNoIntervention("./contacts/results_session20_50sims/sims_sc2000.jld2","./contacts/results_session20_50sims/sims_sc2093.jld2",plot_folder)
CT_comparedtoNoIntervention("./contacts/results_session21_50sims/sims_sc2100.jld2","./contacts/results_session21_50sims/sims_sc2193.jld2",plot_folder)
CT_comparedtoNoIntervention("./contacts/results_session22_50sims/sims_sc2200.jld2","./contacts/results_session22_50sims/sims_sc2293.jld2",plot_folder)
CT_comparedtoNoIntervention("./contacts/results_session23_50sims/sims_sc2300.jld2","./contacts/results_session23_50sims/sims_sc2393.jld2",plot_folder)
