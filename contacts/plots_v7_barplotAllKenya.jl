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
ages=[string(i)*"-"*string(i+4) for i=0:5:70];push!(ages,"75+yrs")#"0-4" "5-9" "10-14" "15-19" "20-24" "25-29" "30-34" "35-39" "40-44" "45-49" ]
########
function kenya_barplot(data_files,plot_folder,S0)
    if !isdir(plot_folder)   mkdir(plot_folder)   end
    wa=12;wa_name="Kilifi"

    final_cum0detection=0;    cumIs=[];   cumI_medians=[];    cumI_medians_kenya=[];
    for i=1:size(data_files,1)
        @load data_files[i]  sims_vector
        @load data_files[i]  sessionParams
        final_cum=[];
        for sim=1:size(sims_vector,1) #for each sim
            push!(final_cum,sum(sims_vector[sim][9][wa,:,1:3])) # final cumulative I=IA+ID+IQ, summed for all ages
        end
        push!(cumIs,median(final_cum));

        ##final_cumI_barplot:
        cumI=[[sum(sims_vector[sim][1][end][wa,1:3]) for wa=1:20]   for sim=1:size(sims_vector,1)]
        median_cumI=[median([cumI[sim][wa] for sim=1:size(sims_vector,1)])      for wa=1:20]
        push!(cumI_medians,median_cumI)

        ##All kenya cumI
        median_cumI_kenya=median([sum([cumI[sim][wa]  for wa=1:20])  for sim=1:size(sims_vector,1)])
        push!(cumI_medians_kenya,median_cumI_kenya)

    end
    S02=[]; for wa=1:20 push!(S02,S0[wa]);   end
    barplotALL=bar([1.5:2:40.5;],S02,color=[:green,:chartreuse4],xticks=([1:2:40;],riskregionnames),fillalpha=.5,linewidth=false,bar_width=2,linecolor=false,
                    size=(width1*1.2,length1*.8),rotation=40,label=false,title="Final number of cases in Kenya without intervention")
    cumI_medians2=[];
    for wa=1:20,i=1:1size(data_files,1)    push!(cumI_medians2,cumI_medians[i][wa]);   end

    bar!(barplotALL,cumI_medians2,color=[colors2[1],colors2[6]],linecolor=false,label=false,bar_width=1.005)#,yerrorbar=cumI_allfiles_sd,markercolor=:black)
    savefig(barplotALL,plot_folder*"jl_cumI_barplot.png")
    return cumI_medians_kenya,cumI_medians

end

function kilifi_ages(data_files,plot_folder,S0)
    if !isdir(plot_folder)   mkdir(plot_folder)   end
    wa=12;wa_name="Kilifi"

    @load data_files[1]  sims_vector
    cumIsA=[[] for i=1:16]
    for sim=1:size(sims_vector,1),a=1:16    push!(cumIsA[a],sum(sims_vector[sim][9][wa,a,1:3]));end

    @load data_files[2]  sims_vector# 9 is sc53
    cumIsA_2=[[] for i=1:16]
    for sim=1:size(sims_vector,1),a=1:16    push!(cumIsA_2[a],sum(sims_vector[sim][9][wa,a,1:3])); end

    S012A_2=[]
    for i=1:16  push!(S012A_2,S012A[i]);push!(S012A_2,S012A[i]) ;end
    b=bar(S012A_2,bar_width=1.01,color=[:green,:chartreuse4],xticks=([1:2:32;],ages),fillalpha=.5,linecolor=false,rotation=30,label=false)#bar_edges=false)

    cumIfinal=[];
    for a=1:16  push!(cumIfinal,median(cumIsA[a]));push!(cumIfinal,median(cumIsA_2[a])); end
    bar!(b,cumIfinal,color=color=[colors2[1],colors2[6]],rotation=30,contour_labels=true,size=(width1*1.2,length1*.8),
            linecolor=false,label=false,title="Estimated cases per age in Kilifi without intervention")
    savefig(b,plot_folder*"jl_cumI_ages_kilifi.png")
    return cumIsA
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

function kilifi_ages_α5_90days(data_files,plot_folder,S0)
    if !isdir(plot_folder)   mkdir(plot_folder)   end
    wa=12;wa_name="Kilifi"

    @load data_files[1]  sims_vector
    cumIsA=[[] for i=1:16];    for sim=1:size(sims_vector,1),a=1:16    push!(cumIsA[a],sum(sims_vector[sim][9][wa,a,1:3]));end

    @load data_files[2]  sims_vector#
    cumIsA_2=[[] for i=1:16];    for sim=1:size(sims_vector,1),a=1:16    push!(cumIsA_2[a],sum(sims_vector[sim][9][wa,a,1:3])); end

    @load data_files[3]  sims_vector#
    cumIsA_3=[[] for i=1:16];    for sim=1:size(sims_vector,1),a=1:16    push!(cumIsA_3[a],sum(sims_vector[sim][9][wa,a,1:3])); end

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
    cumIsA=[[] for i=1:16];    for sim=1:size(sims_vector,1),a=1:16    push!(cumIsA[a],sum(sims_vector[sim][9][wa,a,1:3]));end

    @load data_files[2]  sims_vector#
    cumIsA_2=[[] for i=1:16];    for sim=1:size(sims_vector,1),a=1:16    push!(cumIsA_2[a],sum(sims_vector[sim][9][wa,a,1:3])); end

    @load data_files[3]  sims_vector#
    cumIsA_3=[[] for i=1:16];    for sim=1:size(sims_vector,1),a=1:16    push!(cumIsA_3[a],sum(sims_vector[sim][9][wa,a,1:3])); end

    @load data_files[4]  sims_vector#
    cumIsA_4=[[] for i=1:16];    for sim=1:size(sims_vector,1),a=1:16    push!(cumIsA_4[a],sum(sims_vector[sim][9][wa,a,1:3])); end

    @load data_files[5]  sims_vector#
    cumIsA_5=[[] for i=1:16];    for sim=1:size(sims_vector,1),a=1:16    push!(cumIsA_5[a],sum(sims_vector[sim][9][wa,a,1:3])); end


    S012A_2=[]
    for i=1:16  push!(S012A_2,S012A[i]);push!(S012A_2,S012A[i]) ;push!(S012A_2,S012A[i]) ;end
    b=bar(S012A_2,bar_width=1.01,color=[:green,:green,:green,:chartreuse4,:chartreuse4,:chartreuse4],xticks=([1:3:48;],ages),fillalpha=.5,linecolor=false,rotation=30,label=false)#bar_edges=false)

    cumIfinal=[];
    for a=1:16  push!(cumIfinal,median(cumIsA[a]));push!(cumIfinal,median(cumIsA_2[a]));push!(cumIfinal,median(cumIsA_3[a])); end
    bar!(b,cumIfinal,color=color=[colors2[1],colors2[6],colors2[4]],rotation=30,contour_labels=true,size=(width1*1.2,length1*.8),
            linecolor=false,label=false,title="Estimated cases per age in Kilifi R0=2.5,a=50%,90days")

            cumIfinal4=[];
            for a=1:16  push!(cumIfinal4,median(cumIsA_4[a])); end
    bar!(b,[2:3:48;],cumIfinal4,color=colors2[6],rotation=30,bar_width=.5,lw=.5,#contour_labels=true,size=(width1*1.2,length1*.8),
            ls=:dash,label=false)#,title="Estimated cases per age in Kilifi without intervention")

            cumIfinal5=[];
            for a=1:16  push!(cumIfinal5,median(cumIsA_5[a])); end
    bar!(b,[3:3:48;],cumIfinal5,color=colors2[4],rotation=30,bar_width=.5,lw=.5,ls=:dash,label=false)
    display(b)
    savefig(b,plot_folder*"jl_cumI_ages_kilifi_2.png")
    return cumIsA
end
########
S0=[4.138758e6, 867417.0, 2.326182e6, 8.084069e6, 3.229145e6, 459761.0, 999280.0, 1.979082e6, 926952.0, 340661.0, 2.381706e6, 2.126254e6, 2.960717e6, 786461.0, 7.478259e6, 781212.0, 2.114588e6, 4.094022e6, 569586.0, 917976.0]
S012A=[279992, 264585, 246992, 206735, 220761, 211260, 181931, 130091, 111761, 82036, 55221, 42469, 34381, 23781, 16214, 18044]
#folder="./contacts/results_session1-2-3-4_50sims_0_COPY/"
width1=400;length1=400;
width2=1000;length2=800;
#cumI_medians_kenya,cumI_stds_kenya,cumI_medians,cumI_stds,peaks,peak_stds=make_plots_v7(folder,S0)
data_files=["./contacts/results_session10_20sims/sims_sc1000.jld2","./contacts/results_session11_20sims/sims_sc1153.jld2"]
data_files=["./contacts/results_session10_20sims/sims_sc1000.jld2","./contacts/results_session10_20sims/sims_sc1053.jld2","./contacts/results_session11_20sims/sims_sc1153.jld2"]
plot_folder="./contacts/results_sessions10to11_det5_90days_20sims/"
#cumI_medians_kenya,cumI_medians=kenya_barplot(data_files,plot_folder,S0)
cumIfinalKilifiA=kilifi_ages_α5_90days(data_files,plot_folder,S0)

data_files=["./contacts/results_session10_20sims/sims_sc1000.jld2","./contacts/results_session10_20sims/sims_sc1053.jld2","./contacts/results_session10_20sims/sims_sc1093.jld2","./contacts/results_session11_20sims/sims_sc1153.jld2","./contacts/results_session11_20sims/sims_sc1193.jld2"]
plot_folder="./contacts/results_sessions10to11_det5_90days_20sims/"
cumIfinalKilifiA=kilifi_ages_α5_9_90days(data_files,plot_folder,S0)


#sum(cumIfinalKilifiA[13][12]+cumIfinalKilifiA[14][12]+cumIfinalKilifiA[15][12]+cumIfinalKilifiA[16][12])/sum(S012A[13:16])
#cumI_medians[1][12]/S0[12]
#cumI_medians_kenya[1]/sum(S0)
#make_plots_ages(folder,S0)

#cumIKenya=[2.2839415e6, 302660.0, 1.1054025e6, 4.437384e6, 1.455775e6, 191055.0, 411214.0, 946572.0, 391486.0, 178220.0, 1.1861265e6, 1.0601145e6, 1.3771935e6, 356295.0, 3.38641e6, 299167.0, 1.0465525e6, 1.783778e6, 251199.0, 414011.0]
#cumIKenya_percentage=[cumIKenya[wa]*100/S0[wa]  for wa=1:20]

#sum(cumIKenya)/sum(S0)
#cumIKenya[12]
#cumIKenya_sd[12]/cumIKenya[12]


#@load "./contacts\\results_session12_20sims\\sims_sc1293.jld2" sims_vector
