using Plots,MAT, StatsPlots, Statistics,JLD2, Images, ImageView,ImageDraw

colors=[:blue, :orange, :purple3, :maroon, :gold, :orangered, :grey, :purple]
markershapes=[:cross,:star4, :vcross, :star6, :hline]
riskregionnames = ["Machakos/Muranga" "Mandera" "Baringo/Nakuru" "Nairobi" "Uasin Gishu" "Marsabit" "Garissa" "Nakuru/Narok" "Turkana" "Taita Taveta" "Kitui/Meru" "Kilifi/Mombasa" "Kericho/Kisumu" "Kilifi/Lamu" "Kakamega/Kisumu" "Wajir" "Kajiado/Kisumu" "Homa bay/Migori" "Samburu/Laikipia" "Kilifi/Kwale" "Total"]
wa_coords=[[300,450], [515,85], [165,360], [235,465], [115,300], [300,140], [510,380], [180, 430], [120, 130], [355, 615], [340, 375], [465, 630], [100, 380], [495, 530], [40, 355], [490, 255], [155, 495], [30, 440], [250, 290], [400, 670]]

########
#=
function one_sim_wa(τₚ_list,results_folder,data_files,n_traj,wa_coords,wa,wa_name)
    p1=plot(title="1sim example I=A+D+IQ "*wa_name,legend=:topright);
    p2=plot(title="1sim example Cumulative I=A+D+IQ "*wa_name,legend=:topright);
    p3=plot(title="1sim example Q "*wa_name,legend=:topright);
    ##p4=plot(title="1sim example S, R Kilifi",legend=:topleft);
    for i=1:size(τₚ_list,1)
        @load results_folder*data_files[i]  sims_vector
        I=[sims_vector[1][3][t][wa][1] for t=1:size(sims_vector[1][2],1)]
        plot!(p1,I,label="detection"*string(τₚ_list[i]*100)*"%-I", color=colors[i])
        cumIA=[sims_vector[1][1][t][wa,1] for t=1:size(sims_vector[1][1],1)]
        cumID=[sims_vector[1][1][t][wa,2] for t=1:size(sims_vector[1][1],1)]
        cumIQ=[sims_vector[1][1][t][wa,3] for t=1:size(sims_vector[1][1],1)]
        cumI=cumIA .+ cumID .+cumIQ #OR I=[sum(sims_vector[1][1][t][wa,1:3]) for t=1:size(sims_vector[1][1],1)]
        plot!(p2,cumI,label="detection"*string(τₚ_list[i]*100)*"%-cumI", color=colors[i])
        plot!(p2,cumIA,label="detection"*string(τₚ_list[i]*100)*"%-IA", color=colors[i], line=:dash)
        plot!(p2,cumID,label="detection"*string(τₚ_list[i]*100)*"%-ID", color=colors[i], line=:dot)
        plot!(p2,cumIQ,label="detection"*string(τₚ_list[i]*100)*"%-IQ", color=colors[i], line=(2,:dashdotdot))
        Q=[sims_vector[1][2][t][wa][1] for t=1:size(sims_vector[1][2],1)]
        plot!(p3,Q,label="detection"*string(τₚ_list[i]*100)*"%-Q", color=colors[i])

        ##TO DELETE WHEN RUNNING LARGE SIMS:
        #=S=[sims_vector[1][5][t][wa][1] for t=1:size(sims_vector[1][1],1)]
        plot!(p4,S,label="detection"*string(τₚ_list[i]*100)*"%-S", color=colors[i])
        R=[sims_vector[1][6][t][wa][1] for t=1:size(sims_vector[1][1],1)]
        plot!(p4,R,label="detection"*string(τₚ_list[i]*100)*"%-R", color=colors[i])=#
    end
    savefig(p1,results_folder*"jl_ONE_SIM_I_"*wa_name*".png")
    savefig(p2,results_folder*"jl_ONE_SIM_CumI_"*wa_name*".png")
    savefig(p3,results_folder*"jl_ONE_SIM_Q_"*wa_name*".png")
    ##savefig(p4,results_folder*"jl_ONE_SIM_SR_kilifi.png")
end
function time_intro2peak_wa(τₚ_list,results_folder,data_files,n_traj,wa_coords,wa,wa_name)
    b=boxplot(title="Time from introduction into "*wa_name*" to peak (in days)")
    for i=1:size(τₚ_list,1)
        @load results_folder*data_files[i]  sims_vector

        ##Time from intro to peak
        intro2peak=[]
        for sim=1:size(sims_vector,1) #for each sim
            ##Calculating for each sim, the time of introduction:
            Q=[sims_vector[sim][2][t][wa][1] for t=1:size(sims_vector[1][2],1)]
            indices=findall(x -> x>0 , Q)
            tintro=0
            if size(indices,1)>0
                tintro=indices[1]
            end
            ##calculating for each sim, the time of peak
            I=[sims_vector[sim][2][t][wa][1] for t=1:size(sims_vector[1][2],1)]
            tpeak=findall(x->x==maximum(I),I)[1]
            push!(intro2peak,tpeak-tintro)
        end
        boxplot!(b,intro2peak,color=colors[i],label="detection"*string(τₚ_list[i]*100)*"%",legend=:bottomright)
    end
    savefig(b,results_folder*"jl_intro2peak_"*wa_name*".png")
end
function final_cumI_wa(τₚ_list,results_folder,data_files,n_traj,wa_coords,wa,wa_name)
    b=boxplot(title="Final total cumulative Infecteds in "*wa_name*" (all ages)")
    final_cum0detection=0
    cumI_diff=[]
    cumIs=[]
    for i=1:size(τₚ_list,1)
        @load results_folder*data_files[i]  sims_vector
        final_cum=[]
        for sim=1:size(sims_vector,1) #for each sim
            push!(final_cum,sum(sims_vector[sim][4][wa,:,1:3])) # final cumulative I=IA+ID+IQ, summed for all ages
        end
        push!(cumIs,median(final_cum))
        boxplot!(b,final_cum,color=colors[i],label="detection"*string(τₚ_list[i]*100)*"%",legend=:bottomright)
        if i==1     final_cum0detection=median(final_cum)       end
        push!(cumI_diff,final_cum0detection-median(final_cum))
    end
    println("cumI_diff=",cumI_diff,"   cumIs=",cumIs)
    b2=bar(cumI_diff,color=colors,legend=false,title="Gain in total Infecteds (cum) in "*wa_name*" (all ages)",
            xticks = ([1:1:size(τₚ_list,1);], [string(e*100)*"%" for e in τₚ_list]))
    b3=bar(cumIs,color=colors,legend=false,title="Total Infecteds (cum) in "*wa_name*" (all ages)",
            xticks = ([1:1:size(τₚ_list,1);], [string(e*100)*"%" for e in τₚ_list]))
    #display(b)
    savefig(b,results_folder*"jl_cumI_"*wa_name*"_box.png")
    savefig(b2,results_folder*"jl_cumI_gain_"*wa_name*"_bar.png")
    savefig(b3,results_folder*"jl_cumI_"*wa_name*"_bar.png")
end
function final_cumI_map(τₚ_list,results_folder,data_files,n_traj,wa_coords,S0)
    imgs=[load("./contacts/kenya.jpg") for τₚ in τₚ_list]
    CList=reshape(range(colorant"green",stop=colorant"red",length=20),1,20)
    for i=1:size(τₚ_list,1)
        @load results_folder*data_files[i]  sims_vector
        cumI=[[sum(sims_vector[sim][1][end][wa,1:3]) for wa=1:20]   for sim=1:size(sims_vector,1)]
        median_cumI=[median([cumI[sim][wa] for sim=1:size(sims_vector,1)])      for wa=1:20]
        for wa=1:20
            draw!(imgs[i], CirclePointRadius(wa_coords[wa][1],wa_coords[wa][2],35),CList[Int64(ceil(median_cumI[wa]*20/S0[wa]))])
        end
    end
    p=hcat(imgs[1],imgs[2],imgs[3],imgs[4],imgs[5])
    #display(p)
    save(results_folder*"jl_cumI_MAP.png",p)
end
function final_cumI_barplot(τₚ_list,results_folder,data_files,n_traj,wa_coords,S0)
    imgs=[load("./contacts/kenya.jpg") for τₚ in τₚ_list]

    S02=[]
    for wa=1:20,i=1:size(τₚ_list,1)
        push!(S02,S0[wa])
    end
    barplotALL=bar(S02,color=:green,xticks = ([1:size(τₚ_list,1):20*size(τₚ_list,1);], riskregionnames),rotation = 60,legend=false)
    CList=reshape(range(colorant"green",stop=colorant"red",length=100),1,100)
    medians_cumI=[]
    for i=1:size(τₚ_list,1)
        @load results_folder*data_files[i]  sims_vector
        cumI=[[sum(sims_vector[sim][1][end][wa,1:3]) for wa=1:20]   for sim=1:size(sims_vector,1)]
        median_cumI=[median([cumI[sim][wa] for sim=1:size(sims_vector,1)])      for wa=1:20]
        push!(medians_cumI,median_cumI)
        for wa=1:20
            draw!(imgs[i], CirclePointRadius(wa_coords[wa][1],wa_coords[wa][2],35),CList[Int64(ceil(median_cumI[wa]*100/S0[wa]))])
        end
    end
    p=hcat(imgs[1],imgs[2],imgs[3],imgs[4],imgs[5])
    #display(p)
    save(results_folder*"jl_cumI_map.png",p)

    medians_cumI2=[]
    colors2=[]
    for wa=1:20,i=1:size(τₚ_list,1)
        push!(medians_cumI2,medians_cumI[i][wa])
        push!(colors2,colors[i])
    end
    bar!(barplotALL,medians_cumI2,color=colors2)
    savefig(barplotALL,results_folder*"jl_cumI_barplot.png")
end=#

########
function get_cumIgain_wa(τₚ_list,results_folder,data_files,n_traj,wa_coords,wa,wa_name)
    final_cum0detection=0
    cumI_diff=[]
    for i=1:size(τₚ_list,1)
        @load results_folder*data_files[i]  sims_vector
        final_cum=[]
        for sim=1:size(sims_vector,1) #for each sim
            push!(final_cum,sum(sims_vector[sim][4][wa,:,1:3])) # final cumulative I=IA+ID+IQ, summed for all ages
        end
        if i==1     final_cum0detection=median(final_cum)       end
        push!(cumI_diff,final_cum0detection-median(final_cum))
    end
    return cumI_diff
end
function cumIgain_many_sessions(cumI_diffs,Κ_max_capacity_Kilifi_string,sessions,τₚ_list,wa,wa_name)
    bar_vector=[]
    colors2=[]
    ticks2=["detection&capacity"]
    for j=2:size(τₚ_list,1),i=1:size(sessions,1)
        push!(bar_vector,cumI_diffs[i][j])
        push!(colors2,colors[j])
        push!(ticks2,string(Int(floor(τₚ_list[j]*100)))*"% & "*Κ_max_capacity_Kilifi_string[i])
    end
    b=bar(bar_vector,color=colors2,legend=false,title="Gain in final cumulative Infecteds in Kilifi (all ages)",
            xticks = ([0:1:size(bar_vector,1);], ticks2),rotation = 60)
    #display(b)
    for session in sessions
        results_folder=".\\contacts\\results_session"*string(Int(floor(session/10)))*"0s\\results_session"*string(session)*"\\"
        savefig(b,results_folder*"jl_cumIgain_scenarios.png")
    end
end
function get_peaks_wa(τₚ_list,results_folder,data_files,n_traj,wa_coords,wa,wa_name)
    #final_cum0detection=0
    #cumI_diff=[]
    peaks=[]
    peak_diffs=[]
    median_peak0=0
    for i=1:size(τₚ_list,1)
        @load results_folder*data_files[i]  sims_vector
        #final_cum=[]
        peak=[]
        for sim=1:size(sims_vector,1) #for each sim
            I=[sims_vector[sim][3][t][wa][1] for t=1:size(sims_vector[sim][3],1)]
            tpeak=findall(x->x==maximum(I),I)[1]
            push!(peak,tpeak)
        end
        push!(peaks,median(peak))
        if i==1     median_peak0=median(peaks)       end
        push!(peak_diffs,median(peaks)-median_peak0)
    end
    return peaks,peak_diffs
end
function peaks_many_sessions(peakss,peak_diffss,Κ_max_capacity_Kilifi_string,sessions,τₚ_list,wa,wa_name)
    bar_vector=[];    bar_vectorDiffs=[]
    colors2=[]
    ticks2=["detection&capacity"]
    for j=1:size(τₚ_list,1),i=1:size(sessions,1)
        push!(bar_vector,peakss[i][j]);push!(bar_vectorDiffs,peak_diffss[i][j])
        push!(colors2,colors[j])
        push!(ticks2,string(Int(floor(τₚ_list[j]*100)))*"% & "*Κ_max_capacity_Kilifi_string[i])
    end
    b=bar(bar_vector,color=colors2,legend=false,title="Peaks in Kilifi (all ages)",
            xticks = ([0:1:size(bar_vector,1);], ticks2),rotation = 60)
    #display(b)
    b2=bar(bar_vectorDiffs,color=colors2,legend=false,title="Peak gain in Kilifi (all ages)",
            xticks = ([0:1:size(bar_vectorDiffs,1);], ticks2),rotation = 60)
    for session in sessions
        results_folder=".\\contacts\\results_session"*string(Int(floor(session/10)))*"0s\\results_session"*string(session)*"\\"
        savefig(b,results_folder*"jl_peaks_scenarios.png")
        savefig(b2,results_folder*"jl_peakdiffs_scenarios.png")
    end
end
########
function make_plots_v7(folder,S0)
    data_files=readdir(folder)
    data_files=[folder*s for s in data_files if endswith(s,".jld2")&&s!="stats.jld2"]

    cumI_diffs=[];peakss=[];peak_diffss=[];
    wa=12;wa_name="Kilifi"

    b=boxplot(title="Final total cumulative Infecteds in "*wa_name*" (all ages)")
    final_cum0detection=0;        cumI_diffs=[];        cumIs=[];    medians_cumI=[]
    xticks_labels=[];colors2=[]

    for i=1:size(data_files,1)#for i=1:size(τₚ_list,1)
        @load data_files[i]  sims_vector
        @load data_files[i]  sessionParams
        push!(xticks_labels,"sc"*string(sessionParams.sc_nb))
        push!(colors2,rand(colors))
        ##final_cumI_wa:
        final_cum=[]
        for sim=1:size(sims_vector,1) #for each sim
            push!(final_cum,sum(sims_vector[sim][4][wa,:,1:3])) # final cumulative I=IA+ID+IQ, summed for all ages
        end
        push!(cumIs,median(final_cum))
        #boxplot!(b,final_cum,color=colors2[i],label="detection"*string(sessionParams.τₚ*100)*"%",legend=:bottomright)
        if i==1     final_cum0detection=median(final_cum)       end
        push!(cumI_diffs,final_cum0detection-median(final_cum))

        ##final_cumI_barplot:
        cumI=[[sum(sims_vector[sim][1][end][wa,1:3]) for wa=1:20]   for sim=1:size(sims_vector,1)]
        median_cumI=[median([cumI[sim][wa] for sim=1:size(sims_vector,1)])      for wa=1:20]
        push!(medians_cumI,median_cumI)
    end
    ##final_cumI_wa:
    println("cumI_diff=",cumI_diffs,"   cumIs=",cumIs)
    @save folder*"stats.jld2"   cumI_diffs,cumIs

    b2=bar(cumI_diffs,color=colors2,legend=false,title="Gain in total Infecteds (cum) in "*wa_name*" (all ages)",
            xticks = ([1:1:size(xticks_labels,1);], xticks_labels))
    b3=bar(cumIs,color=colors2,legend=false,title="Total Infecteds (cum) in "*wa_name*" (all ages)",
            xticks = ([1:1:size(xticks_labels,1);], xticks_labels))
    savefig(b,folder*"jl_cumI_"*wa_name*"_box.png")
    savefig(b2,folder*"jl_cumI_gain_"*wa_name*"_bar.png")
    savefig(b3,folder*"jl_cumI_"*wa_name*"_bar.png")
    #=##final_cumI_barplot:
    S02=[]; for wa=1:20,i=1:size(τₚ_list,1) push!(S02,S0[wa]);   end
    barplotALL=bar(S02,color=:green,xticks = ([1:size(τₚ_list,1):20*size(τₚ_list,1);], riskregionnames),rotation = 60,legend=false)

    medians_cumI2=[];
    for wa=1:20,i=1:size(data_files,1)    push!(medians_cumI2,medians_cumI[i][wa]);   end
    bar!(barplotALL,medians_cumI2,color=colors2)
    savefig(barplotALL,folder*"jl_cumI_barplot.png")=#

        ##FOR cumIgain_many_sessions and peaks_many_sessions
        #push!(cumI_diffs,get_cumIgain_wa(τₚ_list,results_folder,data_files,n_traj,wa_coords,12,"Kilifi"))
        #peaks,peak_diffs=get_peaks_wa(τₚ_list,results_folder,data_files,n_traj,wa_coords,12,"Kilifi")
        #push!(peakss,peaks);push!(peak_diffss,peak_diffs)

    #cumIgain_many_sessions(cumI_diffs,["500","1e3","5e3","1e4"],sessions,τₚ_list,12,"Kilifi")
    #peaks_many_sessions(peakss,peak_diffss,["500","1e3","5e3","1e4"],sessions,τₚ_list,12,"Kilifi")
    #return peakss,peak_diffss
end
function make_plots_oneExample(sessions,τₚ_list,S0,Κ_max_capacity_Kilifi)
    wa=12;wa_name="Kilifi"
    p1=plot(title="1sim example I=A+D+IQ "*wa_name,legend=:topright);
    p2=plot(title="1sim example Cumulative I=A+D+IQ "*wa_name,legend=:topright);
    p2Q=plot(title="1sim example Cumulative IQ "*wa_name,legend=:topright);
    p3=plot(title="1sim example Q "*wa_name,legend=:topright);
    p4=plot(title="1sim example S, R Kilifi",legend=:topleft);
    p5=plot(title="1sim example infecteds quarantined after contact Kilifi",legend=:topleft);
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
            I=[sims_vector[1][3][t][wa][1] for t=1:size(sims_vector[1][2],1)]
            plot!(p1,I,label="session"*string(session)*"-"*string(τₚ_list[i]*100)*"%-I", color=colors[j])
            cumIA=[sims_vector[1][1][t][wa,1] for t=1:size(sims_vector[1][1],1)]
            cumID=[sims_vector[1][1][t][wa,2] for t=1:size(sims_vector[1][1],1)]
            cumIQ=[sims_vector[1][1][t][wa,3] for t=1:size(sims_vector[1][1],1)]
            cumI=cumIA .+ cumID .+cumIQ #OR I=[sum(sims_vector[1][1][t][wa,1:3]) for t=1:size(sims_vector[1][1],1)]
            plot!(p2,cumI),label="session"*string(session)*"-"*string(τₚ_list[i]*100)*"%-cumI", color=colors[j])
            plot!(p2,cumIA,label="session"*string(session)*"-"*string(τₚ_list[i]*100)*"%-IA", color=colors[j], line=:dash)
            plot!(p2,cumID,label="session"*string(session)*"-"*string(τₚ_list[i]*100)*"%-ID", color=colors[j], line=:dot)
            plot!(p2,cumIQ,label="session"*string(session)*"-"*string(τₚ_list[i]*100)*"%-IQ", color=colors[j], line=(2,:dashdotdot))
            plot!(p2Q,cumIQ,label="session"*string(session)*"-"*string(τₚ_list[i]*100)*"%-IQ", color=colors[j], line=(2,:dashdotdot))
            Q=[sims_vector[1][2][t][wa][1] for t=1:size(sims_vector[1][2],1)]
            plot!(p3,Q,label="detection"*string(τₚ_list[i]*100)*"%-Q", color=colors[j])

            ##TO DELETE WHEN RUNNING LARGE SIMS:
            S=[sims_vector[1][5][t][wa][1] for t=1:size(sims_vector[1][1],1)]
            plot!(p4,S,label="detection"*string(τₚ_list[i]*100)*"%-S", color=colors[j])
            R=[sims_vector[1][6][t][wa][1] for t=1:size(sims_vector[1][1],1)]
            plot!(p4,R,label="detection"*string(τₚ_list[i]*100)*"%-R", color=colors[j])
            cumIC=[sims_vector[1][7][t][wa][1] for t=1:size(sims_vector[1][1],1)]
            plot!(p5,cumIC,label="detection"*string(τₚ_list[i]*100)*"%-cumIC", color=colors[j])
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
#folders=["./contacts/results_session0/","./contacts/results_session1/"];
τₚ_list=[.0,.5]#.25,.5,.75,.9]
S0=[4.138758e6, 867417.0, 2.326182e6, 8.084069e6, 3.229145e6, 459761.0, 999280.0, 1.979082e6, 926952.0, 340661.0, 2.381706e6, 2.126254e6, 2.960717e6, 786461.0, 7.478259e6, 781212.0, 2.114588e6, 4.094022e6, 569586.0, 917976.0]
sessions=[74,75,76,77]
Κ_max_capacity_Kilifi=[1e3,1e3,1e4,1e4]

#make_plots_v7(folder,S0)
make_plots_v7("./contacts/results_session0/",S0)
#make_plots_oneExample(sessions,τₚ_list,S0,Κ_max_capacity_Kilifi)
