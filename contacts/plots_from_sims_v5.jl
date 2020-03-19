using Plots,MAT, StatsPlots, Statistics,JLD2, Images, ImageView,ImageDraw

colors=[:blue, :orange, :purple3, :maroon, :gold]
markershapes=[:cross,:star4, :vcross, :star6, :hline]
riskregionnames = ["Machakos/Muranga" "Mandera" "Baringo/Nakuru" "Nairobi" "Uasin Gishu" "Marsabit" "Garissa" "Nakuru/Narok" "Turkana" "Taita Taveta" "Kitui/Meru" "Kilifi/Mombasa" "Kericho/Kisumu" "Kilifi/Lamu" "Kakamega/Kisumu" "Wajir" "Kajiado/Kisumu" "Homa bay/Migori" "Samburu/Laikipia" "Kilifi/Kwale" "Total"]
wa_coords=[[300,450], [515,85], [165,360], [235,465], [115,300], [300,140], [510,380], [180, 430], [120, 130], [355, 615], [340, 375], [465, 630], [100, 380], [495, 530], [40, 355], [490, 255], [155, 495], [30, 440], [250, 290], [400, 670]]


## KILIFI
function one_sim_Kilifi(τₚ_list,results_folder,data_files,n_traj,wa_coords)
    dt=.5
    p1=plot();    p2=plot();    p3=plot();  p4=plot()
    for i=1:size(τₚ_list,1)
        file = matopen(results_folder*data_files[i])
        data_per_taup=read(file, "sims")
        close(file)
        duration=0:1:(data_per_taup[:,1][end]*dt)
        duration=[string(e) for e in duration]
        data_per_taup=data_per_taup[:,2:end]

        S=[sum([data_per_taup[t,1][12,a,1]  for a=1:16])        for t=1:size(data_per_taup,1)]
        R=[sum([data_per_taup[t,1][12,a,6]  for a=1:16])        for t=1:size(data_per_taup,1)]
        plot!(p1,S,label="taup"*string(τₚ_list[i])*"-S", color=colors[i],legend=false)
        plot!(p1,R,label="taup"*string(τₚ_list[i])*"-R", color=colors[i],legend=false,title="1sim S-R Kilifi")


        I=[sum([data_per_taup[t,1][12,a,10]+data_per_taup[t,1][12,a,3]+data_per_taup[t,1][12,a,4]  for a=1:16])        for t=1:size(data_per_taup,1)]
        plot!(p2,I,label="taup"*string(τₚ_list[i])*"-I", color=colors[i], title="1sim I=A+D+IQ Kilifi")

        Q=[sum([data_per_taup[t,1][12,a,5]  for a=1:16])        for t=1:size(data_per_taup,1)]
        plot!(p3,Q,label="taup"*string(τₚ_list[i])*"-Q", color=colors[i],legend=false,title="1sim Q Kilifi")#,markershape=markershapes[i])

        CumC=[sum([data_per_taup[t,1][12,a,12]  for a=1:16])        for t=1:size(data_per_taup,1)]
        plot!(p4,CumC,label="taup"*string(τₚ_list[i])*"-CumC", color=colors[i],title="1sim CumContacts Kilifi")

        #CumIq=[sum([data_per_taup[t,1][12,a,9]  for a=1:16])        for t=1:size(data_per_taup,1)]
        #plot!(p2,CumIq,label="taup"*string(τₚ_list[i])*"-CumI_Q")
    end
    #p=plot(p1,p2,layout=(1,2))#,title=["ONE SIM S - R Kilifi","ONE SIM Infecteds Kilifi"])#,"ONE SIM Q Kilifi","ONE SIM Cum-contacteds Kilifi"])
    savefig(p1,results_folder*"jl_ONE_SIM_SR_kilifi.png")
    savefig(p2,results_folder*"jl_ONE_SIM_I_kilifi.png")
    #p=plot(p3,p4,layout=(1,2))#,title=["ONE SIM Q Kilifi","ONE SIM Cum-contacteds Kilifi"])
    savefig(p3,results_folder*"jl_ONE_SIM_Q_kilifi.png")
    savefig(p4,results_folder*"jl_ONE_SIM_CumC_kilifi.png")
end
function time_from_intro_Kilifi(τₚ_list,results_folder,data_files,n_traj,wa_coords)
    b=boxplot(title="Time of introduction into Kilifi")
    for i=2:size(τₚ_list,1)
        file = matopen(results_folder*data_files[i])
        data_per_taup=read(file, "sims")
        close(file)
        data_per_taup=data_per_taup[:,2:end]

        ##Time to intro
        t2intro=[]
        for sim=1:size(data_per_taup,2) #for each sim
            Q=[sum([data_per_taup[t,sim][12,a,5]  for a=1:16])        for t=1:size(data_per_taup,1)]
            indices=findall(x -> x>0 , Q)
            if size(indices,1)>0
                push!(t2intro,indices[1]*dt)    #rescaled by dt
            end
        end
        boxplot!(b,t2intro,color=colors[i],label="taup"*string(τₚ_list[i]))
    end
    savefig(b,results_folder*"jl_time2intro_kilifi.png")
end

function time_from_intro_to_peak_Kilifi(τₚ_list,results_folder,data_files,n_traj,wa_coords)

    b=boxplot(title="Time from introduction into Kilifi to peak (in days)")
    for i=2:size(τₚ_list,1)
        file = matopen(results_folder*data_files[i])
        data_per_taup=read(file, "sims")
        close(file)
        data_per_taup=data_per_taup[:,2:end]

        ##Time from intro to peak
        intro2peak=[]
        for sim=1:size(data_per_taup,2) #for each sim
            ##Calculating for each sim, the time of introduction:
            Q=[sum([data_per_taup[t,sim][12,a,5]  for a=1:16])        for t=1:size(data_per_taup,1)]
            indices=findall(x -> x>0 , Q)
            tintro=0
            if size(indices,1)>0
                tintro=indices[1]
            end
            tintro=tintro*dt #rescale by dt
            ##calculating for each sim, the time to peak
            I=[sum([data_per_taup[t,sim][12,a,10]+data_per_taup[t,sim][12,a,3]+data_per_taup[t,sim][12,a,4]  for a=1:16])        for t=1:size(data_per_taup,1)]
            tpeak=findall(x->x==maximum(I),I)[1]*dt     #rescaled by dt
            push!(intro2peak,tpeak-tintro)
        end
        boxplot!(b,intro2peak,color=colors[i],label="taup"*string(τₚ_list[i]),legend=:bottomright)
    end
    savefig(b,results_folder*"jl_intro2peak_kilifi.png")
end

function ALL_TAUP_cumulative_symptomatics_maps(τₚ_list,results_folder,data_files,n_traj,wa_coords)
    imgs_perc=[load("./contacts/kenya.jpg") for τₚ in τₚ_list]
    imgs=[load("./contacts/kenya.jpg") for τₚ in τₚ_list]
    #CList=reshape(range(colorant"yellow",stop=colorant"blue",length=100),1,100)
    CList=reshape(range(colorant"green",stop=colorant"red",length=100),1,100)
    #CList=reshape(range(HSL(colorant"red"), stop=HSL(colorant"green"), length=100),1,100)
    list_Ids_perc=[]
    list_Ids=[]
    S=[]
    for i=1:size(τₚ_list,1)
        file = matopen(results_folder*data_files[i])
        dataf=read(file, "sims")
        close(file)
        dataf=dataf[:,2:end]
        S=[sum(dataf[1,1][wa,:,1])     for wa=1:20]
        Ids_perc=[median([sum(dataf[end,sim][wa,:,8]) for sim=1:size(dataf,2)])*100/S[wa]   for wa=1:20]
        Ids=[median([sum(dataf[end,sim][wa,:,8]) for sim=1:size(dataf,2)])   for wa=1:20]
        max_Ids=maximum(Ids)
        Ids=[median([sum(dataf[end,sim][wa,:,8]) for sim=1:size(dataf,2)])   for wa=1:20]
        for l=1:size(Ids_perc,1)  ##correcting
            if Ids_perc[l]<0        Ids_perc[l]=0
            elseif Ids_perc[l]>100   Ids_perc[l]=100
            end
        end
        push!(list_Ids_perc,Ids_perc)
        push!(list_Ids,Ids)
    end
    #max_perc=[maximum(list_Ids_perc[i]) for i=1:size(τₚ_list,1)]
    maxl=convert(Int64,floor(maximum([maximum(list_Ids[i]) for i=1:size(τₚ_list,1)])))
    CList2=reshape(range(colorant"green",stop=colorant"red",length=maxl),1,maxl)
    for i=1:size(τₚ_list,1)
        for wa=1:20
            draw!(imgs_perc[i], CirclePointRadius(wa_coords[wa][1],wa_coords[wa][2],30),CList[Int64(floor(list_Ids_perc[i][wa]))])
            draw!(imgs[i], CirclePointRadius(wa_coords[wa][1],wa_coords[wa][2],30),CList2[Int64(floor(list_Ids[i][wa]))])
        end
    end
    img0=load("./contacts/kenya.jpg")#[[1 for j=1:628] for i=1:735]#rand(Gray(1.0), 735,628)
    draw!(img0, CirclePointRadius(Point(trunc(Int, size(img0)[1]/2),trunc(Int, size(img0)[2]/2)),500))
    #save(results_folder*"jl_cumulative_symptomatics_map.png",hcat(vcat(imgs[1],img0),vcat(imgs[2],imgs[5]),vcat(imgs[3],imgs[4])))
    #save(results_folder*"jl_cumulative_symptomatics_map_perc.png",hcat(vcat(imgs[1],img0),vcat(imgs_perc[2],imgs_perc[5]),vcat(imgs_perc[3],imgs_perc[4])))
    save(results_folder*"jl_cumulative_ID_map.png",vcat(imgs[1],imgs[2]))
    #save(results_folder*"jl_cumulative_symptomatics_map_perc.png",hcat(vcat(imgs[1],img0),vcat(imgs_perc[2],imgs_perc[5]),vcat(imgs_perc[3],imgs_perc[4])))
    return maxl,S;
end

########
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
    for i=1:size(τₚ_list,1)
        @load results_folder*data_files[i]  sims_vector
        final_cum=[]
        for sim=1:size(sims_vector,1) #for each sim
            push!(final_cum,sum(sims_vector[sim][4][wa,:,1:3])) # final cumulative I=IA+ID+IQ, summed for all ages
        end
        boxplot!(b,final_cum,color=colors[i],label="detection"*string(τₚ_list[i]*100)*"%",legend=:bottomright)
        if i==1     final_cum0detection=median(final_cum)       end
        push!(cumI_diff,final_cum0detection-median(final_cum))
    end
    println(cumI_diff)
    b2=bar(cumI_diff,color=colors,legend=false,title="Gain in final cumulative Infecteds in "*wa_name*" (all ages)",
            xticks = ([1:1:size(τₚ_list,1);], [string(e*100)*"%" for e in τₚ_list]))
    #display(b)
    savefig(b,results_folder*"jl_cumI_"*wa_name*"_box.png")
    savefig(b2,results_folder*"jl_cumI_"*wa_name*"_bar.png")
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
end
########
function make_plots(sessions,τₚ_list,S0)
    for session in sessions
        println("Working plots for session ",session)
        results_folder=".\\contacts\\results_session"*string(Int(floor(session/10)))*"0s\\results_session"*string(session)*"\\"
        data_files=readdir(results_folder)
        data_files=[s for s in data_files if endswith(s,".jld2")]
        @load results_folder*data_files[1]  sims_vector
        n_traj=size(sims_vector,1)

        one_sim_wa(τₚ_list,results_folder,data_files,n_traj,wa_coords,12,"Kilifi")
        one_sim_wa(τₚ_list,results_folder,data_files,n_traj,wa_coords,4,"Nairobi")
        time_intro2peak_wa(τₚ_list,results_folder,data_files,n_traj,wa_coords,12,"Kilifi")
        final_cumI_wa(τₚ_list,results_folder,data_files,n_traj,wa_coords,12,"Kilifi")
        final_cumI_wa(τₚ_list,results_folder,data_files,n_traj,wa_coords,4,"Nairobi")
        final_cumI_map(τₚ_list,results_folder,data_files,n_traj,wa_coords,S0)
        final_cumI_barplot(τₚ_list,results_folder,data_files,n_traj,wa_coords,S0)
    end
end

########
τₚ_list=[.0,.25,.5,.75,.9]
S0=[4.138758e6, 867417.0, 2.326182e6, 8.084069e6, 3.229145e6, 459761.0, 999280.0, 1.979082e6, 926952.0, 340661.0, 2.381706e6, 2.126254e6, 2.960717e6, 786461.0, 7.478259e6, 781212.0, 2.114588e6, 4.094022e6, 569586.0, 917976.0]
#sessions=[51,52,53]
sessions=[75]
make_plots(sessions,τₚ_list,S0)
