using Plots,MAT, StatsPlots, Images, Statistics, ImageView,ImageDraw

results_folder=".\\contacts\\results_session20s\\results_session22\\"
τₚ_list=[0.0,0.25,0.5,0.75,0.9]
data_files=readdir(results_folder)#["sims100_taup0.0_capacity10000.mat","sims100_taup0.25_capacity10000.mat","sims100_taup0.5_capacity10000.mat","sims100_taup0.75_capacity10000.mat","sims100_taup0.9_capacity10000.mat"];
data_files=[s for s in data_files if endswith(s,".mat")]
n_traj=200
riskregionnames = ["Machakos/Muranga" "Mandera" "Baringo/Nakuru" "Nairobi" "Uasin Gishu" "Marsabit" "Garissa" "Nakuru/Narok" "Turkana" "Taita Taveta" "Kitui/Meru" "Kilifi/Mombasa" "Kericho/Kisumu" "Kilifi/Lamu" "Kakamega/Kisumu" "Wajir" "Kajiado/Kisumu" "Homa bay/Migori" "Samburu/Laikipia" "Kilifi/Kwale" "Total"]
wa_coords=[[300,450], [515,85], [165,360], [235,465], [115,300], [300,140], [510,380], [180, 430], [120, 130], [355, 615], [340, 375], [465, 630], [100, 380], [495, 530], [40, 355], [490, 255], [155, 495], [30, 440], [250, 290], [400, 670]]
## Plotting functions IN NAIROBI
function ALL_TAUP_cumulative_symptomatics_sumages(τₚ_list,results_folder,data_files,n_traj)
    p=plot(title="Final cumulative Infecteds(ID+IQ) "*string(n_traj)*" sims (sum all ages) Nairobi")
    for i=1:size(τₚ_list,1)
        file = matopen(results_folder*data_files[i])
        data_per_taup=read(file, "sims")
        close(file)
        data_per_taup=data_per_taup[:,2:end]

        Id_per_taup_nairobi=[]   # Infecteds=I_D+I_Q
        for j=1:size(data_per_taup,2)
            push!(Id_per_taup_nairobi,sum(data_per_taup[end,j][4,:,8]))
        end
        plot!(p,Id_per_taup_nairobi, label="taup"*string(τₚ_list[i]))
    end
    #display(p)
    savefig(p,results_folder*"jl_cumulative_symptomatics_allages_nairobi.png")
end
function ALL_TAUP_cumulative_symptomatics_sumages_BOXPLOT(τₚ_list,results_folder,data_files,n_traj)
    p=boxplot(title="Final cumulative Infecteds(ID+IQ) "*string(n_traj)*" sims (sum all ages) Nairobi")
    data=[]
    for i=1:size(τₚ_list,1)
        file = matopen(results_folder*data_files[i])
        data_per_taup=read(file, "sims")
        close(file)
        data_per_taup=data_per_taup[:,2:end]

        Id_per_taup_nairobi=[];   # Infecteds=I_D+I_Q
        for j=1:size(data_per_taup,2)
            push!(Id_per_taup_nairobi,sum(data_per_taup[end,j][4,:,8]))
        end
        boxplot!(Id_per_taup_nairobi,label="taup"*string(τₚ_list[i]))
    end
    #display(p)
    savefig(p,results_folder*"jl_cumulative_symptomatics_allages_BOXPLOT_nairobi.png")
end
function ALL_TAUP_cumulative_symptomatics_perage_BOXPLOT(τₚ_list,results_folder,data_files,n_traj)
    p=boxplot(layout=16)
    for a=1:16
        data=[]
        for i=1:size(τₚ_list,1)
            file = matopen(results_folder*data_files[i])
            data_per_taup=read(file, "sims")
            close(file)
            data_per_taup=data_per_taup[:,2:end]

            Id_per_taup_nairobi=[];   # Infecteds=I_D+I_Q
            for j=1:size(data_per_taup,2)
                push!(Id_per_taup_nairobi,data_per_taup[end,j][4,a,8])
            end
            boxplot!(Id_per_taup_nairobi,subplot=a,legend=false)#,label="taup"*string(τₚ_list[i])*" a"*string(a),subplot=a)
        end
    end
    #display(p)
    savefig(p,results_folder*"jl_cumulative_symptomatics_perage_BOXPLOT_nairobi.png")
end
function ALL_TAUP_cumulative_symptomatics_peragegroup_BOXPLOT(τₚ_list,results_folder,data_files,n_traj)
    p=boxplot(layout=4)
    for g=1:4#a=1:16
        data=[]
        for i=1:size(τₚ_list,1)
            file = matopen(results_folder*data_files[i])
            data_per_taup=read(file, "sims")
            close(file)
            data_per_taup=data_per_taup[:,2:end]

            Id_per_taup_nairobi=[];   # Infecteds=I_D+I_Q
            for j=1:size(data_per_taup,2)
                a=
                push!(Id_per_taup_nairobi,sum([data_per_taup[end,j][4,g*4-3,8],data_per_taup[end,j][4,g*4-2,8],data_per_taup[end,j][4,g*4-1,8],data_per_taup[end,j][4,g*4,8]]))
            end
            boxplot!(Id_per_taup_nairobi,subplot=g,legend=false,title="Cumulative symptomatics in Nairobi for ages "*string(g*4-3)*"-"*string(g*4),titlefontsize=8)
        end
    end
    #display(p)
    savefig(p,results_folder*"jl_cumulative_symptomatics_peragegroup_BOXPLOT_nairobi.png")
end

## Calling the functions , IN NAIROBI
#=
ALL_TAUP_cumulative_symptomatics_sumages(τₚ_list,results_folder,data_files,n_traj)
ALL_TAUP_cumulative_symptomatics_sumages_BOXPLOT(τₚ_list,results_folder,data_files,n_traj)
ALL_TAUP_cumulative_symptomatics_perage_BOXPLOT(τₚ_list,results_folder,data_files,n_traj)
ALL_TAUP_cumulative_symptomatics_peragegroup_BOXPLOT(τₚ_list,results_folder,data_files,n_traj)
=#

## Plotting functions in all areas
function ALL_TAUP_cumulative_symptomatics_maps(τₚ_list,results_folder,data_files,n_traj,wa_coords)
    imgs_perc=[load("./contacts/kenya.jpg") for τₚ in τₚ_list]
    imgs=[load("./contacts/kenya.jpg") for τₚ in τₚ_list]
    #CList=reshape(range(colorant"yellow",stop=colorant"blue",length=100),1,100)
    CList=reshape(range(colorant"green",stop=colorant"red",length=100),1,100)
    #CList=reshape(range(HSL(colorant"red"), stop=HSL(colorant"green"), length=100),1,100)
    list_Ids_perc=[]
    list_Ids=[]
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
    save(results_folder*"jl_cumulative_symptomatics_map.png",hcat(vcat(imgs[1],img0),vcat(imgs[2],imgs[5]),vcat(imgs[3],imgs[4])))
    save(results_folder*"jl_cumulative_symptomatics_map_perc.png",hcat(vcat(imgs[1],img0),vcat(imgs_perc[2],imgs_perc[5]),vcat(imgs_perc[3],imgs_perc[4])))
end

## Calling the functions in all areas
ALL_TAUP_cumulative_symptomatics_maps(τₚ_list,results_folder,data_files,n_traj,wa_coords)
