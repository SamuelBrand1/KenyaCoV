using Plots,MAT, StatsPlots, Statistics#, Images, ImageView,ImageDraw

results_folder=".\\contacts\\results_session40s\\results_session43\\"
τₚ_list=[.0,.25,.5,.75,.9]#[.0,.5]
n_traj=500
dt=.5

colors=[:blue,:orange, :purple3, :maroon, :gold]
markershapes=[:circle,:ltriangle, :hexagon, :star, :plus]
data_files=readdir(results_folder)#["sims100_taup0.0_capacity10000.mat","sims100_taup0.25_capacity10000.mat","sims100_taup0.5_capacity10000.mat","sims100_taup0.75_capacity10000.mat","sims100_taup0.9_capacity10000.mat"];
data_files=[s for s in data_files if endswith(s,".mat")]
riskregionnames = ["Machakos/Muranga" "Mandera" "Baringo/Nakuru" "Nairobi" "Uasin Gishu" "Marsabit" "Garissa" "Nakuru/Narok" "Turkana" "Taita Taveta" "Kitui/Meru" "Kilifi/Mombasa" "Kericho/Kisumu" "Kilifi/Lamu" "Kakamega/Kisumu" "Wajir" "Kajiado/Kisumu" "Homa bay/Migori" "Samburu/Laikipia" "Kilifi/Kwale" "Total"]
wa_coords=[[300,450], [515,85], [165,360], [235,465], [115,300], [300,140], [510,380], [180, 430], [120, 130], [355, 615], [340, 375], [465, 630], [100, 380], [495, 530], [40, 355], [490, 255], [155, 495], [30, 440], [250, 290], [400, 670]]
## Plotting functions IN NAIROBI
function ALL_TAUP_cumulative_symptomatics_sumages(τₚ_list,results_folder,data_files,n_traj)
    p=plot(title="Final cumulative ID and IQ "*string(n_traj)*" sims (sum all ages) Nairobi")
    for i=1:size(τₚ_list,1)
        file = matopen(results_folder*data_files[i])
        data_per_taup=read(file, "sims")
        close(file)
        data_per_taup=data_per_taup[:,2:end]

        Id_per_taup_nairobi=[]
        Iq_per_taup_nairobi=[]
        for j=1:size(data_per_taup,2)
            push!(Id_per_taup_nairobi,sum(data_per_taup[end,j][4,:,8]))
            push!(Iq_per_taup_nairobi,sum(data_per_taup[end,j][4,:,9]))
        end
        plot!(p,Id_per_taup_nairobi, label="ID, taup"*string(τₚ_list[i]),markershape=:circle)
        plot!(p,Iq_per_taup_nairobi, label="IQ, taup"*string(τₚ_list[i]),markershape=:ltriangle)
    end
    #display(p)
    savefig(p,results_folder*"jl_cumulative_ID_IQ_allages_nairobi.png")
end
function ALL_TAUP_cumulative_infecteds_sumages(τₚ_list,results_folder,data_files,n_traj)
    p=plot(title="Final cumulative Ia,Id,Iq "*string(n_traj)*" sims (sum all ages) Nairobi")
    for i=1:size(τₚ_list,1)
        file = matopen(results_folder*data_files[i])
        data_per_taup=read(file, "sims")
        close(file)
        data_per_taup=data_per_taup[:,2:end]

        Ia_per_taup_nairobi=[]
        Id_per_taup_nairobi=[]
        Iq_per_taup_nairobi=[]
        #dead=[]
        for j=1:size(data_per_taup,2)
            push!(Ia_per_taup_nairobi,sum(data_per_taup[end,j][4,:,7]))
            push!(Id_per_taup_nairobi,sum(data_per_taup[end,j][4,:,8]))
            push!(Iq_per_taup_nairobi,sum(data_per_taup[end,j][4,:,9]))
        end
            plot!(p,Ia_per_taup_nairobi, label="IA - taup"*string(τₚ_list[i]),markershape = markershapes[i])
            plot!(p,Id_per_taup_nairobi, label="ID - taup"*string(τₚ_list[i]),markershape = markershapes[i])
            plot!(p,Iq_per_taup_nairobi, label="IQ - taup"*string(τₚ_list[i]),markershape = markershapes[i])
            plot!(p,Ia_per_taup_nairobi .+ Id_per_taup_nairobi .+ Iq_per_taup_nairobi, label="ID+IQ+IA - taup"*string(τₚ_list[i]),markershape=markershapes[i])
    end
    savefig(p,results_folder*"jl_cumulative_Ia, Id, Iq_allages_nairobi.png")
end
function ALL_TAUP_cumulative_infecteds_sumages_ONESIM_GROUPEDBAR(τₚ_list,results_folder,data_files,n_traj)

    histo=Int[]
    L=repeat(["IA","ID","IQ"], inner = size(τₚ_list,1))
    #gro=["taup=" .* string(e)   for e in τₚ_list,j=1:3]  #repeat("taup" .* string.(1:5), outer = 2)
    gro=["taup=0.0","taup=0.5","taup=0.0","taup=0.5","taup=0.0","taup=0.5"]
    for i=1:size(τₚ_list,1)
        file = matopen(results_folder*data_files[i])
        data_per_taup=read(file, "sims")
        close(file)
        data_per_taup=data_per_taup[:,2:end]

        push!(histo,sum(data_per_taup[end,1][4,:,7]))
    end
    for i=1:size(τₚ_list,1)
        file = matopen(results_folder*data_files[i])
        data_per_taup=read(file, "sims")
        close(file)
        data_per_taup=data_per_taup[:,2:end]

        push!(histo,sum(data_per_taup[end,1][4,:,8]))
    end
    for i=1:size(τₚ_list,1)
        file = matopen(results_folder*data_files[i])
        data_per_taup=read(file, "sims")
        close(file)
        data_per_taup=data_per_taup[:,2:end]

        push!(histo,sum(data_per_taup[end,1][4,:,9]))
    end

    p=groupedbar(L,histo,group=gro,title="Final cumulative Ia,Id,Iq ONE SIM (sum all ages) Nairobi")
    savefig(p,results_folder*"jl_cumulative_Ia-Id-Iq_ONE_SIM_allages_nairobi.png")
end
function ALL_TAUP_cumulative_IQ_sumages_BOXPLOT(τₚ_list,results_folder,data_files,n_traj)
    p=boxplot(title="Final cumulative Infecteds IQ "*string(n_traj)*" sims (sum all ages) Nairobi")
    data=[]
    for i=1:size(τₚ_list,1)
        file = matopen(results_folder*data_files[i])
        data_per_taup=read(file, "sims")
        close(file)
        data_per_taup=data_per_taup[:,2:end]

        Id_per_taup_nairobi=[];   # Infecteds=I_D+I_Q
        for j=1:size(data_per_taup,2)
            push!(Id_per_taup_nairobi,sum(data_per_taup[end,j][4,:,9]))
        end
        boxplot!(Id_per_taup_nairobi,label="taup"*string(τₚ_list[i]),ylim=(0,1.5e6))
    end
    #display(p)
    savefig(p,results_folder*"jl_cumulative_IQ_allages_BOXPLOT_nairobi.png")
end
#=function ALL_TAUP_cumulative_symptomatics_perage_BOXPLOT(τₚ_list,results_folder,data_files,n_traj) --> to modify to add cum IQ state9
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
end=#
function ALL_TAUP_sumages_CURVES(τₚ_list,results_folder,data_files,n_traj)
    pI=plot(title="E,IA,ID,IQ "*string(n_traj)*" sims (sum all ages) Nairobi")
    pIdIq=plot(title="ID,IQ,ID+IQ "*string(n_traj)*" sims (sum all ages) Nairobi")
    pSR=plot(title="S and R "*string(n_traj)*" sims (sum all ages) Nairobi")
    pQ=plot(title="Q "*string(n_traj)*" sims (sum all ages) Nairobi")
    pIa=plot(title="IA "*string(n_traj)*" sims (sum all ages) Nairobi")
    data=[]
    for i=1:size(τₚ_list,1)
        file = matopen(results_folder*data_files[i])
        data_per_taup=read(file, "sims")
        close(file)
        data_per_taup=data_per_taup[:,2:end]

        Id_per_taup_nairobi=[];   # Infecteds=I_A+I_D+I_Q
        for j=1:size(data_per_taup,2) # for each sim
            #push!(Id_per_taup_nairobi,sum(data_per_taup[end,j][4,:,8]))
            Sj=[sum([data_per_taup[t,j][4,a,1]  for a=1:16])        for t=1:size(data_per_taup,1)]
            Rj=[sum([data_per_taup[t,j][4,a,6]  for a=1:16])        for t=1:size(data_per_taup,1)]
            plot!(pSR,Sj,color=colors[i],legend=false)
            plot!(pSR,Rj,color=colors[i],legend=false)
            IAj=[sum([data_per_taup[t,j][4,a,3]  for a=1:16])        for t=1:size(data_per_taup,1)]
            IDj=[sum([data_per_taup[t,j][4,a,4]  for a=1:16])        for t=1:size(data_per_taup,1)]
            IQj=[sum([data_per_taup[t,j][4,a,10]  for a=1:16])        for t=1:size(data_per_taup,1)]
            Ej=[sum([data_per_taup[t,j][4,a,2]  for a=1:16])        for t=1:size(data_per_taup,1)]
            plot!(pI,IAj .+ IDj .+ IQj .+ Ej,color=colors[i],legend=false)
            plot!(pIa,IAj,color=colors[i],legend=false)
            plot!(pIdIq,IDj .+ IQj, color=colors[i],legend=false)
            plot!(pIdIq,IDj, color=colors[i],legend=false,markershape=markershapes[i])
            plot!(pIdIq,IQj, color=colors[i],legend=false)#,markershape=markershapes[i])
            Qj=[sum([data_per_taup[t,j][4,a,5]  for a=1:16])        for t=1:size(data_per_taup,1)]
            plot!(pQ,Qj)#, color=colors[i])#,legend=false)
        end
    end
    #display(p)
    savefig(pSR,results_folder*"jl_SR_nairobi.png")
    savefig(pI,results_folder*"jl_EIaIdIq_nairobi.png")
    savefig(pIa,results_folder*"jl_Ia_nairobi.png")
    savefig(pIdIq,results_folder*"jl_IdIq_nairobi.png")
    #savefig(pIdIq_peak,results_folder*"jl_IdIq_PEAK_nairobi.png")
    savefig(pQ,results_folder*"jl_Q_nairobi.png")
end

function ALL_TAUP_sumages_CONTACT_CURVES(τₚ_list,results_folder,data_files,n_traj)
    pCumC=plot(title="CumC "*string(n_traj)*" sims (sum all ages) Nairobi")
    pC=plot(title="Contacts "*string(n_traj)*" sims (sum all ages) Nairobi")
    data=[]
    for i=1:size(τₚ_list,1)
        file = matopen(results_folder*data_files[i])
        data_per_taup=read(file, "sims")
        close(file)
        data_per_taup=data_per_taup[:,2:end]

        for j=1:size(data_per_taup,2) # for each sim
            CumCj=[sum([data_per_taup[t,j][4,a,12]  for a=1:16])        for t=1:size(data_per_taup,1)]
            plot!(pCumC,CumCj, color=colors[i])#,markershape=markershapes[i])
            Cj=[sum([data_per_taup[t,j][4,a,11]  for a=1:16])        for t=1:size(data_per_taup,1)]
            plot!(pC,Cj, color=colors[i])#,markershape=markershapes[i])
        end
    end
    savefig(pCumC,results_folder*"jl_CumC_nairobi.png")
    #savefig(pC,results_folder*"jl_Contacts_nairobi.png")
end
function ALL_TAUP_ONESIM(τₚ_list,results_folder,data_files,n_traj)
    p1=plot()#title="Q - ONE SIM (sum all ages) Nairobi")
    p2=plot()#title="I_Q - ONE SIM (sum all ages) Nairobi")
    for i=1:size(τₚ_list,1)
        file = matopen(results_folder*data_files[i])
        data_per_taup=read(file, "sims")
        close(file)
        data_per_taup=data_per_taup[:,2:end]
        Q=[sum([data_per_taup[t,1][4,a,5]  for a=1:16])        for t=1:size(data_per_taup,1)]
        plot!(p1,Q,label="taup"*string(τₚ_list[i])*"-Q")# color=colors[i],markershape=markershapes[i])
        CumC=[sum([data_per_taup[t,1][4,a,12]  for a=1:16])        for t=1:size(data_per_taup,1)]
        plot!(p1,CumC,label="taup"*string(τₚ_list[i])*"-CumC")
        Iq=[sum([data_per_taup[t,1][4,a,10]  for a=1:16])        for t=1:size(data_per_taup,1)]
        plot!(p2,Iq,label="taup"*string(τₚ_list[i])*"-I_Q")
        CumIq=[sum([data_per_taup[t,1][4,a,9]  for a=1:16])        for t=1:size(data_per_taup,1)]
        plot!(p2,CumIq,label="taup"*string(τₚ_list[i])*"-CumI_Q")
    end
    p=plot(p1,p2,layout=(1,2),title="ONE SIM (sum all ages) Nairobi")
    savefig(p,results_folder*"jl_ONE_SIM_Q_nairobi.png")
end
function peak_times_Nairobi(τₚ_list,results_folder,data_files,n_traj)
    peaks=[]
    median_time_diff=[]
    #p=plot(title="I_A+I_D+I_Q Naitobi")
    for i=1:size(τₚ_list,1)
        peaks_per_tau=[]
        file = matopen(results_folder*data_files[i])
        data_per_taup=read(file, "sims")
        close(file)
        data_per_taup=data_per_taup[:,2:end]
        for j=1:size(data_per_taup,2) # for each sim
            I_A=[sum([data_per_taup[t,j][4,a,3]  for a=1:16])        for t=1:size(data_per_taup,1)]
            I_D=[sum([data_per_taup[t,j][4,a,4]  for a=1:16])        for t=1:size(data_per_taup,1)]
            I_Q=[sum([data_per_taup[t,j][4,a,10]  for a=1:16])        for t=1:size(data_per_taup,1)]
            I=I_A .+ I_D .+ I_Q
            #plot!(p,I,color=colors[i])
            push!(peaks_per_tau,findall(x->x==maximum(I),I)[1] *.5)  # rescale by dt= .5
        end
        push!(peaks,peaks_per_tau)
        push!(median_time_diff,median(peaks_per_tau)-median(peaks[1]))
    end
    #global peaks=[e .* .5 for e in peaks] # rescale by dt= .5
    #median_time_diff=[e*.5 for e in median_time_diff] # rescale by dt= .5
    savefig(boxplot(peaks,title="peak times per tau",legend=false),results_folder*"jl_peaks_nairobi.png")
    display(boxplot(peaks,title="peak times per tau",legend=false))
    savefig(bar(median_time_diff),results_folder*"jl_time_diff_days_nairobi.png")
end
function peak_times_allWA(τₚ_list,results_folder,data_files,n_traj)
    peaks=[]
    #p=plot(title="I_A+I_D+I_Q Naitobi")
    for i=1:size(τₚ_list,1)
        peaks_per_tau_wa=[]
        file = matopen(results_folder*data_files[i])
        data_per_taup=read(file, "sims")
        close(file)
        data_per_taup=data_per_taup[:,2:end]
        for wa=1:20
            for j=1:size(data_per_taup,2) # for each sim
                I_A=[sum([data_per_taup[t,j][wa,a,3]  for a=1:16])        for t=1:size(data_per_taup,1)]
                I_D=[sum([data_per_taup[t,j][wa,a,4]  for a=1:16])        for t=1:size(data_per_taup,1)]
                I_Q=[sum([data_per_taup[t,j][wa,a,10]  for a=1:16])        for t=1:size(data_per_taup,1)]
                I=I_A .+ I_D .+ I_Q
                #plot!(p,I,color=colors[i])
                push!(peaks_per_tau_wa,findall(x->x==maximum(I),I)[1])
            end
            push!(peaks,peaks_per_tau_wa)
        end
    end
    #display(p)
    b=boxplot(peaks,legend=false)
    savefig(b,results_folder*"jl_peaks_allWA.png")

end

## Calling the functions , IN NAIROBI

#ALL_TAUP_cumulative_symptomatics_sumages_BOXPLOT(τₚ_list,results_folder,data_files,n_traj)
#ALL_TAUP_cumulative_symptomatics_perage_BOXPLOT(τₚ_list,results_folder,data_files,n_traj)
#ALL_TAUP_cumulative_symptomatics_peragegroup_BOXPLOT(τₚ_list,results_folder,data_files,n_traj)

#ALL_TAUP_cumulative_infecteds_sumages_ONESIM_GROUPEDBAR(τₚ_list,results_folder,data_files,n_traj)
#ALL_TAUP_cumulative_IQ_sumages_BOXPLOT(τₚ_list,results_folder,data_files,n_traj)
#ALL_TAUP_cumulative_symptomatics_sumages(τₚ_list,results_folder,data_files,n_traj)
#ALL_TAUP_sumages_CURVES(τₚ_list,results_folder,data_files,n_traj)
#ALL_TAUP_sumages_CONTACT_CURVES(τₚ_list,results_folder,data_files,n_traj)
#ALL_TAUP_ONESIM(τₚ_list,results_folder,data_files,n_traj)
#peak_times_Nairobi(τₚ_list,results_folder,data_files,n_traj)
#peak_times_allWA(τₚ_list,results_folder,data_files,n_traj)

## Plotting functions in all areas
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

## Calling the functions in all areas
#maxl,S=ALL_TAUP_cumulative_symptomatics_maps(τₚ_list,results_folder,data_files,n_traj,wa_coords)



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
one_sim_Kilifi(τₚ_list,results_folder,data_files,n_traj,wa_coords)
time_from_intro_Kilifi(τₚ_list,results_folder,data_files,n_traj,wa_coords)
time_from_intro_to_peak_Kilifi(τₚ_list,results_folder,data_files,n_traj,wa_coords)
