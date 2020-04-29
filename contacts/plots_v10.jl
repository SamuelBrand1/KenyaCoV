using Plots,MAT, StatsPlots, Statistics,JLD2,CSV,DataFrames,Printf,ColorBrewer,Colors#, Images, ImageView,ImageDraw
    Plots.default(grid=false,legendfontsize=6, legendfontcolor=:gray30,titlefontsize=9, titlefontcolor=:gray30,xtickfontsize=9, ytickfontsize=9,tickfontcolor=:black,titlefontfamily="Cambria")
    StatsPlots.default(grid=false,legendfontsize=6, legendfontcolor=:gray30,titlefontsize=9, titlefontcolor=:gray30,xtickfontsize=9, ytickfontsize=9,tickfontcolor=:black,linewidth=false)
    colors=[:orange,:purple3,:maroon,:gold,:orangered,:grey,:purple,:ivory3,:chocolate1,:tan1,:rosybrown,:rosybrown2,:brown2,:brown3,:brown4,:deeppink3,:deeppink4]
        #colors2=ColorBrewer.palette("Pastel2",8);colors2=repeat(colors2,outer=[10])
        colors2=[ColorBrewer.palette("Set1",7);ColorBrewer.palette("Set2",7);ColorBrewer.palette("Set3",12);ColorBrewer.palette("Pastel2",8)]
        colors2_subset=[colors2[15],colors2[#=16=#26],colors2[9],colors2[#=15/25=#30],colors2[29],colors2[18]]
        colors2_subset2=[colors2[8],colors2[9],colors2[26],colors2[30]]#,colors2[12]]#,colors2[17],colors2[19]]
        c_R01_I=colors2_subset[2];c_R02_I=colors2_subset[3]
    markershapes=[:cross,:star4, :vcross, :star6, :hline]
    riskregionnames = ["Machakos/Muranga" "Mandera" "Baringo/Nakuru" "Nairobi" "Uasin Gishu" "Marsabit" "Garissa" "Nakuru/Narok" "Turkana" "Taita Taveta" "Kitui/Meru" "Kilifi/Mombasa" "Kericho/Kisumu" "Kilifi/Lamu" "Kakamega/Kisumu" "Wajir" "Kajiado/Kisumu" "Homa bay/Migori" "Samburu/Laikipia" "Kilifi/Kwale" "Total"]
    counties_names=["Baringo","Bomet","Bungoma","Busia","Elgeyo-Marakwet","Embu","Garissa","Homa Bay","Isiolo","Kajiado","Kakamega","Kericho","Kiambu","Kilifi","Kirinyaga","Kisii","Kisumu","Kitui","Kwale","Laikipia","Lamu","Machakos","Makueni","Mandera","Marsabit","Meru","Migori","Mombasa","Murang'a","Nairobi","Nakuru","Nandi","Narok","Nyamira","Nyandarua","Nyeri","Samburu","Siaya","Taita Taveta","Tana River","Tharaka-Nithi","Trans Nzoia","Turkana","Uasin Gishu","Vihiga","Wajir","West Pokot"]
        counties_names=[s[1:min(4,length(s))]   for s in counties_names]
        counties_names=["Barin","Bom","Bung","Busi","Elge","Emb","Garis","Hom","Isiol","Kajia","Kaka","Keric","Kiam","Kilifi","Kirin","Kisii","Kisu","Kitui","Kwal","Laiki","Lam","Mach","Maku","Mand","Mars","Meru","Migo","Mom","Mura","Nairo","Naku","Nand","Naro","Nya","Nyan","Nyer","Sam","Siay","Taita","Tana","Thar","Tran","Turk","Uasi","Vihig","Wajir","West"]
    wa_coords=[[300,450], [515,85], [165,360], [235,465], [115,300], [300,140], [510,380], [180, 430], [120, 130], [355, 615], [340, 375], [465, 630], [100, 380], [495, 530], [40, 355], [490, 255], [155, 495], [30, 440], [250, 290], [400, 670]]
    ages=[string(i)*"-"*string(i+4) for i=0:5:70];push!(ages,"75+")#"0-4" "5-9" "10-14" "15-19" "20-24" "25-29" "30-34" "35-39" "40-44" "45-49" ]

    S0=[4.138758e6, 867417.0, 2.326182e6, 8.084069e6, 3.229145e6, 459761.0, 999280.0, 1.979082e6, 926952.0, 340661.0, 2.381706e6, 2.126254e6, 2.960717e6, 786461.0, 7.478259e6, 781212.0, 2.114588e6, 4.094022e6, 569586.0, 917976.0]
    S012A=[279992, 264585, 246992, 206735, 220761, 211260, 181931, 130091, 111761, 82036, 55221, 42469, 34381, 23781, 16214, 18044]
    #@load "data/data_for_age_structuredmodel_16n_a.jld2" N_region_age;    S0perAge=[sum(N_region_age[:,a]) for a=1:16]; print(S0perAge) ;
    S0perAge=[5993115, 6202467, 6345902, 5285706, 4447468, 3854402, 3570565, 2650024, 2259166, 1786205, 1308564, 1118068, 869994, 658162, 514521, 697759]
    width1=400;length1=400;     width2=1000;length2=800;        ylims_boxplot=(6e5,Inf)
    Rotation=40;    xRotationV=90;   limit=1e2;             δs=[.3,.7]

######## Regions to counties
function regions2counties(data_per_region)
    file = matopen("./contacts\\conversion_matrix.mat")
    conversion_matrix=read(file, "conversion_matrix")
    close(file)

    M=[conversion_matrix[:,wa] .* data_per_region[wa]     for wa=1:20]
    A=[sum([M[i][r]    for i=1:20])     for r=1:47]
    return A
end
#[regions2counties(cumI_medians[i])  for i=1:size(cumI_medians,1)]

####### Cleaning functions
function clean(sims_vector,limit)
    i=1
    while(i<=size(sims_vector,1))
        if sum(sims_vector[i][10])<=limit
            deleteat!(sims_vector,i)
            println("deleted sim=",i)
        else
            i+=1
        end
    end
end
#=function clean_JLDs(folders,limit)
    data_files=[]
    for folder in folders
        data_files2=readdir(folder);
        data_files2=[folder*s for s in data_files2 if endswith(s,".jld2")&&s!="stats.jld2"]
        data_files=[data_files ;data_files2]
    end
    for file in data_files
        @load file sims_vector
        clean(sims_vector)
    end
end=#

#######
######## No intervention plots
#=function kenya_barplot(data_files,plot_folder,S0)
    if !isdir(plot_folder)   mkdir(plot_folder)   end
    wa=12;wa_name="Kilifi"

    final_cum0detection=0;    #=cumIs=[];=#   cumI_medians=[];    cumI_medians_kenya=[];
    cumDeaths_medians=[];cumDeaths_medians_kenya=[];cumDeaths_stds_kenya=[]
    for i=1:size(data_files,1)
        @load data_files[i]  sims_vector
        @load data_files[i]  sessionParams

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

    bar!(barplotALL,cumI_medians2,color=[colors2[1],colors2[6]],linecolor=false,label=false,bar_width=1.005)
    bar!(barplotALL,cumDeaths_medians2,color=:red,linecolor=false,label=false,bar_width=1.005)
    savefig(barplotALL,plot_folder*"jl_cumI_barplot.png")

    bdeaths=bar(cumDeaths_medians2,color=:red,yflip=:true,xaxis=:false,size=(width1*.95,length1*.4),label=false);#display(bdeaths)
    savefig(bdeaths,plot_folder*"jl_cumDeaths_Kenya_barplot.png")
    return cumI_medians_kenya,cumDeaths_stds_kenya,cumI_medians,cumDeaths_medians
end=#
function cumIandDeaths_kenya(data_files,plot_folder,S0)
    if !isdir(plot_folder)   mkdir(plot_folder)   end
    #wa=12;wa_name="Kilifi"

    final_cum0detection=0;    #=cumIs=[];=#   cumI_medians=[];    cumI_medians_kenya=[];cumI_stds_kenya=[]
    cumDeaths_medians=[];cumDeaths_medians_kenya=[];cumDeaths_stds_kenya=[];
    for i=1:size(data_files,1)
        @load data_files[i]  sims_vector;clean(sims_vector,limit)
        @load data_files[i]  sessionParams
        #if i==1      CSV.write(plot_folder*"statsCumIAllSims.csv", DataFrame(cumI=[sum(sims_vector[sim][1][end][:,1:3])   for sim=1:size(sims_vector,1)]))        end

        ##final_cumI_barplot:
        cumI=[[sum(sims_vector[sim][1][end][wa,1:3]) for wa=1:20]   for sim=1:size(sims_vector,1)]
        median_cumI=[median([cumI[sim][wa] for sim=1:size(sims_vector,1)    ])    for wa=1:20]
        push!(cumI_medians,median_cumI)

        ##deaths for barplot:
        cumDeaths=[[sum(sims_vector[sim][4][end][wa]) for wa=1:20]   for sim=1:size(sims_vector,1)]
        median_cumDeaths=[median([cumDeaths[sim][wa] for sim=1:size(sims_vector,1)    ])      for wa=1:20]
        push!(cumDeaths_medians,median_cumDeaths)

        ##All kenya cumI
        cumI_median_kenya=median([sum([cumI[sim][wa]  for wa=1:20])  for sim=1:size(sims_vector,1)   ])
        push!(cumI_medians_kenya,cumI_median_kenya)
        cumI_std_kenya=std([sum([cumI[sim][wa]  for wa=1:20])  for sim=1:size(sims_vector,1) ])
        push!(cumI_stds_kenya,cumI_std_kenya)

        ## All Kenya cumDeaths
        cumDeaths_median_kenya=median([sum([cumDeaths[sim][wa]  for wa=1:20])  for sim=1:size(sims_vector,1)  ])
        cumDeaths_std_kenya=std([sum([cumDeaths[sim][wa]  for wa=1:20])  for sim=1:size(sims_vector,1)  ])
        push!(cumDeaths_medians_kenya,cumDeaths_median_kenya);push!(cumDeaths_stds_kenya,cumDeaths_std_kenya)
    end
    S02=[]; for wa=1:20 push!(S02,S0[wa]);   end
    barplotALL=bar([1.5:2:20*2;],S02,color=[:green,:chartreuse4],xticks=([1:2:20*2;],riskregionnames),fillalpha=.5,linewidth=false,bar_width=2,linecolor=false,
                    size=(width1*1.1,length1),xrotation=xRotationV,label=false,title="Estimated cases and deaths in Kenya without intervention")

    cumI_medians2=[];cumDeaths_medians2=[]# No intervention R01 and δ2 + R02 and δ2
    for wa=1:20,i in [2,4]                push!(cumI_medians2,cumI_medians[i][wa]);  push!(cumDeaths_medians2,cumDeaths_medians[i][wa]);    end
    #for wa=1:20,i=1:size(data_files,1)    push!(cumDeaths_medians2,cumDeaths_medians[i][wa]);end

    bar!(barplotALL,cumI_medians2,color=[c_R01_I,c_R02_I],linecolor=false,label=false,bar_width=1.005)#,fillalpha=.8);
    bar!(barplotALL,cumDeaths_medians2,color=colors2[1],linecolor=false,label=false,bar_width=1.005);#display(barplotALL)
    savefig(barplotALL,plot_folder*"jl_cumI_Kenya_2R0-d10.png")

    cumDeaths_medians3=[]
    for wa=1:20,i=1:size(data_files,1)        push!(cumDeaths_medians3,cumDeaths_medians[i][wa]);    end
    #bdeaths=bar(cumDeaths_medians3,color=[c_R01_I,c_R02_I],yflip=:true,xaxis=:false,size=(width1*.95,length1*.4),label=false);#display(bdeaths)
    #savefig(bdeaths,plot_folder*"jl_cumDeaths_Kenya_barplot2R02d.png")
        #x=[];for i=1:40 push!(x,i);push!(x,i);end
        bdeaths_v2=bar([cumDeaths_medians3[i] for i=2:2:80],color=[c_R01_I,c_R02_I],ls=:solid,yflip=:true,xaxis=:false,label=false,fillalpha=.6,bar_width=1.005)
            bar!(bdeaths_v2,[cumDeaths_medians3[i] for i=1:2:80],color=[c_R01_I,c_R02_I],ls=:dash,xaxis=:false,size=(width1,length1*.4),bar_widths=.6,label=false);#display(bdeaths_v2)
        savefig(bdeaths_v2,plot_folder*"jl_cumDeaths_Kenya_2R0-2d_v2.png")

    # REGIONS TO COUNTIES
    S02_counties=regions2counties(S02)
    cumI_medians_counties=[regions2counties(cumI_medians[i])  for i=1:size(cumI_medians,1)]
    cumDeaths_medians_counties=[regions2counties(cumDeaths_medians[i])  for i=1:size(cumDeaths_medians,1)]
    cumI_medians2=[];cumDeaths_medians2=[]# No intervention R01 and δ2 + R02 and δ2
    for wa=1:47,i in [2,4]                push!(cumI_medians2,cumI_medians_counties[i][wa]);  push!(cumDeaths_medians2,cumDeaths_medians_counties[i][wa]);    end
    barplotALL2=bar([1.5:2:47*2;],S02_counties,color=[:green,:chartreuse4],xticks=([1:2:47*2;],counties_names),fillalpha=.5,linewidth=false,bar_width=2,linecolor=false,
                    widen=true,xrotation=xRotationV,label=false,title="Estimated cases and deaths in Kenya without intervention",size=(width1*1.5,length1*.8))
    bar!(barplotALL2,cumI_medians2,color=[c_R01_I,c_R02_I],linecolor=false,label=false,bar_width=1.005)#,fillalpha=.8);
    bar!(barplotALL2,cumDeaths_medians2,color=colors2[1],linecolor=false,label=false,bar_width=1.005);#display(barplotALL)

    cumDeaths_medians_counties=[regions2counties(cumDeaths_medians[i])  for i=1:size(cumI_medians,1)]
    cumDeaths_medians3=[]
    for wa=1:47,i=1:size(data_files,1)        push!(cumDeaths_medians3,cumDeaths_medians_counties[i][wa]);    end
    bdeaths_v3=bar([cumDeaths_medians3[i] for i=2:2:47*4],color=[c_R01_I,c_R02_I],ls=:solid,yflip=:true,xaxis=:false,label=false,fillalpha=.6,bar_width=1.005)
        bar!(bdeaths_v3,[cumDeaths_medians3[i] for i=1:2:47*4],color=[c_R01_I,c_R02_I],ls=:dash,xaxis=:false,size=(width1*1.48,length1*.4),bar_widths=.6,label=false);#display(bdeaths_v2)

    savefig(barplotALL2,plot_folder*"jl_cumI_Kenya_2R0-d10_counties.png")
    savefig(bdeaths_v3,plot_folder*"jl_cumDeaths_Kenya_2R0-2d_v2_counties.png")
    #display(plot(barplotALL2,bdeaths_v3, layout=(2,1)))

    CSV.write(plot_folder*"statsCumIcounties.csv", DataFrame(#=xticks_labels=xticks_labels[2:end],=#cumI_medians_counties=cumI_medians_counties))
    #@save plot_folder*"statsCumIcounties.jld2" cumI_medians_counties;    @save plot_folder*"statsCumDeathsCounties.jld2" cumDeaths_medians_counties
    return cumI_medians_kenya,cumI_stds_kenya,cumI_medians,cumDeaths_medians_kenya,cumDeaths_stds_kenya,cumDeaths_medians
end

function cumIandDeaths_ages_kenya(data_files,plot_folder,S0)
    if !isdir(plot_folder)   mkdir(plot_folder)   end
    #wa=12;wa_name="Kilifi"

    @load data_files[1]  sims_vector;clean(sims_vector,limit) # No detection R01 and δ1 (only calculating deaths)
    cumDeaths_1=[[] for i=1:16];
    for sim=1:size(sims_vector,1),a=1:16    push!(cumDeaths_1[a],sum(sims_vector[sim][11][:,a])); end

    @load data_files[2]  sims_vector;clean(sims_vector,limit) # No detection R01 and δ2 (calculating cumI and deaths)
    cumDeaths_2=[[] for i=1:16];    cumIsA_2=[[] for i=1:16];
    for sim=1:size(sims_vector,1),a=1:16    push!(cumIsA_2[a],sum(sims_vector[sim][10][:,a,1:3]));end
    for sim=1:size(sims_vector,1),a=1:16    push!(cumDeaths_2[a],sum(sims_vector[sim][11][:,a])); end

    @load data_files[3]  sims_vector;clean(sims_vector,limit) # No detection R02 and δ1 (only calculating deaths)
    cumDeaths_3=[[] for i=1:16];
    for sim=1:size(sims_vector,1),a=1:16    push!(cumDeaths_3[a],sum(sims_vector[sim][11][:,a])); end

    @load data_files[4]  sims_vector;clean(sims_vector,limit) ## No detection R01 and δ2 (calculating cumI and deaths)
    cumDeaths_4=[[] for i=1:16];    cumIsA_4=[[] for i=1:16];
    for sim=1:size(sims_vector,1),a=1:16    push!(cumIsA_4[a],sum(sims_vector[sim][10][:,a,1:3])); end
    for sim=1:size(sims_vector,1),a=1:16    push!(cumDeaths_4[a],sum(sims_vector[sim][11][:,a])); end

    b=bar([1.5:2:16*2;],S0perAge,color=[:green,:chartreuse4],xticks=([1.5:2:16*2;],ages),fillalpha=.5,linewidth=false,bar_width=2,linecolor=false,
                    size=(width1*1.1,length1),xrotation=xRotationV,label=false,title="Estimated cases and deaths per age in Kenya without intervention")

    cumIfinal=[];    for a=1:16  push!(cumIfinal,median(cumIsA_2[a]));push!(cumIfinal,median(cumIsA_4[a])); end # collecting only for δ2 (10%)
    bar!(b,cumIfinal,color=[c_R01_I,c_R02_I],xrotation=xRotationV,contour_labels=true,size=(width1*1.2,length1*.8),linecolor=false,label=false);#display(b)

    cumDeathsfinal=[];    for a=1:16  push!(cumDeathsfinal,median(cumDeaths_2[a]));push!(cumDeathsfinal,median(cumDeaths_4[a])); end  # collecting only for δ2 (10%)
    fillalphas=[cumDeathsfinal[a]==0 ? 1 : 0  for a=1:16]
    bar!(b,cumDeathsfinal,color=:red,xrotation=xRotationV,size=(width1*1.2,length1*.8),linecolor=false,label=false);#display(b)
    savefig(b,plot_folder*"jl_cumI_ages_Kenya.png")

    #plotting deaths
    cumDeathsfinal_4files=[];    for a=1:16  push!(cumDeathsfinal_4files,median(cumDeaths_1[a]));push!(cumDeathsfinal_4files,median(cumDeaths_2[a]));push!(cumDeathsfinal_4files,median(cumDeaths_3[a]));push!(cumDeathsfinal_4files,median(cumDeaths_4[a])); end
    bdeaths=bar([cumDeathsfinal_4files[i] for i=2:2:16*4],color=[c_R01_I,c_R02_I],ls=:solid,yflip=:true,xaxis=:false,fillalpha=.6,label=false);
    bar!(bdeaths,[cumDeathsfinal_4files[i] for i=1:2:16*4],color=[c_R01_I,c_R02_I],ls=:dash,xaxis=:false,size=(width1*1.22,length1*.4),bar_widths=.6,label=false);#display(bdeaths)
    savefig(bdeaths,plot_folder*"jl_cumDeaths_ages_Kenya.png")

    return cumIsA_2,cumIsA_4,cumDeaths_1,cumDeaths_2,cumDeaths_3,cumDeaths_4
end

######## Comparing different interventions
function cumIandDeaths_kenya_COMPAREINTERVENTIONS_counties(data_files,plot_folder,titleaddition)
    if !isdir(plot_folder)   mkdir(plot_folder)   end

    final_cum0detection=0;    #=cumIs=[];=#   cumI_medians=[];    cumI_medians_kenya=[];cumI_stds_kenya=[]
    cumDeaths_medians=[];cumDeaths_medians_kenya=[];cumDeaths_stds_kenya=[]
    for i=1:size(data_files,1)
        @load data_files[i]  sims_vector;clean(sims_vector,limit)
        @load data_files[i]  sessionParams

        ##final_cumI_barplot:
        cumI=[[sum(sims_vector[sim][1][end][wa,1:3]) for wa=1:20]   for sim=1:size(sims_vector,1)]
        median_cumI=[median([cumI[sim][wa] for sim=1:size(sims_vector,1)])      for wa=1:20]
        push!(cumI_medians,median_cumI)

        ##deaths for barplot:
        cumDeaths=[[sum(sims_vector[sim][4][end][wa]) for wa=1:20]   for sim=1:size(sims_vector,1)]
        median_cumDeaths=[median([cumDeaths[sim][wa] for sim=1:size(sims_vector,1)])      for wa=1:20]
        push!(cumDeaths_medians,median_cumDeaths)

        ##All kenya cumI
        cumI_median_kenya=median([sum([cumI[sim][wa]  for wa=1:20])  for sim=1:size(sims_vector,1)])
        push!(cumI_medians_kenya,cumI_median_kenya)
        cumI_std_kenya=std([sum([cumI[sim][wa]  for wa=1:20])  for sim=1:size(sims_vector,1)])
        push!(cumI_stds_kenya,cumI_std_kenya)

        ## All Kenya cumDeaths
        cumDeaths_median_kenya=median([sum([cumDeaths[sim][wa]  for wa=1:20])  for sim=1:size(sims_vector,1)])
        cumDeaths_std_kenya=std([sum([cumDeaths[sim][wa]  for wa=1:20])  for sim=1:size(sims_vector,1)])
        push!(cumDeaths_medians_kenya,cumDeaths_median_kenya);push!(cumDeaths_stds_kenya,cumDeaths_std_kenya)
    end

    # REGIONS TO COUNTIES And Plotting
    cumI_medians_counties=[regions2counties(cumI_medians[i])  for i=1:size(cumI_medians,1)]
        cumDeaths_medians_counties=[regions2counties(cumDeaths_medians[i])  for i=1:size(cumDeaths_medians,1)]

        gainI=[];gainDeaths=[]
        for i=1:size(cumI_medians_counties,1)
            push!(gainI,[cumI_medians_counties[1][r]-cumI_medians_counties[i][r]  for r=1:47])
            push!(gainDeaths,[cumDeaths_medians_counties[1][r]-cumDeaths_medians_counties[i][r]  for r=1:47])
        end

        colors3=repeat([colors2_subset[3],colors2[8]],outer=3);colors3=[colors2[1]; colors3]
        markershapes3=repeat([:diamond,:circle,:xcross],inner=2);markershapes3=[:square ;markershapes3]
        legend3=["No intervention", "90%det in Nairobi","90%det+3mCT in Nairobi","90%det in Nairobi+50%det elsewhere","90%det in Nairobi+50%det elsewhere+3mCT","90%det in all Kenya","90%det+3mCT in all Kenya"]
        legend3=["No intervention", "Detection in Nairobi","Detection and CT in Nairobi","Stronger detection in Nairobi","Stronger detection in Nairobi with CT","Detection in all Kenya","Detection and CT in all Kenya"]
        s=scatter(xrotation=70,xticks=([1:1:size(counties_names,1);],counties_names),title="Estimated cases in Kenya for different interventions "*titleaddition)
        sgain=scatter(xrotation=70,xticks=([1:1:size(counties_names,1);],counties_names),title="Estimated gain in cases in Kenya for different interventions "*titleaddition)
        for i=1:size(cumI_medians_counties,1)
            scatter!(s,[1:1:size(counties_names,1);],cumI_medians_counties[i],color=colors3[i],size=(width1*1.5,length1),
                    markersize=4,markershape=markershapes3[i],label=legend3[i],legend=:topleft,markerstrokewidth=0);
            if i!=1      scatter!(sgain,[1:1:size(counties_names,1);],gainI[i],color=colors3[i],size=(width1*1.5,length1),
                            markersize=4,markershape=markershapes3[i],label=legend3[i],legend=:topleft,markerstrokewidth=0);
            end
        end
        savefig(s,plot_folder*"cumI.png");savefig(sgain,plot_folder*"cumIgain.png");savefig(plot!(sgain,ylims=(0,1e4),legend=false),plot_folder*"cumIgain_0focus.png");
        d=scatter(xrotation=70,xticks=([1:1:size(counties_names,1);],counties_names),title="Estimated death counts in Kenya for different interventions "*titleaddition,legend=false)
        dgain=scatter(xrotation=70,xticks=([1:1:size(counties_names,1);],counties_names),title="Estimated gain in death counts in Kenya for different interventions "*titleaddition,legend=false)
        for i=1:size(cumDeaths_medians_counties,1)
            scatter!(d,[1:1:size(counties_names,1);],cumDeaths_medians_counties[i],color=colors3[i],size=(width1*1.5,length1),
                    markersize=4,markershape=markershapes3[i],label=legend3[i],markerstrokewidth=0);
            if i!=1      scatter!(dgain,[1:1:size(counties_names,1);],gainDeaths[i],color=colors3[i],size=(width1*1.5,length1),
                            markersize=4,markershape=markershapes3[i],label=legend3[i],legend=:topleft,markerstrokewidth=0);
            end
        end
        savefig(d,plot_folder*"cumDeaths.png");savefig(dgain,plot_folder*"cumDeathsgain.png")
    return cumI_medians_counties,gainI,cumDeaths_medians_counties,gainDeaths
end
######## Cases per age in Kilifi : Estimated cases per age in Kilifi R0=2.5,a=50%,90days
function kilifi_ages_α5_90days(data_files,plot_folder,S0)
    if !isdir(plot_folder)   mkdir(plot_folder)   end
    wa=12;wa_name="Kilifi"

    @load data_files[1]  sims_vector;clean(sims_vector,limit)
    cumIsA=[[] for i=1:16];    for sim=1:size(sims_vector,1),a=1:16    push!(cumIsA[a],sum(sims_vector[sim][10][wa,a,1:3]));end

    @load data_files[2]  sims_vector;clean(sims_vector,limit)#
    cumIsA_2=[[] for i=1:16];    for sim=1:size(sims_vector,1),a=1:16    push!(cumIsA_2[a],sum(sims_vector[sim][10][wa,a,1:3])); end

    @load data_files[3]  sims_vector;clean(sims_vector,limit)#
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
    bar!(b,[2:2:32;],cumIfinal3,color=color=colors2[6],rotation=30,bar_width=.5,lw=.5,ls=:dash,label=false)
    savefig(b,plot_folder*"jl_cumI_ages_kilifi.png")
    return cumIsA
end

function kilifi_ages_α5_9_90days_over60s(data_files,plot_folder,S0)
    if !isdir(plot_folder)   mkdir(plot_folder)   end
    wa=12;wa_name="Kilifi"

    @load data_files[1]  sims_vector;clean(sims_vector,limit)
    cumIsA=[[] for i=1:16];    for sim=1:size(sims_vector,1),a=1:16    push!(cumIsA[a],sum(sims_vector[sim][10][wa,a,1:3]));end

    @load data_files[2]  sims_vector;clean(sims_vector,limit)#
    cumIsA_2=[[] for i=1:16];    for sim=1:size(sims_vector,1),a=1:16    push!(cumIsA_2[a],sum(sims_vector[sim][10][wa,a,1:3])); end

    @load data_files[3]  sims_vector;clean(sims_vector,limit)#
    cumIsA_3=[[] for i=1:16];    for sim=1:size(sims_vector,1),a=1:16    push!(cumIsA_3[a],sum(sims_vector[sim][10][wa,a,1:3])); end

    @load data_files[4]  sims_vector;clean(sims_vector,limit)#
    cumIsA_4=[[] for i=1:16];    for sim=1:size(sims_vector,1),a=1:16    push!(cumIsA_4[a],sum(sims_vector[sim][10][wa,a,1:3])); end

    @load data_files[5]  sims_vector;clean(sims_vector,limit)#
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
    bar!(b,[2:3:48;],cumIfinal4,color=colors2[6],rotation=30,bar_width=.5,lw=.5,ls=:dash,label=false)

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
    bar!(b,[2:3:12;],cumIfinal4,color=colors2[6],rotation=30,bar_width=.5,lw=.5,ls=:dash,label=false)

            cumIfinal5=[];
            for a=13:16  push!(cumIfinal5,median(cumIsA_5[a])); end
    bar!(b,[3:3:12;],cumIfinal5,color=colors2[4],rotation=30,bar_width=.5,lw=.5,ls=:dash,label=false)
    savefig(b,plot_folder*"jl_cumI_ages_kilifi_over60s.png")

    ##### Over 60s GAIN
    cumIfinalgain=[];
    for a=13:16  push!(cumIfinalgain,median(cumIsA[a])-median(cumIsA_2[a]));push!(cumIfinalgain,median(cumIsA[a])-median(cumIsA_3[a])); end
    b=bar(cumIfinalgain,color=[colors2[6],colors2[4]],rotation=30,contour_labels=true,size=(width1*1.5,length1*.8),
            ls=:solid,bar_width=1.01,linecolor=[false,:black,:black],label=false,title="Gain in death numbers per age in Kilifi R0=2.5,90days")

            cumIfinal4gain=[];
            for a=13:16  push!(cumIfinal4gain,median(cumIsA[a])-median(cumIsA_4[a])); end
    bar!(b,[1:2:8;],cumIfinal4gain,color=colors2[6],rotation=30,bar_width=.5,lw=.5,ls=:dash,label=false)

            cumIfinal5gain=[];
            for a=13:16  push!(cumIfinal5gain,median(cumIsA[a])-median(cumIsA_5[a])); end
    bar!(b,[2:2:8;],cumIfinal5gain,color=colors2[4],rotation=30,bar_width=.5,lw=.5,ls=:dash,label=false)
    savefig(b,plot_folder*"jl_cumI_ages_kilifi_over60sGAIN.png")
    return cumIsA
end

######## INTERVENTION
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
function mysort(folders)
    data_files=[]
    for folder in folders
        data_files2=readdir(folder);
        data_files2=[folder*s for s in data_files2 if endswith(s,".jld2")&&s!="stats.jld2"]
        data_files=[data_files ;data_files2]
    end
    A=[]
    for i=1:size(data_files,1)
        @load data_files[i]  sessionParams
        A=[A; [[data_files[i] sessionParams.sc_nb]]]
    end
    sort!(A, by = x -> x[2]);
    return [A[i][1] for i=1:size(A,1)]
end
function make_plots_1R0_twoδs_Kenya(folders,plot_folder,S0,δs,reference_files)
    if !isdir(plot_folder)   mkdir(plot_folder)   end

    data_files=mysort(folders)
    n_scenarios=Int(size(data_files,1)/size(folders,1))
    @load data_files[1]  sessionParams
    b7=boxplot(title="Peak delays in Kenya")
    cumI_diffs=[];  cumI_diffs=[];  cumIs=[];  cumIs_std=[];  final_cum0detection=0; #medians_cumI=[];
        peakss=[];peak_diffss=[];    peaks=[];    peak_diffs=[];    median_peak0=0;
        Cpeaks=[];      cumCs=[];   cumICs=[];     cumDeaths=[];    cumDeaths_std=[]
    xticks_labels=["det, CT dur"];
    for i=1:size(data_files,1)
        @load data_files[i]  sims_vector; ;clean(sims_vector,limit);       @load data_files[i]  sessionParams
        if sessionParams.CT_dur[12]!=0          push!(xticks_labels,string(Int(sessionParams.α[4]*100))*"% "*string(Int(sessionParams.CT_dur[12]))*"days")
        else                                    push!(xticks_labels,string(Int(sessionParams.α[4]*100))*"% no CT")    end
        ##final_cumI_wa:
        final_cum=[];   peak=[];    cumC=[];    cumIC=[];   cumDeath=[]
        for sim=1:size(sims_vector,1) #for each sim
            push!(final_cum,sum(sims_vector[sim][10][:,:,1:3])) # final cumulative I=IA+ID+IQ, summed for all ages

            I=[]#[sims_vector[sim][5][t][wa][1] for t=1:size(sims_vector[sim][5],1)]
            for t=1:size(sims_vector[sim][5],1)     push!(I,sum([sims_vector[sim][5][t][wa][1]     for wa=1:20]));  end

            tpeak=findall(x->x==maximum(I),I)[1];   push!(peak,tpeak)

            push!(cumC,sum([sims_vector[sim][2][end][wa][1]    for wa=1:20]));    push!(cumIC,sum([sims_vector[sim][3][end][wa][1]    for wa=1:20]))
            push!(cumDeath,sum([sims_vector[sim][4][end][wa][1]    for wa=1:20]))
        end
        push!(cumIs,median(final_cum)); push!(cumIs_std,std(final_cum)); push!(cumCs,median(cumC));  push!(cumICs,median(cumIC));    push!(cumDeaths,median(cumDeath));  push!(cumDeaths_std,std(cumDeath))

        if i∈reference_files     final_cum0detection=median(final_cum)       end
        push!(cumI_diffs,final_cum0detection-median(final_cum))

        push!(peaks,median(peak))
        if i==1     median_peak0=median(peaks)       end
        push!(peak_diffs,median(peaks)-median_peak0)
        boxplot!(b7,peaks .- median_peak0,legend=false,color=repeat([colors2[i]],5))
    end
    ##Plotting
    boxplot!(b7,xticks=([0:1:size(xticks_labels,1);], xticks_labels),size=(width1,length1),xrotation=Rotation)
        b2_v2=bar(cumI_diffs[n_scenarios+2:end],color=colors2_subset2,ls=:solid,title="Gain in cases in Kenya compared to no intervention - R0="*string(sessionParams.R₀),fillalpha=.8,
                legend=false,xticks=([0:1:size(cumI_diffs[n_scenarios+2:end],1);], [xticks_labels[1];xticks_labels[3:end]]),size=(width1*.4,length1*.4),xrotation=Rotation,ylims=(0,4e6))
        bar!(b2_v2,cumI_diffs[2:n_scenarios],color=colors2_subset2,ls=:dash,lw=0,bar_widths=.5,legend=false,size=(width1,length1),xrotation=Rotation)
        b3_v2=bar(cumIs[1:n_scenarios],color=[colors2_subset[end],colors2_subset[1],colors2_subset[2]],ls=:dash,title="Estimated cases in Kenya - R0="*string(sessionParams.R₀),fillalpha=.8,ylims=(1.5e7,2.3e7),
                    legend=false,xticks=([0:1:size(cumI_diffs[n_scenarios+1:end],1);], [xticks_labels[1];xticks_labels[2:end]]),size=(width1,length1),xrotation=Rotation)
        bar!(b3_v2,cumIs[n_scenarios+1:end],color=[colors2_subset[end],colors2_subset[1],colors2_subset[2]],ls=:solid,lw=0,bar_widths=.5,legend=false,xrotation=Rotation);

        b5=bar(peak_diffs,color=colors2_subset,title="Peak time difference in Kenya",legend=false,xticks=([0:1:size(xticks_labels,1);], xticks_labels),size=(width2,length1),xrotation=Rotation)
        b6=bar(peaks,color=colors2_subset,title="Peak times in Kenya",legend=false,xticks=([0:1:size(xticks_labels,1);], xticks_labels),size=(width2,length1),xrotation=Rotation)
        b8=bar(cumCs[2:n_scenarios],color=colors2_subset,ls=:dash,title="Total number of contacted in Kenya - R0="*string(sessionParams.R₀),fillalpha=.8,
                legend=false,xticks=([0:1:size(cumCs[2:n_scenarios],1);], [xticks_labels[1];xticks_labels[3:end]]),size=(width1,length1),xrotation=Rotation)
            bar!(b8,cumCs[n_scenarios+2:end],color=colors2_subset,ls=:solid,lw=0,bar_widths=.5,legend=false,size=(width1,length1),xrotation=Rotation);
        b9=bar(cumICs[2:n_scenarios],color=colors2_subset,ls=:dash,title="Contacted infecteds (E or Ix) in Kenya - R0="*string(sessionParams.R₀),fillalpha=.8,
                legend=false,xticks=([0:1:size(cumICs[2:n_scenarios],1);], [xticks_labels[1];xticks_labels[3:end]]),size=(width1,length1),xrotation=Rotation)
            bar!(b9,cumICs[n_scenarios+2:end],color=colors2_subset,ls=:solid,lw=0,bar_widths=.5,legend=false,size=(width1,length1),xrotation=Rotation);
        b10=bar(cumDeaths[n_scenarios+1:end],color=[colors2_subset[end],colors2_subset[1],colors2_subset[2]],ls=:solid,title="Estimated death counts in Kenya - R0="*string(sessionParams.R₀),fillalpha=.8,ylims=(1e4,5.1e4),
                legend=false,xticks=([0:1:size(cumDeaths[n_scenarios+1:end],1);], [xticks_labels[1];xticks_labels[2:end]]),size=(width1,length1),xrotation=Rotation)
            bar!(b10,cumDeaths[1:n_scenarios],color=[colors2_subset[end],colors2_subset[1],colors2_subset[2]],ls=:dash,lw=0,bar_widths=.5,legend=false);

        savefig(b2_v2,plot_folder*"jl_cumI_gain_Kenya.png");        savefig(b3_v2,plot_folder*"jl_cumI_Kenya.png");
        savefig(b5,plot_folder*"jl_peak_diffs_Kenya.png");            savefig(b6,plot_folder*"jl_peaks_KenyaB.png")
        savefig(b7,plot_folder*"jl_peaks_Kenya.png");
        savefig(b8,plot_folder*"jl_cumC_Kenya.png");                  savefig(b9,plot_folder*"jl_cumIC_Kenya_bar.png");
        savefig(b10,plot_folder*"jl_cumDeaths_Kenya.png");

        CSV.write(plot_folder*"stats.csv", DataFrame(xticks_labels=xticks_labels[2:end],cumI_diffs=cumI_diffs, cumIs=cumIs,cumIs_std=cumIs_std,
                    peak_diffs=peak_diffs,peaks=peaks#=,Cpeaks=Cpeaks=#,cumDeaths=cumDeaths))
    return cumI_diffs,cumDeaths,cumDeaths_std
end

function make_plots_1R0_twoδs_regionKilifi(folders,plot_folder,S0,δs)
    if !isdir(plot_folder)   mkdir(plot_folder)   end
    wa=12;wa_name="Kilifi"
    data_files=[]
    for folder in folders
        data_files2=readdir(folder);
        data_files2=[folder*s for s in data_files2 if endswith(s,".jld2")&&s!="stats.jld2"]
        data_files=[data_files ;data_files2]
    end
    data_files=mysort(data_files)
    n_scenarios=Int(size(data_files,1)/size(folders,1))
    @load data_files[1]  sessionParams
    #b=boxplot(title="Estimated cases in "*wa_name*" - R0="*string(sessionParams.R₀))
    #b_v2=boxplot(title="Estimated cases in "*wa_name*" - R0="*string(sessionParams.R₀))
    #b4=boxplot(title="Gain in number of cases in "*wa_name)
    b7=boxplot(title="Peak delays in "*wa_name)

    cumI_diffs=[];peakss=[];peak_diffss=[];    peaks=[];    peak_diffs=[];    median_peak0=0;   Cpeaks=[]
    final_cum0detection=0;        cumI_diffs=[];        cumIs=[];    medians_cumI=[]
    cumCs=[];cumICs=[];           cumDeaths=[]

    xticks_labels=["det, CT dur"];
    for i=1:size(data_files,1)
        @load data_files[i]  sims_vector;   ;clean(sims_vector,limit);     @load data_files[i]  sessionParams
        if sessionParams.CT_dur[12]!=0          push!(xticks_labels,string(Int(sessionParams.α[wa]*100))*"% "*string(Int(sessionParams.CT_dur[12]))*"days")
        else                                    push!(xticks_labels,string(Int(sessionParams.α[wa]*100))*"% no CT")    end
        ##final_cumI_wa:#peaks code
        final_cum=[];           peak=[];    #Cpeak=[]
        cumC=[];cumIC=[];cumDeath=[]
        for sim=1:size(sims_vector,1) #for each sim
            push!(final_cum,sum(sims_vector[sim][10][wa,:,1:3])) # final cumulative I=IA+ID+IQ, summed for all ages

            I=[sims_vector[sim][5][t][wa][1] for t=1:size(sims_vector[sim][5],1)]
            tpeak=findall(x->x==maximum(I),I)[1];   push!(peak,tpeak)

            push!(cumC,sims_vector[sim][2][end][wa][1]);    push!(cumIC,sims_vector[sim][3][end][wa][1])
            push!(cumDeath,sims_vector[sim][4][end][wa][1])
        end
        push!(cumIs,median(final_cum)); push!(cumCs,median(cumC));  push!(cumICs,median(cumIC));    push!(cumDeaths,median(cumDeath))

        #boxplot!(b,final_cum,color=repeat([colors2[i]],5),legend=false)
        #if i<=9 boxplot!(b_v2,final_cum,color=repeat([colors2[i]],5),legend=false,ls=:dash)
        #else    boxplot!(b_v2,final_cum,color=repeat([colors2[i]],5),legend=false,ls=:solid);   end

        if i==1     final_cum0detection=median(final_cum)       end
        push!(cumI_diffs,final_cum0detection-median(final_cum))
        #boxplot!(b4,final_cum0detection .- final_cum,color=repeat([colors2[i]],5),legend=false)

        ##final_cumI_barplot:
        cumI=[[sum(sims_vector[sim][1][end][wa,1:3]) for wa=1:20]   for sim=1:size(sims_vector,1)]
        median_cumI=[median([cumI[sim][wa] for sim=1:size(sims_vector,1)])      for wa=1:20]
        push!(medians_cumI,median_cumI)

        push!(peaks,median(peak))
        if i==1     median_peak0=median(peaks)       end
        push!(peak_diffs,median(peaks)-median_peak0)
        boxplot!(b7,peaks .- median_peak0,legend=false,color=repeat([colors2[i]],5))
    end
    ##final_cumI_wa:
    CSV.write(plot_folder*"stats.csv", DataFrame(xticks_labels=xticks_labels[2:end],cumI_diffs=cumI_diffs, cumIs=cumIs,
                peak_diffs=peak_diffs,peaks=peaks#=,Cpeaks=Cpeaks=#,cumDeaths=cumDeaths))
    ##Plotting
    #boxplot!(b, xticks=([0:1:size(xticks_labels,1);], xticks_labels),size=(width1,length1),xrotation=Rotation)#,ylims=(1.04e6,Inf))
            #boxplot!(b_v2, xticks=([0:1:size(xticks_labels,1);], xticks_labels),size=(width1,length1),xrotation=Rotation,ylims=(6.5e5,Inf))

        #boxplot!(b4,xticks=([0:1:size(xticks_labels,1);], xticks_labels),size=(width1,length1),xrotation=Rotation)
        boxplot!(b7,xticks=([0:1:size(xticks_labels,1);], xticks_labels),size=(width1,length1),xrotation=Rotation)

        #b2=bar(cumI_diffs[2:end],color=colors2[1:4],title="Gain in number of cases in "*wa_name*" - R0="*string(sessionParams.R₀),linestyles=[[:dot for i=1:4] ;[:solid for i=1:5]],lw=0,
        #        legend=false,xticks=([0:1:size(xticks_labels,1);], xticks_labels),size=(width1,length1),xrotation=Rotation)
            b2_v2=bar(cumI_diffs[n_scenarios+2:end],color=colors2[1:4],ls=:solid,title="Gain in number of cases in "*wa_name*" - R0="*string(sessionParams.R₀),fillalpha=.8,
                    legend=false,xticks=([0:1:size(cumI_diffs[n_scenarios+2:end],1);], [xticks_labels[1];xticks_labels[3:end]]),size=(width1*.4,length1*.4),xrotation=Rotation)
            bar!(b2_v2,cumI_diffs[2:n_scenarios],color=colors2[1:4],ls=:dash,lw=0,bar_widths=.5,legend=false,size=(width1,length1),xrotation=Rotation)
        #b3=bar(cumIs,color=colors2,title="Total number of cases in "*wa_name,legend=false,xticks=([0:1:size(xticks_labels,1);], xticks_labels),size=(width2,length1),xrotation=Rotation)

            b3_v2=bar(cumIs[2:n_scenarios],color=colors2,ls=:dash,title="Estimated cases in "*wa_name*" - R0="*string(sessionParams.R₀),fillalpha=.8,
                    legend=false,xticks=([0:1:size(cumI_diffs[n_scenarios+2:end],1);], [xticks_labels[1];xticks_labels[3:end]]),size=(width1*.4,length1*.4),xrotation=Rotation)
            bar!(b3_v2,cumIs[n_scenarios+2:end],color=colors2,ls=:solid,lw=0,bar_widths=.5,legend=false,size=(width1,length1),xrotation=Rotation);

        b5=bar(peak_diffs,color=colors2,title="Peak time difference in "*wa_name,legend=false,xticks=([0:1:size(xticks_labels,1);], xticks_labels),size=(width2,length1),xrotation=Rotation)
        b6=bar(peaks,color=colors2,title="Peak times in "*wa_name,legend=false,xticks=([0:1:size(xticks_labels,1);], xticks_labels),size=(width2,length1),xrotation=Rotation)
        b8=bar(cumCs[2:n_scenarios],color=colors2,ls=:dash,title="Total number of contacted in "*wa_name*" - R0="*string(sessionParams.R₀),fillalpha=.8,
                legend=false,xticks=([0:1:size(cumCs[2:n_scenarios],1);], [xticks_labels[1];xticks_labels[3:end]]),size=(width1,length1),xrotation=Rotation)
            bar!(b8,cumCs[n_scenarios+2:end],color=colors2,ls=:solid,lw=0,bar_widths=.5,legend=false,size=(width1,length1),xrotation=Rotation);
        b9=bar(cumICs[2:n_scenarios],color=colors2,ls=:dash,title="Contacted infecteds (E or Ix) in "*wa_name*" - R0="*string(sessionParams.R₀),fillalpha=.8,
                legend=false,xticks=([0:1:size(cumICs[2:n_scenarios],1);], [xticks_labels[1];xticks_labels[3:end]]),size=(width1,length1),xrotation=Rotation)
            bar!(b9,cumICs[n_scenarios+2:end],color=colors2,ls=:solid,lw=0,bar_widths=.5,legend=false,size=(width1,length1),xrotation=Rotation);
        b10=bar(cumDeaths[n_scenarios+1:end],color=colors2,ls=:solid,title="Total deaths in "*wa_name*" - R0="*string(sessionParams.R₀),fillalpha=.8,
                legend=false,xticks=([0:1:size(cumDeaths[n_scenarios+1:end],1);], [xticks_labels[1];xticks_labels[2:end]]),size=(width1,length1),xrotation=Rotation)
            bar!(b10,cumDeaths[1:n_scenarios],color=colors2,ls=:dash,lw=0,bar_widths=.5,legend=false,size=(width1,length1),xrotation=Rotation);

        #=savefig(b,plot_folder*"jl_cumI_"*wa_name*"_box.png");        savefig(b_v2,plot_folder*"jl_cumI_"*wa_name*"_box_v2.png");=#
        #=savefig(b2,plot_folder*"jl_cumI_gain_"*wa_name*"_bar.png");=#  savefig(b2_v2,plot_folder*"jl_cumI_gain_"*wa_name*"_bar_v2.png")
        #=savefig(b3,plot_folder*"jl_cumI_"*wa_name*"_bar.png");=#       savefig(b3_v2,plot_folder*"jl_cumI_"*wa_name*"_bar_v2.png");
        #savefig(b4,plot_folder*"jl_cumI_gain_"*wa_name*"_box.png")
        savefig(b5,plot_folder*"jl_peak_diffs_"*wa_name*"_bar.png"); savefig(b6,plot_folder*"jl_peaks_"*wa_name*"_bar.png")
        savefig(b7,plot_folder*"jl_peaks_"*wa_name*"_box.png");

        savefig(b8,plot_folder*"jl_cumC_"*wa_name*"_bar.png");       savefig(b9,plot_folder*"jl_cumIC_"*wa_name*"_bar.png");
        savefig(b10,plot_folder*"jl_cumDeaths_"*wa_name*"_bar.png");
end

######### INFECTEDS AND CT AREA
function oneCurve_CT_comparedtoNoIntervention(no_int_data_file,data_file,plot_folder)
    if !isdir(plot_folder)   mkdir(plot_folder)   end
    wa=12;wa_name="Kilifi"
    @load data_file  sessionParams
    @load data_file  sims_vector;clean(sims_vector,limit)
    sim=10
    if no_int_data_file==""    ##do not add no intervention red curve to the plots
        IKilifi=[sims_vector[sim][5][t][12][1] for t=1:366]
        p=plot(IKilifi,fillalpha=.7,label="Infected",size=(width1*.8,length1*.7),xlims=(0,300),
                title="Example of an epidemic curve in Kilifi with 3 months \ncontact tracing, R0="*string(sessionParams.R₀)*",det="*string(Int(sessionParams.δ*100))*"%,a="*string(Int(sum(sessionParams.α)*100))*"%")
        CT_start=findfirst(x->x!=0,[sims_vector[sim][2][t][12][1] for t=1:366])
        vspan!(p,[CT_start, CT_start+30], linecolor = false, fillcolor = colors2[5], fillalpha=.4,label="1st month CT")
        vspan!(p,[CT_start+31, CT_start+60], linecolor = false, fillcolor = colors2[4], fillalpha=.4,label="2nd month CT")
        vspan!(p,[CT_start+61, CT_start+90], linecolor = false, fillcolor = colors2[1], fillalpha=.4,label="3rd month CT")
        savefig(p,plot_folder*"jl_exampleI_"*wa_name*"_R0="*string(sessionParams.R₀)*",det="*string(Int(sessionParams.δ*100))*",a="*string(Int(sum(sessionParams.α)*100))*".png");
    else                         ##Add no intervention red curve to the plots
        @load no_int_data_file  sessionParams
        @load no_int_data_file  sims_vector;clean(sims_vector,limit)
        IKilifi_nointervention=[sims_vector[sim][5][t][12][1] for t=1:366]
        p=plot(IKilifi_nointervention,fillalpha=.7,label="No interv",size=(width1*.8,length1*.7),xlims=(0,300),color=:red,
                title="Example of an epidemic curve in Kilifi with 3 months \ncontact tracing, R0="*string(sessionParams.R₀)*",det="*string(Int(sessionParams.δ*100))*"%,a="*string(Int(sum(sessionParams.α)*100))*"%")
        @load data_file  sessionParams
        @load data_file  sims_vector;clean(sims_vector,limit)
        IKilifi=[sims_vector[sim][5][t][12][1] for t=1:366]
        plot!(p,IKilifi,fillalpha=.7,label="Infected",size=(width1*.8,length1*.7),xlims=(0,300),color=:blue)
        CT_start=findfirst(x->x!=0,[sims_vector[sim][2][t][12][1] for t=1:366])
        vspan!(p,[CT_start, CT_start+30], linecolor = false, fillcolor = colors2[5], fillalpha=.4,label="1st month CT")
        vspan!(p,[CT_start+31, CT_start+60], linecolor = false, fillcolor = colors2[4], fillalpha=.4,label="2nd month CT")
        vspan!(p,[CT_start+61, CT_start+90], linecolor = false, fillcolor = colors2[1], fillalpha=.4,label="3rd month CT")
        savefig(p,plot_folder*"jl_exampleI_"*wa_name*"_R0="*string(sessionParams.R₀)*",det="*string(Int(sessionParams.δ*100))*",a="*string(Int(sum(sessionParams.α)*100))*"_v2.png");
    end
end
function CT_comparedtoNoIntervention(no_int_data_file,data_file,plot_folder)
    if !isdir(plot_folder)   mkdir(plot_folder)   end
    #wa=12;wa_name="Kilifi"

    @load no_int_data_file  sessionParams
    @load no_int_data_file  sims_vector;clean(sims_vector,limit)
    p=plot(size=(width1*.8,length1*.75),xlims=(0,300))

    for sim=1:size(sims_vector,1)
        IKenya_nointervention=[sum([sims_vector[sim][5][t][wa][1]   for wa=1:20]) for t=1:366]
        if sim==1       plot!(p,IKenya_nointervention,linealpha=.3,color=:red,label="No interv")
        else            plot!(p,IKenya_nointervention,linealpha=.3,xlims=(0,300),color=:red,label=false)   end
    end
    @load data_file  sessionParams
    @load data_file  sims_vector;clean(sims_vector,limit)
    plot!(p,title="Stochastic epidemic curves in Kenya with 3 months \ncontact tracing, R0="*string(sessionParams.R₀)*",sympt="*string(Int(sessionParams.δ*100))*"%,det="*string(Int(ceil(sum(sessionParams.α)*100)))*"%")
    for sim=1:size(sims_vector,1)
        IKenya=[sum([sims_vector[sim][5][t][wa][1]   for wa=1:20]) for t=1:366]
        if sim==1       plot!(p,IKenya,linealpha=.5,size=(width1*.8,length1*.7),xlims=(0,300),color=:blue,label="Intervention")
        else            plot!(p,IKenya,linealpha=.5,size=(width1*.8,length1*.7),xlims=(0,300),color=:blue,label=false); end
        if sim==1
            CT_start=findfirst(x->x!=0,[sum([sims_vector[sim][2][t][wa][1]   for wa=1:20]) for t=1:366])
            #vspan!(p,[CT_start,     CT_start+30], linecolor = false, fillcolor = colors2[5], fillalpha=.4,label="1st month CT")
            #vspan!(p,[CT_start+31,  CT_start+60], linecolor = false, fillcolor = colors2[4], fillalpha=.4,label="2nd month CT")
            #vspan!(p,[CT_start+61,  CT_start+90], linecolor = false, fillcolor = colors2[1], fillalpha=.4,label="3rd month CT")
        end
    end
    savefig(p,plot_folder*"jl_I_Kenya_R0="*string(sessionParams.R₀)*",det="*string(Int(sessionParams.δ*100))*",a="*string(Int(ceil(sum(sessionParams.α)*100)))*"_v2.png");
end

function getRibbonFromCurves(curves)
    curves_min=[];curves_max=[];curves_median=[];
    for t=1:size(curves[1],1)
        push!(curves_min,   minimum([curves[sim][t]   for sim=1:size(curves,1)]))
        push!(curves_max,   maximum([curves[sim][t]   for sim=1:size(curves,1)]))
        push!(curves_median,median( [curves[sim][t]   for sim=1:size(curves,1)]))
    end
    return curves_min,curves_max,curves_median
end
function CT_comparedtoNoIntervention_area(no_int_data_file,data_file,plot_folder,withLabel)
    if !isdir(plot_folder)   mkdir(plot_folder)   end
    #wa=12;wa_name="Kilifi"

    @load no_int_data_file  sessionParams;    @load no_int_data_file  sims_vector;clean(sims_vector,limit)

    IKenya_nointervention_list=[]
    for sim=1:size(sims_vector,1)        push!(IKenya_nointervention_list,[sum([sims_vector[sim][5][t][wa][1]   for wa=1:20]) for t=1:366]);    end
    IKenya_nointervention_min,IKenya_nointervention_max,IKenya_nointervention_median=getRibbonFromCurves(IKenya_nointervention_list)

    @load data_file  sessionParams;    @load data_file  sims_vector;clean(sims_vector,limit)
    CT_start_list=[];IKenya_list=[]
    for sim=1:size(sims_vector,1)
        push!(IKenya_list,[sum([sims_vector[sim][5][t][wa][1]   for wa=1:20]) for t=1:366])
        push!(CT_start_list,findfirst(x->x!=0,[sum([sims_vector[sim][2][t][wa][1]   for wa=1:20]) for t=1:366]))
    end
    IKenya_min,IKenya_max,IKenya_median=getRibbonFromCurves(IKenya_list)

    #Plotting:
    p=plot(size=(width1*.8,length1*.75),xlims=(0,300),title="Stochastic epidemic curves in Kenya with 3 months \ncontact tracing, R0="*string(sessionParams.R₀)*",sympt="*string(Int(sessionParams.δ*100))*"%,det="*string(Int(ceil(sum(sessionParams.α[4])*100)))*"%")
    vspan!(p,[minimum(CT_start_list),   maximum(CT_start_list)+30],ls=:dash,   lc=colors2_subset[2],linealpha=1,fillcolor=colors2_subset[2],fillalpha=.35,label="1st month CT")
    vspan!(p,[minimum(CT_start_list)+31,maximum(CT_start_list)+60],ls=:dashdot,lc=colors2_subset[4],linealpha=1,fillcolor=colors2_subset[4],fillalpha=.35,label="2nd month CT")
    vspan!(p,[minimum(CT_start_list)+61,maximum(CT_start_list)+90],ls=:dot,    lc=colors2_subset[1],linealpha=1,fillcolor=colors2_subset[1],fillalpha=.35,label="3rd month CT")
    plot!(p,IKenya_nointervention_median,ribbon=(IKenya_nointervention_min,IKenya_nointervention_max),color=:red,label="No intervention")
    plot!(p,IKenya_median,ribbon=(IKenya_min,IKenya_max),color=:blue,label="Intervention",legend=withLabel,ylims=(0,3.1e6))#;display(p)
    savefig(p,plot_folder*"jl_I_Kenya_R0_"*string(sessionParams.R₀)*"_d"*string(Int(sessionParams.δ*100))*"_a"*string(Int(ceil(sum(sessionParams.α[4])*100)))*"_v2.png");
    return IKenya_nointervention_median,IKenya_median
end

########### Verify scenarios
function curves(folders,plot_folder)
    curve_names=["cumIs","cumC","cumIc","cumDeaths","I","Iq","Q","S","R"]
    if !isdir(plot_folder)   mkdir(plot_folder)   end
    for folder in folders
        data_files=readdir(folder)
        data_files=[s for s in data_files if endswith(s,".jld2")&&s!="stats.jld2"]
        for data_file in data_files
            plot_folder2=plot_folder*data_file[1:end-5]*"/"
            if !isdir(plot_folder2)   mkdir(plot_folder2)   end
            @load folder*data_file  sims_vector;clean(sims_vector,limit)
            for wa=1:size(riskregionnames,2)-1
                p_curves=[plot(size=(width1,length1)) for i=1:9]
                for sim=1:size(sims_vector,1),i=1:9
                    d=[sims_vector[sim][i][t][wa][1] for t=1:366]
                    plot!(p_curves[i],d,color=colors2[i])
                end
                for i=1:9   savefig(p_curves[i],plot_folder2*string(i)*"_"*curve_names[i]*"_wa"*string(wa)*"_"*replace(riskregionnames[wa],"/"=>"-")*".png")   end
            end
        end
    end
end
