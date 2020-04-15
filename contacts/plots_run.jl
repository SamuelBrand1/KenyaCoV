include("plots_v10.jl")

####### Kenya per region and with No Intervention
data_files=["./contacts/session30_50sims/sims_sc3000.jld2","./contacts/session31_50sims/sims_sc3100.jld2","./contacts/session32_50sims/sims_sc3200.jld2","./contacts/session33_50sims/sims_sc3300.jld2"]
    plot_folder="./contacts/session30_50sims/plots_sessions30-33_det0/"
    cumI_medians_kenya,cumI_stds_kenya,cumI_medians,cumDeaths_medians_kenya,cumDeaths_stds_kenya,cumDeaths_medians=cumIandDeaths_kenya(data_files,plot_folder,S0)
    println("Cases in all Kenya for R0=2.5 and R0=3, with sympt rate δ=10% is: ",cumI_medians_kenya[2],"⨦",cumI_stds_kenya[2], " (",100*cumI_medians_kenya[2]/sum(S0),"% of the total Kenyan population) and ",cumI_medians_kenya[4],"⨦",cumI_stds_kenya[4], " (",100*cumI_medians_kenya[4]/sum(S0),"%)")

####### Kenya per ages with No Intervention
data_files=["./contacts/session30_50sims/sims_sc3000.jld2","./contacts/session31_50sims/sims_sc3100.jld2","./contacts/session32_50sims/sims_sc3200.jld2","./contacts/session33_50sims/sims_sc3300.jld2"]
    plot_folder="./contacts/session30_50sims/plots_sessions30-33_det0/"
    cumIsA_2,cumIsA_4,cumDeaths_1,cumDeaths_2,cumDeaths_3,cumDeaths_4=cumIandDeaths_ages_kenya(data_files,plot_folder,S0)

    cumIsA_over60s_2=[];cumIsA_over60s_4=[];
    for sim=1:3
        I=0;    for a=15:16     I+=cumIsA_2[a][sim];    end;        push!(cumIsA_over60s_2,I)
        I=0;    for a=15:16     I+=cumIsA_4[a][sim];    end;        push!(cumIsA_over60s_4,I)
    end
    cumDeathsA_over60s_1=[];cumDeathsA_over60s_2=[];cumDeathsA_over60s_3=[];cumDeathsA_over60s_4=[]
    for sim=1:3
        d=0;    for a=15:16     d+=cumDeaths_1[a][sim];    end;        push!(cumDeathsA_over60s_1,d)
        d=0;    for a=15:16     d+=cumDeaths_2[a][sim];    end;        push!(cumDeathsA_over60s_2,d)
        d=0;    for a=15:16     d+=cumDeaths_3[a][sim];    end;        push!(cumDeathsA_over60s_3,d)
        d=0;    for a=15:16     d+=cumDeaths_4[a][sim];    end;        push!(cumDeathsA_over60s_4,d)
    end
    println("For R0=2.5 and 3, and δ=1.%, we estimate ",median(cumIsA_over60s_2),"⨦",std(cumIsA_over60s_2)," (",100*median(cumIsA_over60s_2)/sum(S0perAge[13:end]),"%) and ",median(cumIsA_over60s_4),"⨦",std(cumIsA_over60s_4)," (",100*median(cumIsA_over60s_4)/sum(S0perAge[13:end]),"%)of cases are in the over 60s age groups.")
    println("As for fatalities, deaths in the over 60s is estimated as: \n\tfor R0=2.5 and δ=5% ",median(cumDeathsA_over60s_1),"⨦",std(cumDeathsA_over60s_1)," (",100*median(cumDeathsA_over60s_1)/sum(S0perAge[13:end]),"%)",
                                                                        "\n\tfor R0=2.5 and δ=10% ",median(cumDeathsA_over60s_2),"⨦",std(cumDeathsA_over60s_2)," (",100*median(cumDeathsA_over60s_2)/sum(S0perAge[13:end]),"%)",
                                                                        "\n\tfor R0=3 and δ=5% ",median(cumDeathsA_over60s_3),"⨦",std(cumDeathsA_over60s_3)," (",100*median(cumDeathsA_over60s_3)/sum(S0perAge[13:end]),"%)",
                                                                        "\n\tfor R0=10 and δ=10% ",median(cumDeathsA_over60s_4),"⨦",std(cumDeathsA_over60s_4)," (",100*median(cumDeathsA_over60s_4)/sum(S0perAge[13:end]),"%)")

######## Detection without contact tracing (0CT gain)
#R0=2.5
folders=["./contacts/session30_50sims/CT0/","./contacts/session31_50sims/CT0/"];plot_folder="./contacts/session30_50sims/CT0/";reference_files=[1,4]
    cumI_diffs,cumDeaths,cumDeaths_std=make_plots_1R0_twoδs_Kenya(folders,plot_folder,S0,δs,reference_files)
    println("R0=2.5. Gain in cases with detection and no CT:\n\tδ=5% and α=50% and 90%: ",cumI_diffs[2:3],"\n\tδ=10% and α=50% and 90%: ",cumI_diffs[5:6])
    println("R0=2.5.   Deaths with detection and no CT:\n\tδ=0%, 5% and α=50% and 90%: ",[string(cumDeaths[i])*"⨦"*string(cumDeaths_std[i])  for i=1:3],"\n\tδ=10% and α=0%, 50% and 90%: ",[string(cumDeaths[i])*"⨦"*string(cumDeaths_std[i])  for i=4:6])

#R0=3
folders=["./contacts/session32_50sims/CT0/","./contacts/session33_50sims/CT0/"];plot_folder="./contacts/session32_50sims/CT0/";reference_files=[1,4]
    cumI_diffs,cumDeaths,cumDeaths_std=make_plots_1R0_twoδs_Kenya(folders,plot_folder,S0,δs,reference_files)
    println("R0=3.   Gain in cases with detection and no CT:\n\tδ=5% and α=50% and 90%: ",cumI_diffs[2:3],"\n\tδ=10% and α=50% and 90%: ",cumI_diffs[5:6])
    println("R0=3.   Deaths with detection and no CT:\n\tδ=0%, 5% and α=50% and 90%: ",[string(cumDeaths[i])*"⨦"*string(cumDeaths_std[i])  for i=1:3],"\n\tδ=10% and α=0%, 50% and 90%: ",[string(cumDeaths[i])*"⨦"*string(cumDeaths_std[i])  for i=4:6])

######## Detection AND Contact tracing
#R0=2.5
folders=["./contacts/session30_50sims/","./contacts/session31_50sims/"];plot_folder="./contacts/session30_50sims/plots_sessions30-31/";reference_files=[1,10]
    cumI_diffs,cumDeaths,cumDeaths_std=make_plots_1R0_twoδs_Kenya(folders,plot_folder,S0,δs,reference_files)
    println("R0=2.5. With detection and contact tracing, we expect the best case scenarios in cumulative deaths for R0=2.5 and both δ=5% and 10% to be ",cumDeaths[9]," and ",cumDeaths[18],", corresponding to a 90% detection rate and 90 days of contact tracing.")

#R0=3
folders=["./contacts/session32_50sims/","./contacts/session33_50sims/"];plot_folder="./contacts/session32_50sims/plots_sessions32-33/";reference_files=[1,10]
    cumI_diffs,cumDeaths,cumDeaths_std=make_plots_1R0_twoδs_Kenya(folders,plot_folder,S0,δs,reference_files)
    println("R0=3.  With detection and contact tracing, we expect the best case scenarios in cumulative deaths for R0=3 and both δ=5% and 10% to be ",cumDeaths[9]," and ",cumDeaths[18],", corresponding to a 90% detection rate and 90 days of contact tracing.")
    #=println("Scenario parameters\t\t\tGain compared to no intervention\t\t\tGain compared to detection and no tracing")
    println("R0 =2.5	δ=30%	α=50%\t\t\t",
            R0 =2.5	δ=30%	α=90%\t\t\t
            R0 =2.5	δ=70%	α=50%\t\t\t
            R0 =2.5	δ=70%	α=90%\t\t\t
            R0 =3	δ=30%	α=50%\t\t\t
            R0 =3	δ=30%	α=90%\t\t\t
            R0 =3	δ=70%	α=50%\t\t\t
            R0 =3	δ=70%	α=90%\t\t\t")=#

plot_folder="./contacts/session30_50sims/Icurves30-33/"
    CT_comparedtoNoIntervention_area("./contacts/session30_50sims/sims_sc3000.jld2","./contacts/session30_50sims/sims_sc3093.jld2",plot_folder,false)
    CT_comparedtoNoIntervention_area("./contacts/session31_50sims/sims_sc3100.jld2","./contacts/session31_50sims/sims_sc3193.jld2",plot_folder,false)
    CT_comparedtoNoIntervention_area("./contacts/session32_50sims/sims_sc3200.jld2","./contacts/session32_50sims/sims_sc3293.jld2",plot_folder,false)
    CT_comparedtoNoIntervention_area("./contacts/session33_50sims/sims_sc3300.jld2","./contacts/session33_50sims/sims_sc3393.jld2",plot_folder,true)

###### Ages with intervention
data_files=["./contacts/session30_50sims/sims_sc3000.jld2","./contacts/session30_50sims/sims_sc3053.jld2","./contacts/session31_50sims/sims_sc3153.jld2"]
plot_folder="./contacts/session30_50sims/plots_sessions20-21_det5_90days_50sims/"
#cumIfinalKilifiA=kilifi_ages_α5_90days(data_files,plot_folder,S0)

###### CumIs 2δs
folders=["./contacts/session30_50sims/","./contacts/session31_50sims/"];plot_folder="./contacts/session30_50sims/plots_sessions30-31/";
make_plots_twoδs(folders,plot_folder,S0,δs)
#CT_comparedtoNoIntervention("./contacts/results_session20_50sims/sims_sc2090.jld2","./contacts/results_session20_50sims/sims_sc2094.jld2",plot_folder)

######
data_files=["./contacts/session30_50sims/sims_sc3000.jld2","./contacts/session31_50sims/sims_sc3100.jld2","./contacts/session32_50sims/sims_sc3200.jld2","./contacts/session33_50sims/sims_sc3300.jld2"]
plot_folder="./contacts/session30_50sims/plots_sessions30-33_det0_50sims/"
cumI_medians_kenya,cumDeaths_stds_kenya,cumI_medians,cumDeaths_medians=kenya_barplot2R₀2δ(data_files,plot_folder,S0)

#######
data_files=["./contacts/results_session20_50sims/sims_sc2000.jld2","./contacts/results_session20_50sims/sims_sc2053.jld2","./contacts/results_session20_50sims/sims_sc2093.jld2","./contacts/results_session21_50sims/sims_sc2153.jld2","./contacts/results_session21_50sims/sims_sc2193.jld2"]
plot_folder="./contacts/results_session20_50sims/results_sessions20to21_det5_90days_50sims/"

kilifi_ages_α5_9_90days(data_files,plot_folder,S0)
#kilifi_ages_α5_9_90days_over60s(data_files,plot_folder,S0)

######## I curves compared to no intervention
plot_folder="./contacts/results_session20_50sims/oneExample20-23/"
#=oneCurve_CT_comparedtoNoIntervention("./contacts/results_session20_50sims/sims_sc2000.jld2","./contacts/results_session20_50sims/sims_sc2093.jld2",plot_folder)
oneCurve_CT_comparedtoNoIntervention("./contacts/results_session21_50sims/sims_sc2100.jld2","./contacts/results_session21_50sims/sims_sc2193.jld2",plot_folder)
oneCurve_CT_comparedtoNoIntervention("./contacts/results_session22_50sims/sims_sc2200.jld2","./contacts/results_session22_50sims/sims_sc2293.jld2",plot_folder)
oneCurve_CT_comparedtoNoIntervention("./contacts/results_session23_50sims/sims_sc2300.jld2","./contacts/results_session23_50sims/sims_sc2393.jld2",plot_folder)
=#
plot_folder="./contacts/results_session20_50sims/Icurves20-23/"
#=CT_comparedtoNoIntervention("./contacts/results_session20_50sims/sims_sc2000.jld2","./contacts/results_session20_50sims/sims_sc2093.jld2",plot_folder)
CT_comparedtoNoIntervention("./contacts/results_session21_50sims/sims_sc2100.jld2","./contacts/results_session21_50sims/sims_sc2193.jld2",plot_folder)
CT_comparedtoNoIntervention("./contacts/results_session22_50sims/sims_sc2200.jld2","./contacts/results_session22_50sims/sims_sc2293.jld2",plot_folder)
CT_comparedtoNoIntervention("./contacts/results_session23_50sims/sims_sc2300.jld2","./contacts/results_session23_50sims/sims_sc2393.jld2",plot_folder)
=#

#=
plot([sims_vector[10][2][t][12][1] for t=1:366])
plot([sims_vector[10][1][t][12][1] for t=1:366])

findfirst(x->x>=5,[sims_vector[10][1][t][12][1] for t=1:366])
findfirst(x->x==maximum([sims_vector[10][2][t][12][1] for t=1:366]),[sims_vector[10][2][t][12][1] for t=1:366])
findfirst(x->x>0,[sims_vector[10][2][t][12][1] for t=1:366])
=#

##########  Verify curves
folders=["./contacts/session30_50sims/","./contacts/session31_50sims/","./contacts/session32_50sims/","./contacts/session33_50sims/"];
plot_folder="./contacts/session30_50sims/testcurves/";
curves(folders,plot_folder)
