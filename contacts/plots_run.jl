include("plots_v10.jl")
    push!(LOAD_PATH, "./src")
    push!(LOAD_PATH, "./contacts")
    include("contacttracing_funtionswithdeaths.jl")

sc=30;sims=string(50);
session1=string(sc);session2=string(sc+1);session3=string(sc+2);session4=string(sc+3);

####### Kenya per region and with No Intervention
data_files=["./contacts/session"*session1*"_"*sims*"sims/sims_sc"*session1*"00.jld2","./contacts/session"*session2*"_"*sims*"sims/sims_sc"*session2*"00.jld2","./contacts/session"*session3*"_"*sims*"sims/sims_sc"*session3*"00.jld2","./contacts/session"*session4*"_"*sims*"sims/sims_sc"*session4*"00.jld2"]
    #data_files=["./contacts/session"*session1*"_"*sims*"sims/sims_sc"*session1*"00.jld2","./contacts/session"*session1*"_"*sims*"sims/sims_sc"*session1*"01.jld2","./contacts/session"*session1*"_"*sims*"sims/sims_sc"*session1*"02.jld2","./contacts/session"*session1*"_"*sims*"sims/sims_sc"*session1*"03.jld2"]  #For session 42: no intervention with and without death rate
    plot_folder="./contacts/session"*session1*"_"*sims*"sims/plots_sessions1-4_NoInterv_counties/"
    cumI_medians_kenya,cumI_stds_kenya,cumI_medians,cumDeaths_medians_kenya,cumDeaths_stds_kenya,cumDeaths_medians=cumIandDeaths_kenya(data_files,plot_folder,S0)
    println("\nCases in all Kenya for R0=2.5 and R0=3, with sympt rate δ=5% and δ=10% is: ")
    for i=1:4 println(cumI_medians_kenya[i],"±",cumI_stds_kenya[i], " (",100*cumI_medians_kenya[i]/sum(S0),"% of the total Kenyan population)") end
    println("\nDeaths in all Kenya for R0=2.5 and R0=3, with sympt rate δ=5% and δ=10% is: ")
    for i=1:4 println(cumDeaths_medians_kenya[i],"±",cumDeaths_stds_kenya[i], " (",100*cumDeaths_medians_kenya[i]/sum(S0),"% of the total Kenyan population)") end

####### Kenya per ages with No Intervention
data_files=["./contacts/session"*session1*"_"*sims*"sims/sims_sc"*session1*"00.jld2","./contacts/session"*session2*"_"*sims*"sims/sims_sc"*session2*"00.jld2","./contacts/session"*session3*"_"*sims*"sims/sims_sc"*session3*"00.jld2","./contacts/session"*session4*"_"*sims*"sims/sims_sc"*session4*"00.jld2"]
    plot_folder="./contacts/session"*session1*"_"*sims*"sims/plots_sessions1-4_NoInterv_age/"
    cumIsA_2,cumIsA_4,cumDeaths_1,cumDeaths_2,cumDeaths_3,cumDeaths_4=cumIandDeaths_ages_kenya(data_files,plot_folder,S0)
    cumIsA_over60s_2=[];cumIsA_over60s_4=[];
    for sim=1:parse(Int64, sims)
        I=0;    for a=13:16     I+=cumIsA_2[a][sim];    end;        push!(cumIsA_over60s_2,I)
        I=0;    for a=13:16     I+=cumIsA_4[a][sim];    end;        push!(cumIsA_over60s_4,I)
    end
    cumDeathsA_over60s_1=[];cumDeathsA_over60s_2=[];cumDeathsA_over60s_3=[];cumDeathsA_over60s_4=[]
    for sim=1:parse(Int64, sims)
        d=0;    for a=13:16     d+=cumDeaths_1[a][sim];    end;        push!(cumDeathsA_over60s_1,d)
        d=0;    for a=13:16     d+=cumDeaths_2[a][sim];    end;        push!(cumDeathsA_over60s_2,d)
        d=0;    for a=13:16     d+=cumDeaths_3[a][sim];    end;        push!(cumDeathsA_over60s_3,d)
        d=0;    for a=13:16     d+=cumDeaths_4[a][sim];    end;        push!(cumDeathsA_over60s_4,d)
    end
    println("For R0=2.5 and 3, and δ=10.%, in the over 60s we estimate ",median(cumIsA_over60s_2),"±",std(cumIsA_over60s_2)," (",100*median(cumIsA_over60s_2)/sum(S0perAge[13:end]),"%) and ",median(cumIsA_over60s_4),"±",std(cumIsA_over60s_4)," (",100*median(cumIsA_over60s_4)/sum(S0perAge[13:end]),"%).")
    println("As for fatalities, deaths in the over 60s is estimated as: \n\tfor R0=2.5 and δ=5% ",median(cumDeathsA_over60s_1),"±",std(cumDeathsA_over60s_1),
                                                                        "\n\tfor R0=2.5 and δ=10% ",median(cumDeathsA_over60s_2),"±",std(cumDeathsA_over60s_2),
                                                                        "\n\tfor R0=3 and δ=5% ",median(cumDeathsA_over60s_3),"±",std(cumDeathsA_over60s_3),
                                                                        "\n\tfor R0=10 and δ=10% ",median(cumDeathsA_over60s_4),"±",std(cumDeathsA_over60s_4))
    cumIsA_20to35_2=[];cumIsA_20to35_4=[];
    for sim=1:parse(Int64, sims)
        I=0;    for a=5:7     I+=cumIsA_2[a][sim];    end;        push!(cumIsA_20to35_2,I)
        I=0;    for a=5:7     I+=cumIsA_4[a][sim];    end;        push!(cumIsA_20to35_4,I)
    end
    cumDeathsA_20to35_1=[];cumDeathsA_20to35_2=[];cumDeathsA_20to35_3=[];cumDeathsA_20to35_4=[]
    for sim=1:parse(Int64, sims)
        d=0;    for a=5:7     d+=cumDeaths_1[a][sim];    end;        push!(cumDeathsA_20to35_1,d)
        d=0;    for a=5:7     d+=cumDeaths_2[a][sim];    end;        push!(cumDeathsA_20to35_2,d)
        d=0;    for a=5:7     d+=cumDeaths_3[a][sim];    end;        push!(cumDeathsA_20to35_3,d)
        d=0;    for a=5:7     d+=cumDeaths_4[a][sim];    end;        push!(cumDeathsA_20to35_4,d)
    end
    println("\nFor R0=2.5 and 3, and δ=10.%, cases between 20to35yrs, we estimate ",median(cumIsA_20to35_2),"±",std(cumIsA_20to35_2)," (",100*median(cumIsA_20to35_2)/sum(S0perAge[5:7]),"%) and ",median(cumIsA_20to35_4),"±",std(cumIsA_20to35_4)," (",100*median(cumIsA_20to35_4)/sum(S0perAge[5:7]),"%).")
    println("As for fatalities, deaths between 20to35yrs is estimated as: \n\tfor R0=2.5 and δ=5% ",median(cumDeathsA_20to35_1),"±",std(cumDeathsA_20to35_1),
                                                                        "\n\tfor R0=2.5 and δ=10% ",median(cumDeathsA_20to35_2),"±",std(cumDeathsA_20to35_2),
                                                                        "\n\tfor R0=3 and δ=5% ",median(cumDeathsA_20to35_3),"±",std(cumDeathsA_20to35_3),
                                                                        "\n\tfor R0=10 and δ=10% ",median(cumDeathsA_20to35_4),"±",std(cumDeathsA_20to35_4))

######## Detection without contact tracing (CT0 gain)
#R0=2.5
folders=["./contacts/session"*session1*"_"*sims*"sims/CT0/","./contacts/session"*session2*"_"*sims*"sims/CT0/"];plot_folder="./contacts/session"*session1*"_"*sims*"sims/CT0/";reference_files=[1,4]
    cumI_diffs,cumDeaths,cumDeaths_std=make_plots_1R0_twoδs_Kenya(folders,plot_folder,S0,δs,reference_files)
    println("R0=2.5. Gain in cases with detection and no CT:\n\tδ=5% and α=50% and 90%: ",cumI_diffs[2:3],"\n\tδ=10% and α=50% and 90%: ",cumI_diffs[5:6])
    println("R0=2.5.   Deaths with detection and no CT:\n\tδ=0%, 5% and α=50% and 90%: ",[string(cumDeaths[i])*"±"*string(cumDeaths_std[i])  for i=1:3],"\n\tδ=10% and α=0%, 50% and 90%: ",[string(cumDeaths[i])*"±"*string(cumDeaths_std[i])  for i=4:6])

#R0=3
folders=["./contacts/session"*session3*"_"*sims*"sims/CT0/","./contacts/session"*session4*"_"*sims*"sims/CT0/"];plot_folder="./contacts/session"*session3*"_"*sims*"sims/CT0/";reference_files=[1,4]
    cumI_diffs,cumDeaths,cumDeaths_std=make_plots_1R0_twoδs_Kenya(folders,plot_folder,S0,δs,reference_files)
    println("R0=3.   Gain in cases with detection and no CT:\n\tδ=5% and α=50% and 90%: ",cumI_diffs[2:3],"\n\tδ=10% and α=50% and 90%: ",cumI_diffs[5:6])
    println("R0=3.   Deaths with detection and no CT:\n\tδ=0%, 5% and α=50% and 90%: ",[string(cumDeaths[i])*"±"*string(cumDeaths_std[i])  for i=1:3],"\n\tδ=10% and α=0%, 50% and 90%: ",[string(cumDeaths[i])*"±"*string(cumDeaths_std[i])  for i=4:6])

######## Detection AND Contact tracing in all Kenya
#R0=2.5
folders=["./contacts/session"*session1*"_"*sims*"sims/","./contacts/session"*session2*"_"*sims*"sims/"];plot_folder="./contacts/session"*session1*"_"*sims*"sims/plots_sessions1-2/";reference_files=[1,10]
    cumI_diffs,cumDeaths,cumDeaths_std=make_plots_1R0_twoδs_Kenya(folders,plot_folder,S0,δs,reference_files)
    println("R0=2.5. With detection and contact tracing, we expect the best case scenarios in cumulative deaths for R0=2.5 and both δ=5% and 10% to be ",cumDeaths[9]," and ",cumDeaths[18],", corresponding to a 90% detection rate and 90 days of contact tracing.")

#R0=3
folders=["./contacts/session"*session3*"_"*sims*"sims/","./contacts/session"*session4*"_"*sims*"sims/"];plot_folder="./contacts/session"*session3*"_"*sims*"sims/plots_sessions3-4/";reference_files=[1,10]
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

plot_folder="./contacts/session"*session1*"_"*sims*"sims/Icurves/"
    IKenya_nointervention_median1,IKenya_median1=CT_comparedtoNoIntervention_area("./contacts/session"*session1*"_"*sims*"sims/sims_sc"*session1*"00.jld2","./contacts/session"*session1*"_"*sims*"sims/sims_sc"*session1*"99.jld2",plot_folder,false)
    IKenya_nointervention_median2,IKenya_median2=CT_comparedtoNoIntervention_area("./contacts/session"*session2*"_"*sims*"sims/sims_sc"*session2*"00.jld2","./contacts/session"*session2*"_"*sims*"sims/sims_sc"*session2*"99.jld2",plot_folder,false)
    IKenya_nointervention_median3,IKenya_median3=CT_comparedtoNoIntervention_area("./contacts/session"*session3*"_"*sims*"sims/sims_sc"*session3*"00.jld2","./contacts/session"*session3*"_"*sims*"sims/sims_sc"*session3*"99.jld2",plot_folder,false)
    IKenya_nointervention_median4,IKenya_median4=CT_comparedtoNoIntervention_area("./contacts/session"*session4*"_"*sims*"sims/sims_sc"*session4*"00.jld2","./contacts/session"*session4*"_"*sims*"sims/sims_sc"*session4*"99.jld2",plot_folder,true)

# Peak difference
findall(x->x==maximum(IKenya_median1),IKenya_median1)[1]-findall(x->x==maximum(IKenya_nointervention_median1),IKenya_nointervention_median1)[1]
findall(x->x==maximum(IKenya_median2),IKenya_median2)[1]-findall(x->x==maximum(IKenya_nointervention_median2),IKenya_nointervention_median2)[1]
findall(x->x==maximum(IKenya_median3),IKenya_median3)[1]-findall(x->x==maximum(IKenya_nointervention_median3),IKenya_nointervention_median3)[1]
findall(x->x==maximum(IKenya_median4),IKenya_median4)[1]-findall(x->x==maximum(IKenya_nointervention_median4),IKenya_nointervention_median4)[1]
###### Ages with intervention
data_files=["./contacts/session"*session1*"_"*sims*"sims/sims_sc"*session1*"00.jld2","./contacts/session"*session1*"_"*sims*"sims/sims_sc"*session1*"53.jld2","./contacts/session"*session2*"_"*sims*"sims/sims_sc"*session2*"53.jld2"]
plot_folder="./contacts/session"*session1*"_"*sims*"sims/plots_sessions20-21_det5_90days_"*sims*"sims/"
#cumIfinalKilifiA=kilifi_ages_α5_90days(data_files,plot_folder,S0)

###### CumIs 2δs
folders=["./contacts/session"*session1*"_"*sims*"sims/","./contacts/session"*session2*"_"*sims*"sims/"];plot_folder="./contacts/session"*session1*"_"*sims*"sims/plots_sessions30-31/";
make_plots_twoδs(folders,plot_folder,S0,δs)
#CT_comparedtoNoIntervention("./contacts/results_session20_"*sims*"sims/sims_sc2090.jld2","./contacts/results_session20_"*sims*"sims/sims_sc2094.jld2",plot_folder)

######
data_files=["./contacts/session"*session1*"_"*sims*"sims/sims_sc"*session1*"00.jld2","./contacts/session"*session2*"_"*sims*"sims/sims_sc"*session2*"00.jld2","./contacts/session"*session3*"_"*sims*"sims/sims_sc"*session3*"00.jld2","./contacts/session"*session4*"_"*sims*"sims/sims_sc"*session4*"00.jld2"]
plot_folder="./contacts/session"*session1*"_"*sims*"sims/plots_det0_"*sims*"sims/"
cumI_medians_kenya,cumDeaths_stds_kenya,cumI_medians,cumDeaths_medians=kenya_barplot2R₀2δ(data_files,plot_folder,S0)

#######
data_files=["./contacts/results_session20_"*sims*"sims/sims_sc2000.jld2","./contacts/results_session20_"*sims*"sims/sims_sc2053.jld2","./contacts/results_session20_"*sims*"sims/sims_sc2093.jld2","./contacts/results_session21_"*sims*"sims/sims_sc2153.jld2","./contacts/results_session21_"*sims*"sims/sims_sc2193.jld2"]
plot_folder="./contacts/results_session20_"*sims*"sims/results_det5_90days_"*sims*"sims/"

kilifi_ages_α5_9_90days(data_files,plot_folder,S0)
#kilifi_ages_α5_9_90days_over60s(data_files,plot_folder,S0)

######## I curves compared to no intervention
plot_folder="./contacts/results_session"*session1*"_"*sims*"sims/oneExample/"
#=oneCurve_CT_comparedtoNoIntervention("./contacts/results_session20_"*sims*"sims/sims_sc2000.jld2","./contacts/results_session20_"*sims*"sims/sims_sc2093.jld2",plot_folder)
oneCurve_CT_comparedtoNoIntervention("./contacts/results_session21_"*sims*"sims/sims_sc2100.jld2","./contacts/results_session21_"*sims*"sims/sims_sc2193.jld2",plot_folder)
oneCurve_CT_comparedtoNoIntervention("./contacts/results_session22_"*sims*"sims/sims_sc2200.jld2","./contacts/results_session22_"*sims*"sims/sims_sc2293.jld2",plot_folder)
oneCurve_CT_comparedtoNoIntervention("./contacts/results_session23_"*sims*"sims/sims_sc2300.jld2","./contacts/results_session23_"*sims*"sims/sims_sc2393.jld2",plot_folder)
=#
plot_folder="./contacts/results_session"*session1*"_"*sims*"sims/Icurves/"
#=CT_comparedtoNoIntervention("./contacts/results_session20_"*sims*"sims/sims_sc2000.jld2","./contacts/results_session20_"*sims*"sims/sims_sc2093.jld2",plot_folder)
CT_comparedtoNoIntervention("./contacts/results_session21_"*sims*"sims/sims_sc2100.jld2","./contacts/results_session21_"*sims*"sims/sims_sc2193.jld2",plot_folder)
CT_comparedtoNoIntervention("./contacts/results_session22_"*sims*"sims/sims_sc2200.jld2","./contacts/results_session22_"*sims*"sims/sims_sc2293.jld2",plot_folder)
CT_comparedtoNoIntervention("./contacts/results_session23_"*sims*"sims/sims_sc2300.jld2","./contacts/results_session23_"*sims*"sims/sims_sc2393.jld2",plot_folder)
=#

#plot_folder="./contacts/session"*session1*"_"*sims*"sims/IcurvesNoIntervention/"
#CT_comparedtoNoIntervention("./contacts/session"*session1*"_"*sims*"sims/sims_sc"*session1*"00.jld2","./contacts/session"*session2*"_"*sims*"sims/sims_sc"*session2*"00.jld2",plot_folder)

#=
plot([sims_vector[10][2][t][12][1] for t=1:366])
plot([sims_vector[10][1][t][12][1] for t=1:366])

findfirst(x->x>=5,[sims_vector[10][1][t][12][1] for t=1:366])
findfirst(x->x==maximum([sims_vector[10][2][t][12][1] for t=1:366]),[sims_vector[10][2][t][12][1] for t=1:366])
findfirst(x->x>0,[sims_vector[10][2][t][12][1] for t=1:366])
=#

##########  Verify curves
folders=["./contacts/session"*session1*"_"*sims*"sims/","./contacts/session"*session2*"_"*sims*"sims/"]
    plot_folder="./contacts/session"*session1*"_"*sims*"sims/testcurves/";
    curves(folders,plot_folder)

folders=["./contacts/session"*session3*"_"*sims*"sims/","./contacts/session"*session4*"_"*sims*"sims/"]
    plot_folder="./contacts/session"*session3*"_"*sims*"sims/testcurves/";
    curves(folders,plot_folder)

####### Kenya per region comparing different interventions (no interv/Nairobi only/all kenya/)
#R0=2.5 δ=10%
folder="./contacts/sessionR025delta10/";#data_files=readdir(folder);data_files=[folder*data_files[i] for i=1:size(data_files,1) if endswith(data_files[i],".jld2")];
    data_files=["./contacts/sessionR025delta10/sims_sc3100.jld2", "./contacts/sessionR025delta10/sims_sc3590.jld2", "./contacts/sessionR025delta10/sims_sc3599.jld2", "./contacts/sessionR025delta10/sims_sc3990.jld2", "./contacts/sessionR025delta10/sims_sc3999.jld2", "./contacts/sessionR025delta10/sims_sc3190.jld2", "./contacts/sessionR025delta10/sims_sc3199.jld2"]
    plot_folder=folder*"output/";
    cumI_medians_counties,gainI,cumDeaths_medians_counties,gainDeaths=cumIandDeaths_kenya_COMPAREINTERVENTIONS_counties(data_files,plot_folder,"with R0=2.5 and sympt d=10%")
    for i=1:size(gainI,1)  println(gainI[i][1:10])    end
    for i=1:size(gainDeaths,1)  println(gainDeaths[i][1:10])    end

#R0=3 δ=10%
folder="./contacts/sessionR03delta10/";
    data_files=["./contacts/sessionR03delta10/sims_sc3300.jld2", "./contacts/sessionR03delta10/sims_sc3790.jld2", "./contacts/sessionR03delta10/sims_sc3799.jld2", "./contacts/sessionR03delta10/sims_sc4190.jld2", "./contacts/sessionR03delta10/sims_sc4199.jld2", "./contacts/sessionR03delta10/sims_sc3390.jld2", "./contacts/sessionR03delta10/sims_sc3399.jld2"]
    plot_folder=folder*"output/";
    cumI_medians_counties,gainI,cumDeaths_medians_counties,gainDeaths=cumIandDeaths_kenya_COMPAREINTERVENTIONS_counties(data_files,plot_folder,"with R0=3 and sympt d=10%")
    for i=1:size(gainI,1)  println(gainI[i][1:10])    end
    for i=1:size(gainDeaths,1)  println(gainDeaths[i][1:10])    end
