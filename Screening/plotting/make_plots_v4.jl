include("plotting_functions_v2.jl")

function do_plots(I_folders,titles)
    W_folder="W:/BACKUP_KenyaCoV/2020-08-07_KenyaCoV/Screening/";sub1="output_rand_K_500tests/";sub2="output_rand_K_1000tests/";sims=string(1000)
    for i = 1:size(I_folders,1)
        scenarios=[(W_folder*sub2*"I0_NoInterv/session1_"*sims*"sims/sc1.jld2",0,1,0),
            (W_folder*sub1*I_folders[i]*"session1_"*sims*"sims/sc5.jld2",1,5,500), (W_folder*sub1*I_folders[i]*"session2_"*sims*"sims/sc5.jld2",2,5,500), (W_folder*sub1*I_folders[i]*"session3_"*sims*"sims/sc5.jld2",3,5,500),
            (W_folder*sub2*I_folders[i]*"session1_"*sims*"sims/sc5.jld2",1,5,1000),(W_folder*sub2*I_folders[i]*"session2_"*sims*"sims/sc5.jld2",2,5,1000),(W_folder*sub2*I_folders[i]*"session3_"*"500"*"sims/sc5.jld2",3,5,1000)]
        rez=make_plots_bestOfs(scenarios,"Screening/plotting/plots__2020_"*Dates.format(now(), "mm_dd__HH_MM")*"__fixed2_"*I_folders[i],titles[i])
    end
    for i = 1:size(I_folders,1)
        scenarios=[(W_folder*sub2*"I0_NoInterv/session1_"*sims*"sims/sc1.jld2",0,1,0),
            #(W_folder*sub1*I_folders[i]*"session1_100sims/sc5.jld2",1,5,500), (W_folder*sub1*I_folders[i]*"session2_100sims/sc5.jld2",2,5,500), (W_folder*sub1*I_folders[i]*"session3_100sims/sc5.jld2",3,5,500),
            (W_folder*sub2*I_folders[i]*"session1_"*sims*"sims/sc5.jld2",1,5,1000),(W_folder*sub2*I_folders[i]*"session2_"*sims*"sims/sc5.jld2",2,5,1000),(W_folder*sub2*I_folders[i]*"session3_"*"500"*"sims/sc5.jld2",3,5,1000)]
        rez=make_plots_bestOfs(scenarios,"Screening/plotting/plots__2020_"*Dates.format(now(), "mm_dd__HH_MM")*"__fixed2_"*I_folders[i][1:end-1]*"_2/",titles[i])
    end
end

##
do_plots(["I1_CTH/"],
         ["Contact tracing of hospitalized"])


do_plots(["I2_SympS"],#,"I3_SympSCT/","I4_MS/","I5_MSCT/"],
         ["Screening of symptomatics"])#,"Screening of symptomatics+CT at K","Randomized screening","Randomized screening+CT at K"])

do_plots(["I3_SympSCT/"],
        ["Screening of symptomatics+CT"])

do_plots(["I4_MS/"],
        ["Randomized screening"])

do_plots(["I5_MSCT/"],
        ["Randomized screening+CT"])
