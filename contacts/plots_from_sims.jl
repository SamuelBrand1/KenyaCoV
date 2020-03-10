using Plots,MAT

results_folder="./contacts\\results_session10\\"
taup_list=[0.0,0.25,0.5,0.75,0.9]
data_files=readdir(results_folder)#["sims100_taup0.0_capacity10000.mat","sims100_taup0.25_capacity10000.mat","sims100_taup0.5_capacity10000.mat","sims100_taup0.75_capacity10000.mat","sims100_taup0.9_capacity10000.mat"];
data_files=[s for s in data_files if endswith(s,".mat")]
n_traj=10

## Plotting functions IN NAIROBI
function ALL_TAUP_cumulative_symptomatics_sumages(taup_list,results_folder,data_files)
    p=plot(); #title("Final cumulative Infecteds(ID+IQ) 100 sims (sum all ages)");   hold on
    #Legend=cell(size(taup_list,2),1);
    for i=1:size(taup_list,1)
        file = matopen(results_folder*data_files[i])
        data_per_taup=read(file, "sims")
        close(file)
        data_per_taup=data_per_taup[:,2:end]

        Id_per_taup_nairobi=[]   # Infecteds=I_D+I_Q
        for j=1:size(data_per_taup,2)
            #Id_per_taup_nairobi=[Id_per_taup_nairobi sum(data_per_taup[end,j][4,:,8])]
            push!(Id_per_taup_nairobi,sum(data_per_taup[end,j][4,:,8]))
        end
        #Legend{i}=strcat("taup=",num2str(taup_list(i)));
        plot!(p,Id_per_taup_nairobi, label="taup"*string(taup_list[i]))
    end
    #hold off
    #legend(Legend);
    display(p)
    savefig(p,results_folder*"cumulative_symptomatics_allages.png")
end
#=
function ALL_TAUP_cumulative_symptomatics_sumages_BOXPLOT(taup_list,results_folder,data_files)
    p=figure(); title("Final cumulative Infecteds(ID+IQ) 100 sims (sum all ages)");   hold on
    Legend=cell(size(taup_list,2),1);
    data=[];
    for i=1:size(taup_list,2)
        data_per_taup=load(data_files(i));
        data_per_taup=data_per_taup.sims(:,2:end);

        Id_per_taup_nairobi=[];   % Infecteds=I_D+I_Q
        for j=1:size(data_per_taup,2)
            Id_per_taup_nairobi=[Id_per_taup_nairobi;  sum(data_per_taup{end,j}(4,:,8))];
        end
        Legend{i}=strcat("taup=",num2str(taup_list(i)));
        data=[data Id_per_taup_nairobi];
    end
    boxplot(data);
    hold off
    legend(findobj(gca,'Tag','Box'),Legend);
    saveas(p,results_folder+"cumulative_symptomatics_allages_BOXPLOT.png");
end
=#
## Calling the functions , IN NAIROBI
ALL_TAUP_cumulative_symptomatics_sumages(taup_list,results_folder,data_files)
#ALL_TAUP_cumulative_symptomatics_sumages_BOXPLOT(taup_list,results_folder,data_files)
