addpath '.\';
addpath '.\results_session3\';
results_folder='.\results_session3\';
taup_list=[0.0,0.25,0.5,0.75,0.9];
data_files=["sims100_taup0.0_capacity10000.mat","sims100_taup0.25_capacity10000.mat","sims100_taup0.5_capacity10000.mat","sims100_taup0.75_capacity10000.mat","sims100_taup0.9_capacity10000.mat"];
n_traj=100;



% data_per_taup=load(data_files(1)); % for taup=taup_list(1)
% time=data_per_taup.sims(:,1);
% data_per_taup=data_per_taup.sims(:,2:end);

%% calling the functions , IN NAIROBI
% cumulative_symptomatics_perage(0.0,data_per_taup,results_folder);
% cumulative_symptomatics_sumallages(0.0,data_per_taup,results_folder);

%ALL_TAUP_cumulative_symptomatics_perage(taup_list,data_per_taup,results_folder,data_files);
%ALL_TAUP_cumulative_symptomatics_sumages(taup_list,results_folder,data_files);
data=ALL_TAUP_cumulative_symptomatics_sumages_BOXPLOT(taup_list,results_folder,data_files);

%% plotting functions IN NAIROBI
function cumulative_symptomatics_perage(taup,data_per_taup,results_folder)
    p=figure(); title("Final cumulative Infecteds(ID+IQ) 100 sims (per age) - taup="+num2str(taup));   hold on
    Legend=cell(16,1);
    for a=1:16
        Id_per_taup_nairobi_age=[];   % Infecteds=I_D+I_Q
        for j=1:size(data_per_taup,2)
            Id_per_taup_nairobi_age=[Id_per_taup_nairobi_age; sum(data_per_taup{end,j}(4,a,8))];
            %Id_per_taup_nairobi=[Id_per_taup_nairobi; cell(4,:,8)];
        end
        Legend{a}=strcat("a=",num2str(a));
        plot(Id_per_taup_nairobi_age); %legend("a"+a);
    end
    legend(Legend);
    
    saveas(p,results_folder+"cumulative_symptomatics_perage_"+num2str(taup)+".png");
    hold off
end
function cumulative_symptomatics_sumallages(taup,data_per_taup,results_folder)
    p=figure(); title("Final cumulative Infecteds(ID+IQ) 100 sims (sum all ages) - taup="+num2str(taup));
    Id_per_taup_nairobi=[];   % Infecteds=I_D+I_Q
    for j=1:size(data_per_taup,2)
        Id_per_taup_nairobi=[Id_per_taup_nairobi; sum(data_per_taup{end,j}(4,:,8))];
    end
    plot(Id_per_taup_nairobi);
    saveas(p,results_folder+"cumulative_symptomatics_allages_"+num2str(taup)+".png");
end
function ALL_TAUP_cumulative_symptomatics_perage(taup_list,results_folder,data_files)
    for i=1:size(taup_list,2)
        data_per_taup=load(data_files(i));
        data_per_taup=data_per_taup.sims(:,2:end);
        cumulative_symptomatics_perage(taup_list(i),data_per_taup,results_folder);
    end
end
function ALL_TAUP_cumulative_symptomatics_sumages(taup_list,results_folder,data_files)
    p=figure(); title("Final cumulative Infecteds(ID+IQ) 100 sims (sum all ages)");   hold on
    Legend=cell(size(taup_list,2),1);
    for i=1:size(taup_list,2)
        data_per_taup=load(data_files(i));
        data_per_taup=data_per_taup.sims(:,2:end);
        
        Id_per_taup_nairobi=[];   % Infecteds=I_D+I_Q
        for j=1:size(data_per_taup,2)
            Id_per_taup_nairobi=[Id_per_taup_nairobi; sum(data_per_taup{end,j}(4,:,8))];
        end
        Legend{i}=strcat("taup=",num2str(taup_list(i)));
        plot(Id_per_taup_nairobi);
    end
    hold off
    legend(Legend);
    saveas(p,results_folder+"cumulative_symptomatics_allages.png");%_"+num2str(taup)+".png");
end
function data=ALL_TAUP_cumulative_symptomatics_sumages_BOXPLOT(taup_list,results_folder,data_files)
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
    saveas(p,results_folder+"cumulative_symptomatics_allages_BOXPLOT.png");%_"+num2str(taup)+".png");
end