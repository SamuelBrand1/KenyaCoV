using MAT, Plots,StatsPlots

##Plotting from mat files containing the final cumulative sum of Infecteds in Nairobi
function boxplot_τₚ_list(mat_files,τₚ_list)
    bplt=boxplot(title="Cumulative infecteds (all ages) in Nairobi")
    for i=1:size(mat_files,1)
        file = matopen(mat_files[i])
        finalCumINairobiSum=read(file, "finalCumINairobiSum") # note that this does NOT introduce a variable ``varname`` into scope
        close(file)
        boxplot!(bplt,finalCumINairobiSum, x=1,label="max_capacity=1e3,tau_p="*string(τₚ_list[i]))
    end
    display(bplt)
    savefig("CumulativeInfectedsNairobi1e3.png")
end

boxplot_τₚ_list(["finalCumIData_sims1000_taup0.0.mat","finalCumIData_sims1000_taup0.25.mat",
                "finalCumIData_sims1000_taup0.5.mat","finalCumIData_sims1000_taup0.75.mat",
                "finalCumIData_sims1000_taup0.9.mat"],
                [0.0,0.25,0.5,0.75,0.9])
