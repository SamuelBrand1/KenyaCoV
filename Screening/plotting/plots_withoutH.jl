push!(LOAD_PATH, "./src")
    push!(LOAD_PATH, "./Screening")
    using DifferentialEquations,JLD2, CSV, DataFrames,Dates,Statistics

####### get cumI and cumQ functions
function final_cumulatives_CSV(folders)
    cumIs=[];cumQs=[];files=[]
    cumIs_sd=[];cumQs_sd=[]
    for folder in folders
        data_files=readdir(folder)
        data_files=[s for s in data_files if endswith(s,".jld2")]
        for data_file in data_files
            @load folder*data_file  R
            push!(files,folder*data_file)
            push!(cumIs,median( [sum(R.u[sim][1][:,:])+sum(R.u[sim][2][:,:])+sum(R.u[sim][3][:,:])   for sim=1:size(R.u,1)]))
            push!(cumQs,median( [sum(R.u[sim][5][:,:])   for sim=1:size(R.u,1)]))
            push!(cumIs_sd,std( [sum(R.u[sim][1][:,:])+sum(R.u[sim][2][:,:])+sum(R.u[sim][3][:,:])   for sim=1:size(R.u,1)]))
            push!(cumQs_sd,std( [sum(R.u[sim][5][:,:])   for sim=1:size(R.u,1)]))
        end
    end
    df=DataFrame(files=files,cumIs=cumIs,cumIs_sd=cumIs_sd,cumQs=cumQs,cumQs_sd=cumQs_sd)
    CSV.write(folders[1]*"final_cumulatives__2020_"*Dates.format(now(), "mm_dd__HH_MM") *".csv",df)
    return cumIs,cumQs
end

######## final cumulatives to CSV
#with β=2.5
folders=["./Screening/output_withoutH/0_NoInterv/session2_1000sims/",
         "./Screening/output_withoutH/1_CTH/session2_1000sims/","./Screening/output_withoutH/1_CTH/session4_1000sims/",
         "./Screening/output_withoutH/2_SympS/session3_1000sims/","./Screening/output_withoutH/2_SympS/session4_1000sims/",
         "./Screening/output_withoutH/3_SympSCT/session3_1000sims/","./Screening/output_withoutH/3_SympSCT/session4_1000sims/","./Screening/output_withoutH/3_SympSCT/session7_1000sims/","./Screening/output_withoutH/3_SympSCT/session8_1000sims/",
         "./Screening/output_withoutH/4_MS/session3_1000sims/","./Screening/output_withoutH/4_MS/session4_1000sims/",
         "./Screening/output_withoutH/5_MSCT/session3_1000sims/","./Screening/output_withoutH/5_MSCT/session4_1000sims/","./Screening/output_withoutH/5_MSCT/session7_1000sims/","./Screening/output/5_MSCT/session8_1000sims/"
         ]
    cumIs,cumQs=final_cumulatives_CSV(folders)

############### Compare with No Interv
@load "./Screening/output/0_NoInterv/session2_2sims/noInterv_sc1_2sims.jld2" R
    R_0=deepcopy(R)
    cumI_0=median([sum(R_0.u[sim][1][:,:])+sum(R_0.u[sim][2][:,:])+sum(R_0.u[sim][3][:,:])   for sim=1:size(R_0.u,1)])

using Plots;p=plot();for sim=1:size(R_0.u,1)   plot!(p,R_0.u[sim][end][:],color=:red);end;display(p)

@load "./Screening/output/1_CTH/session2_2sims/intervCTH_sc1_2sims.jld2" R
    cumI=median( [sum(R.u[sim][1][:,:])+sum(R.u[sim][2][:,:])+sum(R.u[sim][3][:,:])   for sim=1:size(R.u,1)])
    cumI_0-cumI,#=cumQ==#median( [sum(R.u[sim][5][:,:])   for sim=1:size(R.u,1)])

#=using Plots;b=bar();@load "./Screening/6_Soc/session2_100sims/intervSoc_sc1_100sims.jld2" R#@load "./Screening/6_Lock/session2_1sims/lock_sc1_1000sims.jld2" R
    cumI=median( [sum(R.u[sim][1][:,:])+sum(R.u[sim][2][:,:])+sum(R.u[sim][3][:,:])   for sim=1:size(R.u,1)])
    bar!(b,[cumI_0,cumI],legend=false,ylims=(3e7,Inf));display(b)
    cumI_0-cumI,#=cumQ==#median( [sum(R.u[sim][5][:,:])   for sim=1:size(R.u,1)])=#

@load "./Screening/6_Soc/session2_100sims/intervSoc_sc1_100sims.jld2" R
    cumI=median( [sum(R.u[sim][1][:,:])+sum(R.u[sim][2][:,:])+sum(R.u[sim][3][:,:])   for sim=1:size(R.u,1)])
    cumI_0-cumI,#=cumQ==#median( [sum(R.u[sim][5][:,:])   for sim=1:size(R.u,1)])

@load "./Screening/7_CTHSoc/session2_100sims/intervCTHSoc_sc1_100sims.jld2" R
    cumI=median( [sum(R.u[sim][1][:,:])+sum(R.u[sim][2][:,:])+sum(R.u[sim][3][:,:])   for sim=1:size(R.u,1)])
    cumI_0-cumI,#=cumQ==#median( [sum(R.u[sim][5][:,:])   for sim=1:size(R.u,1)])

@load "./Screening/8_SympSSoc/session3_100sims/intervSympSSoc_sc15_100sims.jld2" R
    cumI=median( [sum(R.u[sim][1][:,:])+sum(R.u[sim][2][:,:])+sum(R.u[sim][3][:,:])   for sim=1:size(R.u,1)])
    cumI_0-cumI,#=cumQ==#median( [sum(R.u[sim][5][:,:])   for sim=1:size(R.u,1)])
p2=plot();for sim=1:size(R.u,1)   plot!(p2,R.u[sim][end][:],color=:blue,legend=false);end;display(p2)

using Plots;p=plot();for sim=1:size(R_0.u,1)   plot!(p,R_0.u[sim][end][:],color=:red);end;
                    for sim=1:size(R.u,1)   plot!(p,R.u[sim][end][:],color=:blue,legend=false);end;display(p)


## Bar plots for comparison

function final_cumulatives(files)
    cumIs=[];cumQs=[]#;files=[]
    cumIs_sd=[];cumQs_sd=[]
    #for folder in folders
    #    data_files=readdir(folder)
    #    data_files=[s for s in data_files if endswith(s,".jld2")]
        for data_file in files#data_files
            @load #=folder*=#data_file  R
            push!(files,#=folder*=#data_file)
            push!(cumIs,median( [sum(R.u[sim][1][:,:])+sum(R.u[sim][2][:,:])+sum(R.u[sim][3][:,:])   for sim=1:size(R.u,1)]))
            push!(cumQs,median( [sum(R.u[sim][5][:,:])   for sim=1:size(R.u,1)]))
            push!(cumIs_sd,std( [sum(R.u[sim][1][:,:])+sum(R.u[sim][2][:,:])+sum(R.u[sim][3][:,:])   for sim=1:size(R.u,1)]))
            push!(cumQs_sd,std( [sum(R.u[sim][5][:,:])   for sim=1:size(R.u,1)]))
        end
    #end
    #df=DataFrame(files=files,cumIs=cumIs,cumIs_sd=cumIs_sd,cumQs=cumQs,cumQs_sd=cumQs_sd)
    #CSV.write(folders[1]*"final_cumulatives__2020_"*Dates.format(now(), "mm_dd__HH_MM") *".csv",df)
    return cumIs,cumQs
end

files=["./Screening/output_withoutH/0_NoInterv/session2_1000sims/noInterv_sc1_1000sims.jld2",
        "./Screening/output_withoutH/1_CTH/session2_1000sims/intervCTH_sc1_1000sims.jld2",
        "./Screening/output_withoutH/1_CTH/session4_1000sims/intervCTH_sc1_1000sims.jld2"#,
        ]
    cumIs,cumQs=final_cumulatives(files)
