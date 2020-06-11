push!(LOAD_PATH, "./src")
    push!(LOAD_PATH, "./Screening")
    using JLD2#DifferentialEquations,JLD2, CSV, DataFrames,Dates,Statistics
    #include("../forecast_functions_screening.jl")
    import KenyaCoV_screening

function make_reports(folders)
    for folder in folders
        data_files=readdir(folder)
        data_files=[s for s in data_files if endswith(s,".jld2")]
        for data_file in data_files
            @load folder*data_file  results
            scenariodata = KenyaCoV_screening.generate_report_screening(folder,results,data_file[1: findfirst(".",data_file)[1]-1],KenyaCoV_screening.counties.county);
        end
    end
end


######

folders=["./Screening/output/0_NoInterv/session1_1sims/",
         "./Screening/output/1_CTH/session1_1sims/",
         "./Screening/output/2_SympS/session1_1sims/","./Screening/output/2_SympS/session2_1sims/",
         #"./Screening/3_SympSCT/session3_1000sims/",#"./Screening/3_SympSCT/session4_1000sims/",
         #"./Screening/4_MS/session3_1000sims/","./Screening/4_MS/session4_1000sims/"#,
         #"./Screening/5_MSCT/session1_500sims/","./Screening/5_MSCT/session2_500sims/"
         ]
    make_reports(folders)
