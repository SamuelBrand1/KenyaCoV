push!(LOAD_PATH, joinpath(homedir(),"GitHub/KenyaCoV/src"))
using JLD2,MAT
@load joinpath(homedir(),"Github/KenyaCoVOutputs/results_1.jld2") results_1
@load joinpath(homedir(),"Github/KenyaCoVOutputs/results_3.jld2") results_3

incidence_epsilon_0_no_control = results_1[1][1][1:20,:,1]
incidence_epsilon_0_25_no_control = results_3[1][1][1:20,:,1]
incidence_epsilon_0_25_CI = results_3[6][1][1:20,:,1]

matwrite("data/selected_incidence.mat",Dict("inc_0_nc"=>incidence_epsilon_0_no_control,
                                            "inc_025_nc"=>incidence_epsilon_0_25_no_control,
                                            "inc_025_CI"=>incidence_epsilon_0_25_CI))
