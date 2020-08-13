heading_folder="/home/gemvi/rabia_aziza_midasnetwork_us/julia_projects/2020-07-22_KenyaCoV/Screening/simFiles_4/rand_K_1000tests/heading_1.jl"
#heading_folder="C:/Users/rabia/Documents/GitHub/KenyaCoV_cleaned/Screening/simFiles_4/rand_K_1000tests/heading_1.jl"
include(heading_folder)


run_intervention_session(1,[1],"I0_NoInterv",1,CallbackSet(),2)

run_intervention_session(1,[1],"I0_NoInterv",1,CallbackSet(),n_traj)
