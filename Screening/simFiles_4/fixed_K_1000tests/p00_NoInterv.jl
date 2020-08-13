include("/home/gemvi/rabia_aziza_midasnetwork_us/julia_projects/2020-07-22_KenyaCoV/Screening/simFiles_4/distribute.jl")
include("/home/gemvi/rabia_aziza_midasnetwork_us/julia_projects/2020-07-22_KenyaCoV/Screening/simFiles_4/fixed_K_1000tests/heading_1.jl")


run_intervention_session(1,[1],"I0_NoInterv",1,CallbackSet(),2)

run_intervention_session(1,[1],"I0_NoInterv",1,CallbackSet(),n_traj)
