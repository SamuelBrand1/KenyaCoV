include("/home/gemvi/rabia_aziza_midasnetwork_us/julia_projects/2020-07-22_KenyaCoV/Screening/simFiles_4/distribute.jl")
include("/home/gemvi/rabia_aziza_midasnetwork_us/julia_projects/2020-07-22_KenyaCoV/Screening/simFiles_4/fixed_K_1000tests_30/heading_1.jl")

run_intervention_session(1,scenarios,"I4_MS",1,CallbackSet(callback_selection_and_mvt_probabilities,callback_MS),2)

for m=1:3
	run_intervention_session(m,scenarios,"I4_MS",m,CallbackSet(callback_selection_and_mvt_probabilities,callback_MS),n_traj)
end
