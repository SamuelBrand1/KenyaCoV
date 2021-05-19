heading_folder="/home/gemvi/rabia_aziza_midasnetwork_us/julia_projects/2020-07-22_KenyaCoV/ScreeningR/simFiles_4/rand_K_1000tests/heading_1.jl"
#heading_folder="C:/Users/rabia/Documents/GitHub/KenyaCoV_cleaned/Screening/simFiles_4/rand_K_1000tests/heading_1.jl"
include(heading_folder)

run_intervention_session(1,scenarios,"I3_SympSCT",1,CallbackSet(callback_selection_and_mvt_probabilities,callback_selection_and_mvt_probabilities_symptomatics,callback_SympS,callback_CT),2)

for m=2
	run_intervention_session(m,scenarios,"I3_SympSCT",m,CallbackSet(callback_selection_and_mvt_probabilities,callback_selection_and_mvt_probabilities_symptomatics,callback_SympS,callback_CT),n_traj)
end
