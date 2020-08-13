

include("/home/gemvi/rabia_aziza_midasnetwork_us/julia_projects/2020-07-22_KenyaCoV/Screening/simFiles_4/distribute.jl")
include("/home/gemvi/rabia_aziza_midasnetwork_us/julia_projects/2020-07-22_KenyaCoV/Screening/simFiles_4/fixed_K_500tests/heading_1.jl")

run_intervention_session(1,[1],"I1_CTH",1,CallbackSet(callback_selection_and_mvt_probabilities,callback_detect_H,callback_CT),2)

for m=1:3
	run_intervention_session(m,scenarios,"I1_CTH",m,CallbackSet(callback_selection_and_mvt_probabilities,callback_detect_H,callback_CT),n_traj)
end
