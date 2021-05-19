#include("/home/gemvi/rabia_aziza_midasnetwork_us/julia_projects/2020-07-22_KenyaCoV/Screening/simFiles_4/distribute.jl")
#include("/home/gemvi/rabia_aziza_midasnetwork_us/julia_projects/2020-07-22_KenyaCoV/Screening/simFiles_4/fixed_K_1000tests/heading_1.jl")
include("heading_1.jl")
#run_intervention_session(1,scenarios,"I3_SympSCT",1,CallbackSet(callback_selection_and_mvt_probabilities,callback_selection_and_mvt_probabilities_symptomatics,callback_SympS,callback_CT),2)

for m=1#:3
	run_intervention_session(m,scenarios,"I7_SympSCT_CTC",m,
		CallbackSet(callback_selection_and_mvt_probabilities,callback_selection_and_mvt_probabilities_symptomatics,
					callback_SympS,callback_CT,callback_screen_population,callback_CT2),n_traj;test_and_trace_contacts=true)
end
