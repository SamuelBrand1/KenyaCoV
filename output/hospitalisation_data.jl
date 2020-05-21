#Estimate from CDC COVID-19 Response Team. Severe Outcomes Among Patients with Coronavirus Disease 2019 (COVID-19) - United States, February 12-March 16, 2020. MMWR Morb. Mortal. Wkly. Rep. 69, 343–346 (2020).

using CSV

# hosp_capacity = CSV.read(joinpath(homedir(),"Documents/Covid-19/jl_models/KenyaCoVOutputs/Health_system_capacity_data_Kenya.csv"))
hosp_capacity = CSV.read(joinpath(homedir(),"Github/KenyaCoV/data/Health_system_capacity_data_Kenya.csv"))

spare_capacity_H_by_county = (hosp_capacity[:,2].*hosp_capacity[:,5])[2:end]
spare_capacity_ICU_by_county = (hosp_capacity[:,3].*hosp_capacity[:,6])[2:end]


age_cats = ["0-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39","40-44","45-49","50-54","55-59","60-64","65-69","70-74","75-79","80+"]

verity_IFR = [0.00161,0.00695,0.0309,
            0.0844,
            0.161,
            0.595,
            1.93,
            4.28,
            7.8]./100


hosp_rate_CDC = [mean([1.6,2.5]),#0-19 year olds
                 mean([14.3,20.8]),#20-44 yos
                 mean([21.2,28.3]),#45-54
                 mean([20.5,30.1]),#55-64
                 mean([28.6,43.5]),#65-74
                 mean([30.5,58.7]),#75-84
                 mean([31.3,70.3])]./100 #85+

ICU_rate_CDC = [0.,#0-19 year olds
              mean([2.,4.2]),#20-44 yos
              mean([5.4,10.4]),#45-54
              mean([4.7,11.2]),#55-64
              mean([8.1,18.8]),#65-74
              mean([10.5,31.0]),#75-84
              mean([6.3,29.])]./100 #85+


hosp_rate_by_age = [hosp_rate_CDC[1],hosp_rate_CDC[1],hosp_rate_CDC[1],hosp_rate_CDC[1],
                    hosp_rate_CDC[2],hosp_rate_CDC[2],hosp_rate_CDC[2],hosp_rate_CDC[2],hosp_rate_CDC[2],
                    hosp_rate_CDC[3],hosp_rate_CDC[3],
                    hosp_rate_CDC[4],hosp_rate_CDC[4],
                    hosp_rate_CDC[5],hosp_rate_CDC[5],
                    hosp_rate_CDC[6],hosp_rate_CDC[7]]

ICU_rate_by_age = [ICU_rate_CDC[1],ICU_rate_CDC[1],ICU_rate_CDC[1],ICU_rate_CDC[1],
                    ICU_rate_CDC[2],ICU_rate_CDC[2],ICU_rate_CDC[2],ICU_rate_CDC[2],ICU_rate_CDC[2],
                    ICU_rate_CDC[3],ICU_rate_CDC[3],
                    ICU_rate_CDC[4],ICU_rate_CDC[4],
                    ICU_rate_CDC[5],ICU_rate_CDC[5],
                    ICU_rate_CDC[6],ICU_rate_CDC[7]]


ICU_rate_by_age_cond_hosp = ICU_rate_by_age./hosp_rate_by_age
ICU_rate_by_age_cond_hosp.*0.625


# plt_hosp_estimates = bar(hosp_rate_by_age,lab = "",xticks = (1:2:17,age_cats[1:2:17]),
#     title = "Hospitalisation rate by age (CDC estimate)",xlabel = "Age of infected")
# plt_ICU_estimates = bar(ICU_rate_by_age_cond_hosp,lab = "",xticks = (1:2:17,age_cats[1:2:17]),
#     title = "ICU rate by age, given hospitalisation (CDC estimate)",xlabel = "Age of infected")
# savefig(plt_hosp_estimates,"output/hosp_estimates_by_age.png")
# savefig(plt_ICU_estimates,"output/ICU_estimates_by_age.png")

# IFR = 0.9*d_1.*hosp_rate_by_age.*ICU_rate_by_age_cond_hosp*0.625
# bar(IFR)
# verity_comparison = zeros(17)
# for a = 1:8
#     verity_comparison[2*(a-1)+1] = verity_IFR[a]
#     verity_comparison[2*(a-1)+2] = verity_IFR[a]
# end
# verity_comparison[17] = verity_IFR[9]
# brand_verity_comparison = groupedbar(hcat(IFR,verity_comparison),lab = ["Brand et al." "Verity et al."],
#             legend = :topleft,xticks = (1:2:17,age_cats[1:2:17]),
#             title = "Infection Fatality Ratio predictions",
#             ylabel = "Prob. of mortality per infection",
#             xlabel = "Age group")
#
# brand_verity_rel_comparison = bar(verity_comparison./(IFR.+0.00001),lab = "",
#             legend = :topleft,xticks = (1:2:17,age_cats[1:2:17]),
#             title = "Infection Fatality Ratio predictions",
#             ylabel = "OR (Brand vs Verity)",
#             xlabel = "Age group")
#
# savefig(brand_verity_comparison,"data/IFR_comparison_Brand_vs_verity.png")
