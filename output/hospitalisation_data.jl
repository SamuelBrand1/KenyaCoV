#Estimate from CDC COVID-19 Response Team. Severe Outcomes Among Patients with Coronavirus Disease 2019 (COVID-19) - United States, February 12-March 16, 2020. MMWR Morb. Mortal. Wkly. Rep. 69, 343–346 (2020).
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


bar(hosp_rate_by_age,lab = "",xticks = (1:2:17,))
