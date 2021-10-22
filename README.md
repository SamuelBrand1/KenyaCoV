# KenyaCoV - Forecasting for SARS-CoV-2 transmission in Kenya.

This version of the KenyaCoV model is now depreciated. The goal of this modelling exercise was to forecast the potential size of the SARS-CoV-2 epidemic in Kenya before wide-spread transmission and under a wide range of scenarios. Since, the pandemic has arrived in Kenya we have gone forwards with modelling which has been calibrated to Kenyan data.

* J. Ojal et al., Revealing the extent of the first wave of the COVID-19 pandemic in Kenya based on serological and PCR-test data. Wellcome Open Res. 6, 127 (2021) (https://wellcomeopenresearch.org/articles/6-127). With Github repository here: https://github.com/ojal/KenyaSerology.
* S. P. C. Brand et al., COVID-19 transmission dynamics underlying epidemic waves in Kenya. Science, eabk0414 (2021) (https://www.science.org/doi/epdf/10.1126/science.abk0414). With Github repository here: https://github.com/SamuelBrand1/kenya-covid-three-waves .

**Additionally, we have the following projects ongoing:**

* KenyaCoVMultistrains: Mechanistic modelling of the spread of SARS-CoV-2 variants into and around Kenya.
* KenyaCoVaccines: More detailed and realistic modelling of disease burden with vaccination rollout forecasting. 

**The legacy code in this repository has been updated from the original pre-print version to include extra realism in the possible disease pathways of infected individuals. The modelling results of the pre-print can be recovered by the new version of KenyaCoV by setting the average period of pre-symptomatic transmission to 0, and combining mild and mild-then-severe cases. Therefore, we are not maintaining a legacy version of the KenyaCoV model used in the pre-print.**

This repository contains source code for transmission modelling, forecasting, and visualisation for the SARS-CoV-2/COVID-19 epidemic in Kenya.

The main model is based on the wider population groupings of *Wesolowski et al. (2012)*, with explicit age structure. Within population group age mixing is given by the Prem et al estimate for Kenya. A pre-print giving modelling forecasts for Kenya, and describing the underlying mathematical structure of KenyaCoV can now be found on Medrxiv (https://www.medrxiv.org/content/10.1101/2020.04.09.20059865v1).





