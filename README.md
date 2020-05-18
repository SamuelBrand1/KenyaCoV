# KenyaCoV - Forecasting for SARS-CoV-2 transmission in Kenya.

**We are currently updating the KenyaCoV model to include extra realism (e.g. a pre-symptomatic transmission period for infected individuals, differentiation of "symptomatic" cases into mild or mild-then-severe), more spatially resolved simulation, a regression method for tuning parameters, and forecasts for hospital and ICU usage. This has introduced breaking changes into the code base which will be resolved ASAP.**

**The modelling results of the pre-print can be recovered by the new version of KenyaCoV by setting the average period of pre-symptomatic transmission to 0, and combining mild and mild-then-severe cases. Therefore, we are not maintaining a legacy version of the KenyaCoV model used in the pre-print.**

This repository contains source code for transmission modelling, forecasting, and visualisation for the SARS-CoV-2/COVID-19 epidemic in Kenya.

The main model is based on the wider population groupings of *Wesolowski et al. (2012)*, with explicit age structure. Within population group age mixing is given by the Prem et al estimate for Kenya. A pre-print giving modelling forecasts for Kenya, and describing the underlying mathematical structure of KenyaCoV can now be found on Medrxiv (https://www.medrxiv.org/content/10.1101/2020.04.09.20059865v1).

*src/run_KenyaCoV_agestructuredmodel.jl* demonstrates basic methods for running KenyaCoV simulations. There will be more detailed notebook examples to follow. 




