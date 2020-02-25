# KenyaCoV

Code for simulating Wuhan CoV in Kenya.

There are two main models:
1. A county level model which distinguishes between a rural population (assumed to be immobile) and an urban population (assumed to be mobile around the country). Population flux is parameterised using a fit to the wider level data given in Wesolowski et al.
2. A model based on the wider population groupings of Wesolowski et al. In this model 20-50 year olds are assumed to be mobile, therefore the rural/urban population structure is discarded. However, explicit age structure in the population is introduced. Within population group age mixing is given by the Prem et al estimate for Kenya.

Disease state enumeration:
* 1 -> S
* 2 -> E
* 3 -> I_subclinical
* 4 -> I_diseased
* 5 -> H(ospitalised)
* 6 -> Recovered
* 7 -> Cumulative I_sub
* 8 -> Cumulative I_dis
* 9 -> Cumulative Dead

*src/run_KenyaCoV_basicmodel.jl* demonstrates basic methods for running a KenyaCoV simulation version 1 (spatial + urban/rural). *src/run_KenyaCoV_basicmodel.jl* demonstrates basic methods for running KenyaCoV simulation version 2 (spatial + age structured).

In both simulation versions the state of the epidemic is represented as a multi-dimensional array:
* Simulation 1: $X(n,s,u)$ is the number of individuals (discrete) in the $n$th county, in disease state $s$, in urban/rural environment $u$.
* Simulation 2:
