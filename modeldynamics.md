# Enumeration of the possible social/age/infection states in KenyaCoV

The basic object tracking the state of the age structured model is `u[i,a,state]`, which gives the number of individuals in area/social group `i`, age group `a` and infections state, `state`.

### Enumeration of infection states

States:
1. S(usceptibles)
2. E(xposed/latent infected)
3. P(re-symptomatic infected)
4. A(symptomatic)
5. M(ild) symptomatics
6. first mild then eventually (se)V(ere) symptomatics
7. H(ospitalised)
8. Recovered
9. Cumulative P->A
10. Cumulative P->M
11. Cumulative P->V
12. Cumulative V->H

### Enumeration of events for stochastic version of model

1. S to E
2. E to A
3. E to P
4. P to M
5. P to V
6. V to H
7. M to R
8. A to R
