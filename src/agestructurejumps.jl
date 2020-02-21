"""
Basic representation of the state of the age structured model:
u[wa_index,disease_state,age_group]

This is row major unpacked so
first 1,...,n_wa entries are 0-4 year old susceptibles in wider areas 1,...,n_wa
then n_wa+1,...,2n_wa are 0-4 year old exposeds in urban areas 1,...,n_wa
.
.
.
then n_wa*n_s + 1, ...,  n_wa*n_s + n_wa entries are 5-9 year old susceptibles in wider areas 1,...,n_wa
.
.
.

States:
1 -> S
2 -> E
3 -> I_subclinical
4 -> I_diseased
5 -> H(ospitalised)
6 -> Recovered
7 -> Cumulative I_sub
8 -> Cumulative I_dis
9 -> Cumulative Dead

Events for each wider area :

1-> Urban transmission
2-> Rural transmission
3-> #urban E->A
4-> rural E->A
5-> urban E->D
6-> rural E->D
7-> urban D->H
8-> rural D->H
9-> urban H->R
10-> rural H->R
11-> urban D->R
12-> rural D->R
13-> urban A->R
14-> rural A->R
15-> urban H->death
16-> rural H->death

"""
