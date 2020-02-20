"""
Basic representation of the state of the age 
For the age structured model the 3-dim representation should
be unpacked as a vector

This is row major unpacked so
first 1,...,n entries are susceptibles in urban areas 1,...,n
then n+1,...,2n are exposeds in urban areas 1,...,n
.
.
.
then n*num_states + 1, ..., n*num_states + n entries are susceptibles in rural areas 1,...,n
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

Events for each location/county:

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
