"""
This file defines variables used to implement selection on variability with both heritable mean and variance. It is built upon the file 'var.py' from the default model with only heritable variability 
"""
N = 980 #Total Population Size
S = 1000 #Total number of sites
gen = 100 #Number of generations for selection
rr = 0.002 #Per base recombination rate
st = 20 #Percent strength of selection

#For mass selection
scr = int(N) #all animals are screened in mass selection
kept = int((st/100)*N)

#For family-based selection
fam = 35 #Number of families
ani = 28 #Number of animals per family
fem = int(ani/2)
mal = int(ani-fem)

scr_f = int(fem) #all animals are screened in each family
scr_m = int(mal) #all animals are screened in each family
scr_ani = scr_f + scr_m
kept_fam = int((st/100)*fam)

#For Full-Sib selection
kept_fem = int(fam/kept_fam) #To generate fam families, we keep kept_fem females from each of the kept_fam families
kept_mal = int(fam/kept_fam)