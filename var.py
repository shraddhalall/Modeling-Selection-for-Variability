"""
This file sets up the starting population conditions
"""


# Define input variables

N = 980 #Total population size
S = 500  # Total number of sites
gen = 100  # Number of generations for selection
sex = True  # For deciding if the model is sex-based
rr = 0.002  # Per base recombination rate (i.e., 0.2% chance that a recombination event will occur at any given base pair position in the DNA sequence)
st = 20  # Percent strength of selection
no_trials = 500  # Trials for each individual

# For mass selection
scr = int(N)  # Number of individuals screened from the population (Here, all individuals are screened)
kept = int((st/100)*N)  # Number of individuals kept after selection

# For family-based selection
fam = 35  # Number of families
ani = 28  # Number of animals per family
fem = int(ani/2)  # Number of females per family 
mal = int(ani-fem)  # Number of males per family

scr_f = int(fem)  # Number of females screened per family (Here, all are screened)
scr_m = int(mal)  # Number of males screened per family
scr_ani = scr_f + scr_m  # Total number of individuals screened
kept_fam = int((st/100)*fam)  # Number of families kept

# For Full-Sib selection
kept_fem = int(fam/kept_fam)  #To generate fam families, we keep kept_fem females from each of the kept_fam families
kept_mal = int(fam/kept_fam)

# Depending upon the behavior
mu = 0.5  # The mean of the behavioural phenotype of interest. For locomotor handedness, mean = 0.5
