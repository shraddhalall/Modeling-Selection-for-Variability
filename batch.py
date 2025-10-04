"""
This example script shows how to set up a batch of selection runs - In this case, 250 runs of mass selection with low weight of epistasis

"""

#Import packages, functions and variables:

from var import *
from setup import *
from selection import *
from sites_fixed import *

import numpy as np

ind = 'mel1' #Batch ID
runs = 250 #Number of iterations of selection

variability = np.zeros((runs,gen-1)) #To store population variability at every generation for all runs 
gt_sites = np.zeros((runs,S,N)) #To store final genotype after selection for all runs 


for r in range(runs): # Iterate over runs
	variability[r], gtsites[r] = mass_sel(True, epi_str = 'low') #Implement mass selection with low level of epistasis

per_fixed,fixation,MAF = sites_fixed(fixed_sites,runs) #parse genotypic data to get site fixation metrics

#Save variables of interest as numpy objects

np.save('var_' + ind + '.npy',variability) #Store variability
np.save('per_fixed_' + ind + '.npy', per_fixed) #Store percentage of sites fixed every run
np.save('fixation_' + ind + '.npy', fixation) #Store average fixation probability of all sites
np.save('MAF_' + ind + '.npy', MAF) #Store frequency of variability allele after selection for all sites