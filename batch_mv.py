"""
This example script shows how to set up a batch of selection runs with hertiable mean and variability - In this case, 500 runs of basic mass selection with 15x larger effect size sites for mean

"""

#Import packages, functions and variables:

from var_mv import *
from site_funcs import *
from setup_mv import *
from selection_mv import *
import numpy as np

ind = 'mb_exp' #Batch ID
runs = 1 #Number of iterations of selection

variability = np.zeros((runs,gen-1))  #To store population variability at every generation for all runs 
means = np.zeros((runs,gen-1))  #To store population means at every generation for all runs 

#To store genetic value of mean, genetic value of variability and trait phenotype for every individual in every generation for all runs:
phenos = np.zeros((runs,gen,N,3)) #for mass selection
#phenos = np.zeros((runs,gen,fam,ani,3)) #for family selection

for r in range(runs):
	variability[r], means[r], phenos[r] = mass_sel_basic_mv('exp')

np.save('var_' + ind + '.npy',variability)
np.save('mean_' + ind + '.npy', means)
np.save('phenos_' + ind + '.npy', phenos)
