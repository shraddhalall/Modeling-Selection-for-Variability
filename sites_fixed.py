"""
This file defines a function to extract genotypic parameters such as site fixation after selection, from the final generation genotype matrix

"""

#Import packages and variables:
from var import *
from setup import *
import pickle
import scipy.stats as stats
import numpy as np


def sites_fixed(gt_arr, runs):
    """
    Extracts genotypic parameters after selection from the final generation genotype matrix across all runs, i.e. 'gt_arr'
    
    A site is considered to be fixed if, at the end of the selection, both the alleles at the site are the same across all animals
    
    """

    fixed = np.zeros((runs, S)) #To store run-wise value of fixation for each site
    MAF_mat = np.zeros((runs, S)) #To store run-wise value of minor allele frequency for each site

    #Output variables:
    per_fixed = np.zeros((runs))  # Percent of total sites fixed in every run
    fixation = np.zeros((S))  # for each site, the percentage of the runs in which it gets fixed
    MAF = np.zeros((S))  # for each site, the average frequency of the variability allele across all animals in all runs

    for n, r in enumerate(gt_arr): #Iterate over runs
        for m, s in enumerate(r): #Iterate over sites
            s = np.array(s) #Array of size N, with genotype value of site (1, 0.5, 0) in each animal 
            if np.all(s == s[0]): #If genotype of site is identical across all animals
                if (s[0]!=0.5): #Making sure that the site is not heterozygous for all animals
                    fixed[n, m] = 1 #Site is considered fixed for that run
                    
            freq = (np.sum(s)) / N #Total frequency of variability allele at that site = Sum of genotypes for all animals/ Total number of animals
            MAF_mat[n, m] = freq #Store frequency for that site for that run
        
        per = (np.sum(fixed[n]) / S) * 100 #Percentage of sites fixed = Number of sites fixed for that run/ Total number of sites
        per_fixed[n] = per #Store percentage of sites fixed in that run

    fixed_T = fixed.T #make matrix shape sites x runs 

    for n, s in enumerate(fixed_T): #Iterate over sites
        site_fix = (np.sum(s) / runs) * 100 #Percentage of runs in which site s is fixed = Number of runs in which site 's' is fixed/ Total number of runs
        fixation[n] = site_fix # Store fixation probability for site s

    MAF_mat_T = MAF_mat.T #make matrix shape sites x runs

    for n, s in enumerate(MAF_mat_T): #Iterate over sites
        maf_s = np.average(s) #Average frequency of site s over all runs
        MAF[n] = maf_s #Store average frequency of site s

    return per_fixed, fixation, MAF #Return the 3 variables of interest