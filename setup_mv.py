"""
This file defines functions to set up the starting population with heritable mean and variability in the trait of interest and initialize different selection strategies for selection on variability

It builds upon setup.py from the default model with heritable variability

"""

from var_mv import *
from site_funcs import *
import pickle
import scipy.stats as stats
import numpy as np
rng = np.random.default_rng()


def load_pickle(file_name):
    """ Loads data from .pickle files """
    with open(file_name, 'rb') as fp:
        return pickle.load(fp)

def gen_gamete(parentChr):

    """ 
    Generates a new gamete by simulating recombination events on a parent's two chromosomes;
    Identical to gen_gamete fucntion in the default heritable variability model (script setup.py)
    
    """
    
    recom_bybase = rng.binomial(1,rr,size=(S-1))

    alleles = [0,1]
    a_str = np.random.choice(alleles)
    alleles.remove(a_str)
    for al in alleles:
        b_str = al

    count = 0
    bp = []

    for m,b in enumerate(recom_bybase):
        if b == 1:
            count +=1
            bp.append(m+1)      
    bp.append(S)      

    gam = list(parentChr[a_str][0:bp[0]])

    for c in range(count):
        if (c % 2) == 0:
            items = list(parentChr[b_str][bp[c]:bp[c+1]])
            for x in items:
                gam.append(x)
        else:
            items = list(parentChr[a_str][bp[c]:bp[c+1]])
            for x in items:
                gam.append(x)
                    
    return gam

def gen_0_mass_mv(site_names,var_dict,mean_dict):

    """ 
    Generates the starting population for mass selection with heritable mean and variability 
    Builds upon gen_0_mass() from setup.py in the default model
    """

    animalSites = np.full((N, S, 2), 'null')
    pheno_mat = np.zeros((N, 3))
    vsum = np.zeros(N)
    msum = np.zeros(N)

    for i in range(N):

        temp = np.random.random((S, 1)) < 0.5
        for j, k in enumerate(temp):
            if k:
                animalSites[i, j, 0] = site_names[j]

        temp2 = np.random.random((S, 1)) < 0.5
        for l, m in enumerate(temp2):
            if m:
                animalSites[i, l, 1] = site_names[l]

        indiVar, indiMean = sum_effects(animalSites[i], var_dict, mean_dict)

        indiPheno = np.random.normal(indiMean, indiVar)
        pheno_mat[i] = [indiVar, indiMean, indiPheno]

    animalHistory = np.full((gen, N, S, 2), 'null')
    animalHistory[0] = animalSites

    return animalHistory, pheno_mat

def gen_0_fam_mv(site_names,var_dict,mean_dict):

    """ 
    Generates the starting population for family based selection with heritable mean and variability 
    Builds upon gen_0_fam() from setup.py in the default model
    """

    FemSites_0 = np.full((fam, S, 2), 'null')
    MalSites_0 = np.full((fam, S, 2), 'null')

    animalHistory = np.full((gen, fam, ani, S, 2), 'null')
    pheno_mat = np.zeros((fam, ani, 3))

    FemHist_0 = np.full((fam, fem, S, 2), 'null')
    MalHist_0 = np.full((fam, mal, S, 2), 'null')

    for i in range(fam):
        temp = np.random.random((S, 1)) < 0.5
        for j, k in enumerate(temp):
            if k:
                FemSites_0[i, j, 0] = site_names[j]

        temp2 = np.random.random((S, 1)) < 0.5
        for l, m in enumerate(temp2):
            if m:
                FemSites_0[i, l, 1] = site_names[l]

        temp3 = np.random.random((S, 1)) < 0.5
        for j, k in enumerate(temp3):
            if k:
                MalSites_0[i, j, 0] = site_names[j]

        temp4 = np.random.random((S, 1)) < 0.5
        for l, m in enumerate(temp4):
            if m:
                MalSites_0[i, l, 1] = site_names[l]

        parent1Sites = FemSites_0[i]
        parent2Sites = MalSites_0[i]

        for j in range(fem):
            parent1Chr = parent1Sites.T
            parent2Chr = parent2Sites.T

            g1 = gen_gamete(parent1Chr)
            g2 = gen_gamete(parent2Chr)

            progenySites = np.vstack((g1, g2)).T

            FemHist_0[i, j] = progenySites

        for k in range(mal):
            parent1Chr = parent1Sites.T
            parent2Chr = parent2Sites.T

            g1 = gen_gamete(parent1Chr)
            g2 = gen_gamete(parent2Chr)

            progenySites = np.vstack((g1, g2)).T

            MalHist_0[i, k] = progenySites

    animalHistory[0] = np.hstack((FemHist_0, MalHist_0))
    animalSites = animalHistory[0]

    for i in range(fam):
        for a in range(ani):
            indiVar, indiMean = sum_effects(animalSites[i, a], var_dict, mean_dict)

            indiPheno = np.random.normal(indiMean, indiVar)

            pheno_mat[i, a] = [indiVar, indiMean, indiPheno]

    return animalHistory, pheno_mat
