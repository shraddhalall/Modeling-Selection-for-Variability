"""
This file defines functions to set up the starting population and initialize different selection strategies for selection on mean value of a trait

It builds upon setup.py from the default model of selection on variability, but has no epistatis included

"""

#Import packages and variables:

from var import *  # Import all variables from var.py
import pickle  # For saving and loading data
import scipy.stats as stats  # For using statistical functions
import numpy as np  # For using numpy arrays

rng = np.random.default_rng()  # Initialize the random number generator

def load_pickle(file_name):
    """ Loads data from .pickle files """
    with open(file_name, 'rb') as fp:
        return pickle.load(fp)


#Load data generated for use across runs
sites = load_pickle('sites.pickle')  #The mean determining sites have the same effect sizes as the variability determining sites in the default model
sums_map = load_pickle('sums_map.pickle')

SD_tb = load_pickle('SD_tb.pickle') #Load empirical locomotor handedness variability data -  this is used as mean value mapping data for the trait of interest
map_1 = np.sort(SD_tb) #sort SD_tb in ascending order for mapping

sd_fix = 0.04 #Fixed value of standard deviation of trait for drawing measured phenotype with no heritable variability

#Define setup functions:

def SD_map(x,sums,mapping):
    SD_x = np.interp(x, sums, mapping)
    return SD_x

def gen_gamete(parentChr):
    """ Generates a new gamete by simulating recombination events on a parent's two chromosomes; identical to gen_gamete() in setup.py """
    
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

def gen_0_mass_mean(): 
    
    animalSites = np.zeros((N,S,2)) #Matrix to store site values for each animal

    for i in range(N):
        temp =  np.random.random((S,1)) < 0.5 #Pick S number of numbers randomly between 0 and 1, and check if they are less than 0.5 
        for j,k in enumerate(temp):
            if k == True: #If number picked for a site is randomly less than 0.5
                animalSites[i,j,0] = sites[j] #Assign that animal to have a non-zero effect size value for that site
            
        temp2 = np.random.random((S,1)) < 0.5
        for l,m in enumerate(temp2):
            if m == True: #If number picked for a site is randomly less than 0.5
                animalSites[i,l,1] = sites[l] 
            
    animalHistory = np.zeros((gen,N,S,2)) #Matrix to store genotypes for each animal at every generation

    animalHistory[0]=animalSites #Initialise matrix with parent genotype

    sum_sites_0 = np.sum(np.sum(animalHistory[0], axis = 1), axis = 1) #For calculating parent values of sum of sites

    map_1 = np.sort(SD_tb)

    pheno_V_0 = []
    for i in sum_sites_0:
        v_i = SD_map(i, sums_map, map_1)
        pheno_V_0.append(v_i)            
    
    pheno_0_B = []
    for i in pheno_V_0:
        pheno_i = np.random.normal(i,sd_fix) #Draw phenotype with mean determined by genetic sum of sites and fixed standard deviation
        if pheno_i<0:
            pheno_i = 0
        if pheno_i >1:
            pheno_i = 1 
      
        pheno_j = (np.random.binomial(no_trials,pheno_i))/no_trials
        pheno_0_B.append(pheno_j)
    
    return animalHistory, pheno_0_B

def gen_0_fam_mean():
    
    FemSites_0 = np.zeros((fam,S,2))
    MalSites_0 = np.zeros((fam,S,2))
    
    for i in range(fam):
        temp =  np.random.random((S,1)) < 0.5 #Pick S number of numbers randomly between 0 and 1, and check if they are less than 0.5 
        for j,k in enumerate(temp):
            if k == True: #If number picked for a site is randomly less than 0.5
                FemSites_0[i,j,0] = sites[j] #Assign that animal to have a non-zero effect size value for that site
            
        temp2 = np.random.random((S,1)) < 0.5
        for l,m in enumerate(temp2):
            if m == True: #If number picked for a site is randomly less than 0.5
                FemSites_0[i,l,1] = sites[l] 

    for i in range(fam):
        temp =  np.random.random((S,1)) < 0.5 #Pick S number of numbers randomly between 0 and 1, and check if they are less than 0.5 
        for j,k in enumerate(temp):
            if k == True: #If number picked for a site is randomly less than 0.5
                MalSites_0[i,j,0] = sites[j] #Assign that animal to have a non-zero effect size value for that site
            
        temp2 = np.random.random((S,1)) < 0.5
        for l,m in enumerate(temp2):
            if m == True: #If number picked for a site is randomly less than 0.5
                MalSites_0[i,l,1] = sites[l] 
        
    
    animalHistory = np.zeros((gen,fam,ani,S,2)) #Matrix to store genotypes for each animal in each family at every generation

    #Set up families for gen 0
    FemHist_0 = np.zeros((fam,fem,S,2))
    MalHist_0 = np.zeros((fam,mal,S,2)) 
    
    for i in range(fam):
    
        parent1Sites = FemSites_0[i]
        parent2Sites = MalSites_0[i]
    
        for j in range(fem):
        
            parent1Chr = parent1Sites.T
            parent2Chr = parent2Sites.T
        
            g1 = gen_gamete(parent1Chr)
            g2 = gen_gamete(parent2Chr)
            
            progenySites = np.vstack((g1,g2))
            
            progenySites = progenySites.T
        
            FemHist_0[i,j] = progenySites
        
        for k in range(mal):
        
            parent1Chr = parent1Sites.T
            parent2Chr = parent2Sites.T
        
            g1 = gen_gamete(parent1Chr)
            g2 = gen_gamete(parent2Chr)
            
            progenySites = np.vstack((g1,g2))
            
            progenySites = progenySites.T
            
            MalHist_0[i,k] = progenySites
    
    animalHistory[0] = np.hstack((FemHist_0,MalHist_0))
    sum_sites_0 = np.sum(np.sum(animalHistory[0], axis = 2), axis = 2)
    map_1 = np.sort(SD_tb)
    
    pheno_V_0 = np.zeros((fam, ani))
    for n,i in enumerate(sum_sites_0):
        for m,j in enumerate(i):
            v_i = SD_map(j, sums_map, map_1)
            pheno_V_0[n,m] = v_i
    

    pheno_beh_0 = np.zeros((fam,ani))

    for n,i in enumerate(pheno_V_0):
        for m,j in enumerate(i):
            pheno_j = np.random.normal(j,sd_fix) #Draw phenotype with mean determined by genetic sum of sites and fixed standard deviation
            if pheno_j<0:
                pheno_j = 0
            if pheno_j >1:
                pheno_j = 1 
    
            pheno_mes = (np.random.binomial(no_trials,pheno_j))/no_trials
            pheno_beh_0[n,m] = pheno_mes
    
    return animalHistory, pheno_beh_0
