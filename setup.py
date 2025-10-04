"""
This file defines functions to set up the starting population and initialize different selection strategies

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
sites = load_pickle('sites.pickle') # Loads the data from the file sites.pickle into the variable sites
sums_map = load_pickle('sums_map.pickle') # Loads the data from the file sums_map.pickle into the variable sums_map

SD_tb = load_pickle('SD_tb.pickle') #Load empirical locomotor handedness variability data
map_1 = np.sort(SD_tb) #sort SD_tb in ascending order for mapping

#Define setup functions:

def SD_map(x, sums, mapping):  # Function to map the values of x to the values of SD_x
    SD_x = np.interp(x, sums, mapping)  # takes the input value x, finds its position within the sums array, and interpolates the corresponding SD value from the mapping array.
    return SD_x


def gen_gamete(parentChr):  
    
    """ Generates a new gamete by simulating recombination events on a parent's two chromosomes """
    
    recom_bybase = rng.binomial(1, rr, size=(S - 1))  # generates an array of binary values (0s and 1s) representing whether a recombination event occurs at each site (except the last one)
   
    alleles = [0, 1]  # Two possible chromosome choices for initializing the gamete
    
    a_str = np.random.choice(alleles)  # randomly chooses one of the two alleles 
    alleles.remove(a_str)  # removes the allele that was randomly chosen from the list of possible alleles
    
    for al in alleles:
        b_str = al  # the remaining allele is assigned to b_str

    count = 0  # counter for the number of recombination events that occur between the parent chromosomes
    bp = []  # list of the breakpoints between the parent chromosomes (position at which the recombination event occurs)

    for m, b in enumerate(recom_bybase):  # for each recombination event
        if b == 1:  # if a recombination event occurs
            count += 1  # increment the count of recombination events
            bp.append(m + 1)  # append the position of the recombination event to the list of breakpoints
    
    bp.append(S)  # append the number of sites to the list of breakpoints

    gam = list(parentChr[a_str][0:bp[0]])  #For chosen parent chromosome a_str, splice the base pairs from start to first break point. This gives the first chunk of the gamete

    for c in range(count):  # for each recombination event
        if (c % 2) == 0:  # if the recombination event is even
            items = list(parentChr[b_str][bp[c]:bp[c + 1]])
            # Extract a segment of the second chromosome starting from the current recombinant breakpoint, and ending at the next recombination breakpoint
            for x in items:  # iterate through the extracted segment
                gam.append(x)  # Append the segment to the gamete
        else:  # if the recombination event is odd
            items = list(parentChr[a_str][bp[c]:bp[c + 1]])
            # Extract a segment of the first chromosome, starting from the current recombinant breakpoint, and ending at the next recombination breakpoint
            for x in items:  # iterate through the extracted segment
                gam.append(x)  # Append the segment to the gamete

    return gam  # returns the new gamete


def gen_0_mass(epi, epi_str):

    """ Generates the starting population for mass selection """
    
    animalSites = np.zeros((N, S, 2))  # Matrix to store site values for each animal

    for i in range(N):  # For each animal
        temp = np.random.random((S, 1)) < 0.5  # Pick S numbers randomly between 0 and 1, and check if they are less than 0.5

        for j, k in enumerate(temp):  # For each number picked
            if k == True:  # If number picked for a site is randomly less than 0.5
                animalSites[i, j, 0] = sites[j]  # Assign that animal to have a non-zero effect size value for that site, for that chromosome

        #Repeat the above process for the second chromosome
        temp2 = np.random.random((S, 1)) < 0.5 # Pick S numbers randomly between 0 and 1, and check if they are less than 0.5
        
        for l, m in enumerate(temp2):  # For each number picked
            if m == True:  # If number picked for a site is randomly less than 0.5
                animalSites[i, l, 1] = sites[l]  # Assign that animal to have a non-zero effect size value for that site


    animalHistory = np.zeros((gen, N, S, 2))  # Matrix to store site values for each animal at every generation

    animalHistory[0] = animalSites  # Initialise the first generation (Gen 0) of the animalHistory matrix with simulated genotypes

    sum_sites_0 = np.sum(np.sum(animalHistory[0], axis=1), axis=1)  # Calculating values of sum of sites for all Gen 0 individuals
   
    if epi == True:  # If epistatis is included in the model 

        if epi_str == 'high':  # High level of epistasis
            sums_overall = load_pickle('sums_overall_h.pickle') # Load the sums_overall_h array into the sums_overall variable
            epi_mat = load_pickle('epi_mat_h.pickle') # Load the epi_mat_h array into the epi_mat variable

        if epi_str == 'mid':  # Medium level of epistasis
            sums_overall = load_pickle('sums_overall_m.pickle') # Load the sums_overall_m array into the sums_overall variable
            epi_mat = load_pickle('epi_mat_m.pickle') # Load the epi_mat_m array into the epi_mat variable

        if epi_str == 'low':  # Low level of epistasis
            sums_overall = load_pickle('sums_overall_l.pickle') # Load the sums_overall_l array into the sums_overall variable
            epi_mat = load_pickle('epi_mat_l.pickle') # Load the epi_mat_l array into the epi_mat variable

        Sites_rs = animalSites.transpose(0, 2, 1)  # Reshape the animalSites array.

        epi_all = []  # To store epistatic effect size values for all individuals
        
        for a in Sites_rs:  # For each animal
            epi_a = 0  # to store the epistatic value for the current individual
            
            for b in a:  # iterate over each chromosome for animal'a'
                indices = np.where(b != 0, 1, 0) # returns an array of 1s and 0s, where 1 indicates that the site has a non-zero effect size

                indices = indices.reshape((S, 1))  # reshapes the indices array into a two-dimensional array

                #Extract the epistatic value for animal 'a' for all interactions on chromosome 'b'
                epi_b = np.sum(np.multiply(np.multiply(epi_mat, indices), indices.T)) # if site has a zero additive effect size, its epistatic value is also 0
                
                epi_a += epi_b  # add value for chromosome 'b' to cumulative epistatic value for animal a

            epi_all.append(epi_a)  # store the total epistatic effects for each animals in the epi_all list
        
        epi_all = np.array(epi_all)  # convert the epi_all list to an array

        s_overall_0 = []  # Create an empty list to store overall genotypic value for all individuals
        for (item1, item2) in zip(sum_sites_0, epi_all): # For each pair of values in the sum_sites_0 and epi_all arrays
            s_overall_0.append(item1 + item2)  # Add the values from the pair and append the result to the s_overall list

        pheno_V_0 = []  # Create an empty list to store mapped genotypic values for all individuals
        for i in s_overall_0:  # For each value in the s_overall array
            v_i = SD_map(i, sums_overall, map_1) # map the value of i to the corresponding empirical phenotypic standard deviation value v_i
            pheno_V_0.append(v_i)

    else:  # If epistatis is not included in the model

        pheno_V_0 = []  # Create an empty list to store mapped genotypic values for all individuals
        for i in sum_sites_0:  # For each value in the sum_sites_0 array
            v_i = SD_map(i, sums_map, map_1) # map the value of i to the corresponding empirical phenotypic standard deviation value v_i

            pheno_V_0.append(v_i)

    pheno_0_B = []  # Create an empty list to store turn bias phenotype for all individuals 
    
    for i in pheno_V_0:  # For each mapped genotypic value i in the pheno_V_0 array
        pheno_i = np.random.normal(mu, i) # Generate a random number from a normal distribution with mean mu and standard deviation i
        
        # Constrain the generated turn bias value to be between 0 and 1 (inclusive)
        if pheno_i < 0:  # If the generated number is less than 0
            pheno_i = 0  # Set the number to 0
        if pheno_i > 1:  # If the generated number is greater than 1
            pheno_i = 1  # Set the number to 1
        
        # Estimate the measured behavior for this individual over several trials via binomial sampling
        pheno_j = (np.random.binomial(no_trials, pheno_i)) / no_trials
        
        pheno_0_B.append(pheno_j) # Append the generated value to the phenotype list

    return animalHistory, pheno_0_B  # Return the genotypic history with animalHistory and the measured phenotypes with pheno_0_B


def gen_0_fam(epi, epi_str):
    """ Generates the starting population for family selection. Creates 'fam' number of families from 'fam' monogamous starting pairs """


    FemSites_0 = np.zeros((fam, S, 2))  # Matrix to store site values for each female parent
    MalSites_0 = np.zeros((fam, S, 2))  # Matrix to store site values for each male parent

    for i in range(fam):  # For each family
        temp = np.random.random((S, 1)) < 0.5  # Pick S number of numbers randomly between 0 and 1, and check if they are less than 0.5
       
        for j, k in enumerate(temp):  # For each number picked
            if k == True:  # If number picked for a site is randomly less than 0.5
                FemSites_0[i, j, 0] = sites[j]  # Assign that female animal to have a non-zero effect size value for that site on the first chromosome

        #Repeat the above process for the second chromosome
        temp2 = np.random.random((S, 1)) < 0.5 # Pick S number of numbers randomly between 0 and 1, and check if they are less than 0.5
      
        for l, m in enumerate(temp2):  # For each number picked
            if m == True:  # If number picked for a site is randomly less than 0.5
                FemSites_0[i, l, 1] = sites[l]  # Assign that female animal to have a non-zero effect size value for that site on the second chromosome

    # Repeat the process above for the male parent:
    for i in range(fam):  # For each family
        temp = np.random.random((S, 1)) < 0.5  # Pick S number of numbers randomly between 0 and 1, and check if they are less than 0.5
        for j, k in enumerate(temp):  # For each number picked
            if k == True:  # If number picked for a site is randomly less than 0.5
                MalSites_0[i, j, 0] = sites[j]  # Assign that male animal to have a non-zero effect size value for that site on the first chromosome
                

        #Repeat the above process for the second chromosome
        temp2 = np.random.random((S, 1)) < 0.5 # Pick S number of numbers randomly between 0 and 1, and check if they are less than 0.5
    
        for l, m in enumerate(temp2):  # For each number picked
            if m == True:  # If number picked for a site is randomly less than 0.5
                MalSites_0[i, l, 1] = sites[l]  # Assign that male animal to have a non-zero effect size value for that site
                

    animalHistory = np.zeros((gen, fam, ani, S, 2)) # Matrix to store genotypes for each animal in each family at every generation

    # Set up families for gen 0
    FemHist_0 = np.zeros((fam, fem, S, 2)) # store the effect size values of each site for each female animal in every family in generation 0.

    MalHist_0 = np.zeros((fam, mal, S, 2)) # store the effect size values of each site for each male animal in every family in generation 0.


    for i in range(fam):  # For each family in generation 0

        parent1Sites = FemSites_0[i]  # Takes effect sizes for each site in each female animal in the current family
        # parent1Sites is a 2D array with shape (S, 2)

        parent2Sites = MalSites_0[i]  # Takes effect sizes for each site in each male animal in the current family
        # parent2Sites is a 2D array with shape (S, 2)

        for j in range(fem):  # To generate each female in the current family (generation 0)
            parent1Chr = parent1Sites.T  # Transposes parent1Sites to make it a 2D array with shape (2, S)

            parent2Chr = parent2Sites.T  # Transposes parent2Sites to make it a 2D array with shape (2, S)

            g1 = gen_gamete(parent1Chr)  # Generates gamete 1 from parent 1
            g2 = gen_gamete(parent2Chr)  # Generates gamete 2 from parent 2

            progenySites = np.vstack((g1, g2))  # Stacks gametes 1 and 2 vertically to create progeny's genotype

            progenySites = progenySites.T  # Transposes progenySites so that the sites become the rows and the alleles from each parent become the columns
 
            FemHist_0[i, j] = progenySites  # Store the genetic information for each female offspring in each family in generation 0 into FemHist_0

        
        for k in range(mal):  # To generate each male in the current family (generation 0)
            parent1Chr = parent1Sites.T  # Transposes parent1Sites to make it a 2D array with shape (2, S)

            parent2Chr = parent2Sites.T  # Transposes parent2Sites to make it a 2D array with shape (2, S)

            g1 = gen_gamete(parent1Chr)  # Generates gamete 1 from parent 1
            g2 = gen_gamete(parent2Chr)  # Generates gamete 2 from parent 2

            progenySites = np.vstack((g1, g2))  # Stacks gametes 1 and 2 vertically to create progeny's genotype

            progenySites = progenySites.T  # Transposes progenySites so that the sites become the rows and the alleles from each parent become the columns

            MalHist_0[i, k] = progenySites # Store the genetic information for each male offspring in each family in generation 0 into MalHist_0

    animalHistory[0] = np.hstack((FemHist_0, MalHist_0))  # Stacks male and female genotypes to store information for all animals in each family

    sum_sites_0 = np.sum(np.sum(animalHistory[0], axis=2), axis=2)  # Sums the effect size values for each site in each animal in each family in generation 0

    if epi == True:  # If epistatis is included in the model

        if epi_str == 'high':  # High level of epistasis
            sums_overall = load_pickle('sums_overall_h.pickle') # Load the sums_overall_h array into the sums_overall variable
            epi_mat = load_pickle('epi_mat_h.pickle') # Load the epi_mat_h array into the epi_mat variable

        if epi_str == 'mid':  # Medium level of epistasis
            sums_overall = load_pickle('sums_overall_m.pickle') # Load the sums_overall_m array into the sums_overall variable
            epi_mat = load_pickle('epi_mat_m.pickle') # Load the epi_mat_m array into the epi_mat variable

        if epi_str == 'low':  # Low level of epistasis
            sums_overall = load_pickle('sums_overall_l.pickle') # Load the sums_overall_l array into the sums_overall variable
            epi_mat = load_pickle('epi_mat_l.pickle') # Load the epi_mat_l array into the epi_mat variable


        animalSites = animalHistory[0]  # Assigns the genetic information for each animal in each family in generation 0

        Sites_rs = animalSites.transpose(0, 1, 3, 2)  # Reshapes the genotype array for gen 0 

        epi_all = np.zeros((fam, ani))  # To store epistatic effect size values for all individuals in each family
        
        for n, f in enumerate(Sites_rs):  # For each family f in Sites_rs
            
            for m, a in enumerate(f):  # For each animal a in family f in Sites_rs (generation 0)
                epi_a = 0  # To store epistatic effect size for animal a
                
                for b in a:  # For each chromosome in animal a
                    indices = np.where(b != 0, 1, 0) # returns an array of 1s and 0s, where 1 indicates that the site has a non-zero effect size
                    
                    indices = indices.reshape((S, 1))  # Reshapes the array

                    #Extract the epistatic value for animal 'a' for all interactions on chromosome 'b'
                    epi_b = np.sum(np.multiply(np.multiply(epi_mat, indices), indices.T)) # if site has a zero additive effect size, its epistatic value is also 0
                   
                    epi_a += epi_b  # add value for chromosome 'b' to cumulative epistatic value for animal a

                    epi_all[n, m] = epi_a   # store the total epistatic effects for each animal in the epi_all matrix

        
        s_overall_0 = np.add(sum_sites_0, epi_all)  #adds the corresponding elements of sum_sites_0 and epi_all to calculate the overall genotypic value for each animal in each family

        pheno_V_0 = np.zeros((fam, ani))   # Create an empty array to store mapped genotypic values for all individuals
        
        for n, i in enumerate(s_overall_0):  # For each family i in s_overall_0
            for m, j in enumerate(i):  # For each animal in family i with genotypic value j
                v_i = SD_map(j, sums_overall, map_1)  # map the value of j to the corresponding empirical phenotypic standard deviation value v_i
                pheno_V_0[n, m] = v_i  # Store the mapped value (v_i) in the pheno_V_0 array at the corresponding family (n) and animal (m) indices in generation 0

    else: # If epistatis is included in the model

        pheno_V_0 = np.zeros((fam, ani))  # Create an empty array to store mapped genotypic values for all individuals
        
        for n, i in enumerate(sum_sites_0):  # For each family in sum_sites_0
            for m, j in enumerate(i):  # For each animal in family i with genotypic value j
                v_i = SD_map(j, sums_map, map_1)  # map the value of j to the corresponding empirical phenotypic standard deviation value v_i
                pheno_V_0[n, m] = v_i  # Stores the mapped value (v_i) in the pheno_V_0 array at the  corresponding family (n) and animal (m) indices in generation 0

   
    pheno_beh_0 = np.zeros((fam, ani))  # Create an empty array to store turn bias phenotype for all individuals 

    for n, i in enumerate(pheno_V_0):  # For each family in pheno_V_0
        for m, j in enumerate(i):  # For each animal family i with mapped value j
            
            pheno_j = np.random.normal(mu, j)  # Generate a random number from a normal distribution with mean mu and standard deviation j

             # Constrain the generated turn bias value to be between 0 and 1 (inclusive)
            if pheno_j < 0:  # If the random number is less than 0, set it to 0
                pheno_j = 0 
            if pheno_j > 1:  # If the random number is greater than 1, set it to 1
                pheno_j = 1  

            # Estimate the measured behavior for this individual over several trials via binomial sampling
            pheno_mes = (np.random.binomial(no_trials, pheno_j)) / no_trials
            
            pheno_beh_0[n, m] = pheno_mes   # Append the generated value to the phenotype matrix

    return animalHistory, pheno_beh_0  # Return the genotypic history with animalHistory and the measured phenotypes with pheno_beh_0