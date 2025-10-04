"""
This file defines functions to generate data for use across runs

"""

#Import packages and variables:

from var import *  # importing all variables from var.py
import pickle  # Library that serializes and deserializes python objects
import scipy.stats as stats  # Stats module from scipy library
import numpy as np  # Numpy library

def site_def():
    """ Defines site effect sizes corresponding to an effect size distribution """ 
    
    # S is defined in var
    eff_size_params = np.array((9.76428009e-01, 5.78934275, 2.79397064e-07, 2.71518734e-01))
    # An array of four parameters used in the beta prime distribution
    # Shape parameter, a = 9.76428009e-01
    # Shape parameter, b = 5.78934275
    # Scale parameter, m = 2.79397064e-07
    # Location parameter, loc = 2.71518734e-01
    
    sites = stats.betaprime.rvs(*eff_size_params, size=S) # Generates a random sample containing S random effect sizes from the beta prime distribution defined

    #To ensure that each effect size is less than or equal to 0.5: 
    for n, i in enumerate(sites): # Iterate over each element i in the sites array
        while i > 0.5:  # while the effect size is greater than 0.5
            i = stats.betaprime.rvs(*eff_size_params, size=1) # generate a new random effect size using the same beta prime distribution parameters
            sites[n] = i # Assign the new effect size to the current element of the sites array
            
    with open('sites.pickle', 'wb') as fp:  # Creates and opens a file called sites.pickle in write binary mode
        pickle.dump(sites, fp) # Serializes the sites array and writes it to the sites.pickle file

def sums_map():
    """ Generates a set of potential additive genotypic values for individuals to allow for linear mapping between simulated genotypic values and empirical values of variance in turn bias """
    
    # SD_tb is the mapping data from Ayroles 2015; for a different map, replace the file imported
    with open('SD_tb.pickle', 'rb') as fp:  # Opens the SD_tb.pickle file in read binary mode
        SD_tb = pickle.load(fp)  
    with open('sites.pickle', 'rb') as fp:  # Opens the sites.pickle file in read binary mode
        sites = pickle.load(fp)  
    map_1 = np.sort(SD_tb)  # Sorts the SD_tb NumPy array in ascending order
 
    numAni = 10000 # Sets the number of animals to 10000
    aniPool = np.zeros((numAni, S, 2))  
    
    for i in range(numAni):  # Iterates over each animal
        
        temp = np.random.random((S, 1)) < 0.5 # Pick S numbers randomly between 0 and 1, and check if they are less than 0.5
        
        # Assign non-zero effect sizes to certain animal-site combinations based on the random boolean values generated in the previous step: 
        for j, k in enumerate(temp):  # Iterate over each element k in the temp array, containing boolean values
            if k == True:  # If number picked for a site is randomly less than 0.5
                aniPool[i, j, 0] = sites[j]  # Assign that animal to have a non-zero effect size value for that site, for the first chromosome
                
        #Repeat the steps above for the second chromosome
        temp2 = np.random.random((S, 1)) < 0.5

        for l, m in enumerate(temp2):  # Iterate over each element k in the temp2 array, containing boolean values
            if m == True:  # If number picked for a site is randomly less than 0.5
                aniPool[i, l, 1] = sites[l] # Assign that animal to have a non-zero effect size value for that site
                

    sums_sites = np.sum(np.sum(aniPool, axis=1), axis=1)  # Sums the aniPool array along the first and second axis to obtain total genoypic value for each individual

    sums_sites = np.sort(sums_sites)  # Sorts the sums_S1 array in ascending order

    idx = np.round(np.linspace(0, len(sums_sites) - 1, len(SD_tb))).astype(int)
    # Generates an array of evenly spaced integer numbers over the interval from 0 to len(sums_S1) - 1
    # The number of elements in the array is equal to the length of the SD_tb array

    sample_sums = sums_sites[idx]
    # Creates a new array containing the elements of the sums_S1 array at the indices specified by the idx array
    # This new array has the same length as the SD_tb array and contains elements from sums_S1 at evenly spaced positions

    with open('sums_map.pickle', 'wb') as fp:  # Creates and opens a file called sums_map.pickle in write binary mode
        pickle.dump(sample_sums, fp)  # Dumps the sample_sums array into the sums_map.pickle file


def epi_map(level):  # Level parameter is a string representing the level of epistasis, defined as the ratio of variance in epistatic value to the variance in total genptypic value
    """Generates a set of potential additive + epistatic genotypic values for individuals to allow for linear mapping between simulated genotypic values and empirical values of variance in turn bias"""
    
    if level == 'high':  # If the level is high
        sc = 0.0021  # Set the sc variable to 0.0021. 
        
    if level == 'mid':  # If the level is mid
        seed_run = 430  # Set the seed_run variable to 430
        sc = 0.0021 / np.sqrt(3)  # Set the sc variable to 0.0021 / the square root of 3.

    if level == 'low':  # If the level is low
        seed_run = 268  # Set the seed_run variable to 268
        sc = 0.0021 / np.sqrt(7)  # Set the sc variable to 0.0021 / the square root of 7.

    # SD_tb is the mapping data from Ayroles 2015; for a different map, replace the file imported
    with open('SD_tb.pickle', 'rb') as fp:  # Opens the SD_tb.pickle file in read binary mode
        SD_tb = pickle.load(fp)  # Loads the contents of the SD_tb.pickle file into the SD_tb variable
        
    with open('sites.pickle', 'rb') as fp:  # Opens the sites.pickle file in read binary mode
        sites = pickle.load(fp)  # Loads the contents of the sites.pickle file into the sites variable
    
    numAni = 10000  # Sets the number of animals to 10000
    aniPool = np.zeros((numAni, S, 2))  
    
    #Additive Effects (same procedure as sums_map):
    for i in range(numAni):  # Iterate over each animal
        temp = np.random.random((S, 1)) < 0.5 # Pick S number of numbers randomly between 0 and 1, and check if they are less than 0.5
        
        for j, k in enumerate(temp):  # Iterate over each element k in the temp array, containing boolean values
            if k == True:  # If number picked for a site is randomly less than 0.5
                aniPool[i, j, 0] = sites[j]  # Assign that animal to have a non-zero effect size value for that site, for the first chromosome
 
        #Repeat the steps above for the second chromosome
        temp2 = np.random.random((S, 1)) < 0.5

        for l, m in enumerate(temp2):  # Iterate over each element m in the temp2 array, containing boolean values
            if m == True:  # If number picked for a site is randomly less than 0.5
                aniPool[i, l, 1] = sites[l] # Assign that animal to have a non-zero effect size value for that site

    sums_sites = np.sum(np.sum(aniPool, axis=1), axis=1) # Sums the aniPool array along the first and second axis to obtain additive genoypic value for each individual
    # 1D array where each element represents the total effect size for a specific animal
    sums_sites = np.sort(sums_sites)  # Sorts the elements of the sums_sites array in ascending order

    # Epistatic Effects:
    
    epi_val_1 = -stats.expon.rvs(loc=0, scale=sc, size=int(((S * S) - S) / 4))
    # Generates an array of negative random samples from an exponential distribution with the specified location and scale parameters.
    # The array size is set to generate values for filling half of the non-diagonal elements of a symmetric SxS matrix of pairwise interactions

    epi_val_2 = stats.expon.rvs(loc=0, scale=sc, size=int(((S * S) - S) / 4))
    # Generates an array of positive random samples from an exponential distribution with the specified location and scale parameters.

    epi_val = np.concatenate((epi_val_1, epi_val_2)) # Concatenates the two arrays of random samples into a single array
    np.random.shuffle(epi_val)  # Randomly shuffles the elements of the epi_val array

    
    epi_mat = np.zeros((S, S))  # Creates a 2D array of zeros with dimensions S x S to hold pairwise epistatic interaction values
    
    xs, ys = np.triu_indices(S, 1)  # Returns the row & column indices for the upper-triangle of the (S x S) array
    
    epi_mat[xs, ys] = epi_val  # Assign the shuffled epi_val array values to the upper triangle of the epi_mat matrix
    epi_mat[ys, xs] = epi_val  # Assign the shuffled epi_val array values to the lower triangle of the epi_mat matrix
    # Makes epi_mat a symmetric matrix.

    pool_rs = aniPool.transpose(0, 2, 1)  # Rearrange the elements of the aniPool array

    epi_all = [] # To store epistatic effect sizes for all individuals
    for a in pool_rs:  # iterate over each animal
        epi_a = 0 # to store the epistatic value for the current individual
        for b in a:  # iterate over each chromosome for animal'a'
            
            indices = np.where(b != 0, 1, 0) # returns an array of 1s and 0s, where 1 indicates that the site has a non-zero additive effect size
            
            indices = indices.reshape((S, 1))  # reshapes the indices array into a two-dimensional array

            #Extract the epistatic value for animal 'a' for all interactions on chromosome 'b'
            epi_b = np.sum(np.multiply(np.multiply(epi_mat, indices), indices.T)) # if site has a zero additive effect size, its epistatic value is also 0
        
            epi_a += epi_b # add value for chromosome 'b' to cumulative epistatic value for animal a
            
        epi_all.append(epi_a)  # Append the total epistatic effect size for each animal to the epi_all list
    
    epi_all = np.array(epi_all)  # Converts the epi_all list to a numpy array

    #Combined genotypic effect:
    s_overall = sums_sites + epi_all  # adds the sums_sites array to the epi_all array
    #s_overall array contains the total combined effect (additive and epistasis effects) for each animal.
    
    sums_sort = np.sort(s_overall)  # sorts the s_overall array in ascending order
    
    idx_epi = np.round(np.linspace(0, len(sums_sort) - 1, len(SD_tb))).astype(int)
    # Generates an array of evenly spaced integer numbers over the interval from 0 to len(sums_sort) - 1
    # The number of elements in the array is equal to the length of the SD_tb array
    sums_overall = sums_sort[idx_epi]  # creates a new array containing the elements of the sums_sort array
    # at the indices specified by the idx_epi array

    if level == 'high':  # if the level of epistasis is high
        with open('sums_overall_h.pickle', 'wb') as fp:  # creates and opens a pickle file in write mode
            pickle.dump(sums_overall, fp)  # dumps the sums_overall array into the pickle file
        with open('epi_mat_h.pickle', 'wb') as fp:  # creates and opens a pickle file in write mode
            pickle.dump(epi_mat, fp)  # dumps the epi_mat array into the pickle file
    
    if level == 'mid':  # if the level of epistasis is medium
        with open('sums_overall_m.pickle', 'wb') as fp:  # creates and opens a pickle file in write mode
            pickle.dump(sums_overall, fp)  # dumps the sums_overall array into the pickle file
        with open('epi_mat_m.pickle', 'wb') as fp:  # creates and opens a pickle file in write mode
            pickle.dump(epi_mat, fp)  # dumps the epi_mat array into the pickle file
    
    if level == 'low':  # if the level of epistasis is low
        with open('sums_overall_l.pickle', 'wb') as fp:  # creates and opens a pickle file in write mode
            pickle.dump(sums_overall, fp)  # dumps the sums_overall array into the pickle file
        with open('epi_mat_l.pickle', 'wb') as fp:  # creates and opens a pickle file in write mode
            pickle.dump(epi_mat, fp)  # dumps the epi_mat array into the pickle file
