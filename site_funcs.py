"""
This file generates site IDs for mean and variability sites, and generates a chromosome for an individual with mean and variability determining sites

"""

from var_mv import *
import pickle
import random
import numpy as np

def gen_site_names(id, num):
    """
    Generates unique IDs for every site combining a string (representing site type) and a number representing the unique site of that type

    id: string that uniquely identifies the site type
    num: the number of sites to generate
    """
    num_vec = np.arange(0, num).astype(str)

    site_names = []
    for i in num_vec:
        sname = str(id + i)
        site_names.append(sname)

    return site_names

def site_dict(site_names, site_eff):
    """
    Generates a dictionary to store the effect size associated with each unique site ID

    site_names: vector of string names for sites, used as keys
    site_eff: vector of site effect sizes, used as values
    """
    site_dict = dict(zip(site_names, site_eff))
    return site_dict

def random_merge_ordered(L1, L2):
    """
    Given two lists or arrays, generates a third list with all elements from the first two randomly arranged
    """
    L3 = []
    i, j = 0, 0
    while i < len(L1) and j < len(L2):
        if random.choice([True, False]):
            L3.append(L1[i])
            i += 1
        else:
            L3.append(L2[j])
            j += 1
    # Add remaining elements from either list
    L3.extend(L1[i:])
    L3.extend(L2[j:])
    return L3

def chr_layout(v_effs, m_effs):
    """
    Given effect sizes of variability and mean defining sites, generates a chromosome with both site types randomly arranged and dictionaries with respective effect sizes
    """
    v_names = gen_site_names('V', int(S / 2))  # generate IDs for variability sites
    var_dict = site_dict(v_names, v_effs)     # generate dictionary storing variability site IDs and corresponding effect sizes

    m_names = gen_site_names('M', int(S / 2))  # same as above for mean sites
    mean_dict = site_dict(m_names, m_effs)

    site_names = random_merge_ordered(v_names, m_names)  # creates a chromosome layout with variability and mean sites
    return site_names, var_dict, mean_dict

def sum_effects(indSites, var_dict, mean_dict):
    """
    Calculates the additive mean and variability effects for an individual given their chromosomes

    indSites: size S,2. Two chromosomes, S sites each, where each value is either null, or a site name in the keys of one of the dicts
    var_dict: dict containing effect size values for var sites
    mean_dict: dict containing effect size values for mean sites
    """
    indSites = indSites.flatten()
    var_sum = 0
    mean_sum = 0
    for i in indSites:
        if i in var_dict.keys():
            var_sum += var_dict[i]
        elif i in mean_dict.keys():
            mean_sum += mean_dict[i]
        else:
            pass
    return var_sum, mean_sum
