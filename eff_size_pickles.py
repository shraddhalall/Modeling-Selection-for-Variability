"""
This file generates different effect size distributions for sites that determine the mean and variability of a trait, and stores them as pickles to use for downstream selection runs

"""

import numpy as np
import scipy.stats as stats
import pickle

from var_mv import *

def dump_pickle(data,file_name):
    """ dumps data into .pickle files """
    with open(file_name, 'wb') as fp:
        pickle.dump(data,fp)
        
def bp_site_eff(num,fname):
    """
    This function generates site effect sizes with a betaprime distribution with parameters determined by emperical effect size data for varaibility in turn bias
    
    (This function is identical to sites_eff() from pickles.py, used in the default heritable variability model)
    
    num: the number of sites to generate
    fname: in format 'x.pickle', name of file to store the generated effect sizes
    """
    
    eff_size_params = np.array((9.8e-01,5.8,2.8e-07,2.7e-01))
    
    sites_eff = stats.betaprime.rvs(*eff_size_params, size = num)
    for n,i in enumerate(sites_eff):
        while i>0.5:
            i = stats.betaprime.rvs(*eff_size_params, size = 1)
            sites_eff[n] = i

    dump_pickle(sites_eff,fname) #stores the generated effect size data in .pickle format

def exp_site_eff(num,loc,scale,fname):
    """
    This function generates site effect sizes with an exponential distribution, with user-set parameters
    
    num: the number of sites to generate
    loc: of exponential dist
    scale: of exponential dist
    fname: in format 'x.pickle', name of file to store the generated effect sizes
    """
    
    sites_eff = stats.expon.rvs(loc,scale, size = num)

    dump_pickle(sites_eff,fname)





