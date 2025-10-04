""" This file uses the functions in eff_size_pickles.py to generate new data for mean and variability site effect sizes for use across runs of the model """

# Import packages and variables:

from var_mv import *  # importing all variables from var_mv.py
from eff_size_pickles import * # import all functions from eff_size_pickles.py

# Identical underlying distributions for mean and variability site effects
bp_site_eff(int(S/2),'var_effs_bp.pickle') #generates and stores site effect sizes distributed as a betaprime function, for the half of the total sites (i.e., sites that determine variability)

bp_site_eff(int(S/2),'mean_effs_bp.pickle') #as above for mean determining sites


# Mean site effects are 15x larger than variability site effects

exp_site_eff(int(S/2),0,0.07,'var_effs_exp.pickle') #generates and store site effect sizes determined by an exponential distribution with given parameters for variability sites

exp_site_eff(int(S/2),1,0.15,'mean_effs_exp.pickle') #as above for mean sites, scaled to be ~15x larger effects than variability sites