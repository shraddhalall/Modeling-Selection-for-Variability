""" This file uses the functions in pickles.py to generate new data for use across runs of the model """

#Import packages and variables:

from var import *  # importing all variables from var.py
from pickles import * # import all functions from pickles.py

site_def() #Generate and store site effect sizes

sums_map() #Generates data for linear mapping between simulated additive genotypic values and empirical values of variance in turn bias

epi_map("high") #Generates data for linear mapping between simulated combined (additive + epistatic) genotypic values and empirical values of variance in turn bias, with high level of epistasis

epi_map("mid") #Same as preceeding, with medium level of epistasis

epi_map("low") #Same as preceeding, with low level of epistasis

