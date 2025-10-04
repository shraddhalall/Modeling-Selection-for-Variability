"""
This file defines functions to carry out mass or family selection over multiple generations 

"""

#Import packages and variables:
from var import *
from setup import *
import pickle
import scipy.stats as stats
import numpy as np

rng = np.random.default_rng()  # Initialize a random number generator

#Load data generated for use across runs
sites = load_pickle('sites.pickle') # Loads the data from the file sites.pickle into the variable sites
sums_map = load_pickle('sums_map.pickle') # Loads the data from the file sums_map.pickle into the variable sums_map

#Load epistatic effect size data: 
sums_overall_h = load_pickle('sums_overall_h.pickle')
epi_mat_h = load_pickle('epi_mat_h.pickle')

sums_overall_m = load_pickle('sums_overall_m.pickle')
epi_mat_m = load_pickle('epi_mat_m.pickle')

sums_overall_l = load_pickle('sums_overall_l.pickle')
epi_mat_l = load_pickle('epi_mat_l.pickle')

#Load empirical data:
SD_tb = load_pickle('SD_tb.pickle') #Load empirical locomotor handedness variability data
map_1 = np.sort(SD_tb) #sort SD_tb in ascending order for mapping

fitness_prob = load_pickle('fitness_prob.pickle') #Load empirical male Drosophila fitness data


def mass_sel(epi: bool, epi_str='No'):
    """This function runs mass selection. 'epi' should be Boolean, with epi_str being 'high', 'low' or 'mid'.
    Returns population variability at each generation and matrix with genotype values of sites for each animal after selection"""
    
    if epi == True:  # If epistasis is included in the model

        if epi_str == 'high':  # for high epistasis
            sums_overall = sums_overall_h
            epi_mat = epi_mat_h

        if epi_str == 'mid': # for medium epistasis
            sums_overall = sums_overall_m
            epi_mat = epi_mat_m

        if epi_str == 'low': # for low epistasis
            sums_overall = sums_overall_l
            epi_mat = epi_mat_l

    animalHistory, pheno_0_B = gen_0_mass(epi, epi_str) #Intitialize starting generation for mass selection with gen_0_mass function defined in setup.py

    animalPheno = np.zeros((gen, N))  # Matrix to store measured behavioural phenotype for each animal at every generation

    animalPheno[0] = pheno_0_B # Initialize phenotypic matrix with gen 0 behavior data

    ani_range = np.arange(0, scr, 1) #List of indices for animals screened from population (Here, all animals are screened)

    varia_beh = [] #To store measured population variability in turn bias at every generation

    for i in range(1, gen): #For every generation after the initial (Gen 0)

        genoMeasured = np.zeros((scr, S, 2)) # To store genotypes for every animal screened in this generation
        phenoMeasured = np.zeros(scr) # To store measured turn bias phentoype for every animal screened 

        animals_scr = np.random.choice(range(N), scr, replace=False) #Here, scr = N so all animals are screened

        for j, k in enumerate(animals_scr): #For each screened animal:
            whichAnimal = k #the index of the animal
            
            genoMeasured[j] = (animalHistory[(i - 1), whichAnimal]) #Extract the genotype of this animal from the animalHistory matrix using its index
            phenoMeasured[j] = (animalPheno[(i - 1), whichAnimal]) #Extract the phenotype of this animal from the animalPheno matrix using its index

        stdev_mea = np.std(phenoMeasured) #Calculate the population variability as standard deviation of all measured turn biases 
        varia_beh.append(stdev_mea) #Store the measured variability for this generation

        screened_fit = stats.norm.fit(phenoMeasured) #Fit measured phenotypes to a normal distribution

        PDF_ratio = [] #To store values of the probability of each animal being from a high variability distribution
        
        for p in phenoMeasured: #For each measured phenotype
           
            #Ratio of measured phenotype being from a target high variability distribution (uniform) vs from the normal distribution with fit parameter estimated above from phenoMeasured:
            PDF_rat_i = (stats.uniform.pdf(p)) / (stats.norm.pdf(p, *(screened_fit))) 
            PDF_ratio.append(PDF_rat_i)

        PDF_ratio = PDF_ratio / np.sum(PDF_ratio) #Normalize calculated ratios

        CDF_scr = np.cumsum(PDF_ratio) #Get CDF from PDF

        picked_ani = [] #To store indices of animals selected
        
        for q in range(kept): #For each animal to be kept for the next generation
            r = np.random.uniform(0, 1) #Draw a random sample from a uniform distribution between 0 & 1
            
            picked_r = ani_range[int(np.argwhere(CDF_scr == min(CDF_scr[(CDF_scr - r) >= 0])))] #Pick the animal that has a high probability of being from a unfirom distribution
            
            while picked_r in picked_ani: #Repeat to get unique picked animals 
                r = np.random.uniform(0, 1)
                picked_r = ani_range[int(np.argwhere(CDF_scr == min(CDF_scr[(CDF_scr - r) >= 0])))]
            picked_ani.append(picked_r)
            
        allInd = list(range(kept)) #indices of all selected individuals
        femInd = np.random.choice(allInd, size=int(kept / 2), replace=False) #Assign half of the animals kept to be female
        malInd = np.setdiff1d(allInd, femInd) #Assign all other animals to be male

        fit_prob = np.random.choice(fitness_prob, size=len(malInd)) #Fitness probabilities for all male individuals
        fit_prob = fit_prob / np.sum(fit_prob) #Normalize fitness probabilities

        for k in range(N): #For each animal 'k' to be generated for the next generation
            parent1ind = (np.random.choice(femInd)) #Choose a female parent from females kept
            parent2ind = (np.random.choice(malInd, p=fit_prob)) #Choose a male parent from males kept, weighted by their fitness

            parent1 = picked_ani[parent1ind] #Get original index of first chosen parent
            parent2 = picked_ani[parent2ind] #Get original index of second chosen parent
            
            parent1Chr = genoMeasured[parent1].T #Extract genotype of first parent from all measured genotypes
            parent2Chr = genoMeasured[parent2].T #Extract genotype of second parent from all measured genotypes

            #Generate gametes for both parents:
            g1 = gen_gamete(parent1Chr) 
            g2 = gen_gamete(parent2Chr)

            progenySites = np.vstack((g1, g2)) #Stack gametets to create animal k's gentoype  

            if epi == True: #If epistasis is included in the model

                epi_k = 0 #To store epistatic value for animal k
                
                for b in progenySites: #For each chromosome for animal k
                    indices = np.where(b != 0, 1, 0) # returns an array of 1s and 0s, where 1 indicates that the site has a non-zero effect size
                    
                    indices = indices.reshape((S, 1))  # reshapes the indices array into a two-dimensional array

                    #Extract the epistatic value for animal 'k' for all interactions on chromosome 'b'
                    epi_b = np.sum(np.multiply(np.multiply(epi_mat, indices), indices.T)) # if site has a zero additive effect size, its epistatic value is also 0
                    epi_k += epi_b  # add value for chromosome 'b' to cumulative epistatic value for animal k

                progenySites = progenySites.T

                animalHistory[i, k] = progenySites #Store progeny k's genotype

                pheno_ss_k = np.sum(progenySites) #Calculate additive genotypic effect for animal k

                pheno_s_overall = (pheno_ss_k) + (epi_k) #Calculate overall (additive + epistatic) genotype for animal k

                pheno_V_k = SD_map(pheno_s_overall, sums_overall, map_1) #Map measured genotype to empirical standard deviation value 

            else: #If epistasis is not included in the model
                
                progenySites = progenySites.T

                animalHistory[i, k] = progenySites #Store progeny k's genotype

                pheno_ss_k = np.sum(progenySites)  #Calculate genotypic effect for animal k

                pheno_V_k = SD_map(pheno_ss_k, sums_map, map_1) #Map measured genotype to empirical standard deviation value 

            pheno_k = np.random.normal(mu, pheno_V_k) #Estimate animal k's true turn bias by drawing from a normal distribution with mean mu and standard deviation = mapped genotypic value 

            # Constrain the generated turn bias value to be between 0 and 1 (inclusive)
            if pheno_k < 0:
                pheno_k = 0
            if pheno_k > 1:
                pheno_k = 1

            # Estimate the measured behavior for this individual over several trials via binomial sampling
            pheno_mea_k = (np.random.binomial(no_trials, pheno_k)) / no_trials
            animalPheno[i, k] = pheno_mea_k #Store the animal's turn bias phentoype 
        
        #End of loop for animals 
    
    #End of loop for generations

    #Estimating site genotype values after selection
    gt_sites = [] #To store genotype of sites after selection
    
    hom = np.zeros((N, S)) 
    het = np.zeros((N, S))
    
    for m, i in enumerate(animalHistory[-1]): #For each individual i in the last generation
        for n, j in enumerate(i): #For each site in individual i
            if np.any(j != 0): #If any site has a non-zero value (i.e. any site has a variability allele) on either chromosome
                gt_sites.append([m, n, j]) #Append that individual,site,genotype to the list of genotypes

    for f in gt_sites: #For each site combination with a variability allele 
        if np.all(f[2] != 0): #If both chromosomes have the variability allele
            hom[f[0], f[1]] = 1 #Assign that animal to be homozygous (GT = 1) for that site
        else: #If only one chromosome has the variability allele
            het[f[0], f[1]] = 0.5 #Assign that animal to be heterozygous (GT = 0.5) for that site

    gt_all = hom + het #Combine all information across animals and sites 
    gt_all = gt_all.T #Transpose combined matrix

    #gt_all will be an S x N matrix where for each site s and each animal i, 
        #gt_all[s,i] = 1 if animal i is homozygous for the variability allele at s, 
        #gt_all[s,i] = 0.5 if i is heterozygous and 
        #gt_all[s,i] = 0 if i is homozygous for the non-variability allele

    return varia_beh, gt_all #Return the measured variability phenotype at every generation and the final genotype values after selection

def mass_sel_alt(epi: bool, epi_str = 'No'):
    
    """
    This function runs the alternative method mass selection, by keeping individuals with the most extreme phenotypes every generation. 
    'epi' should be Boolean, with epi_str being 'high', 'low' or 'mid'. 
    Returns variability after selection (Note: does not return site fixation matrix, although that functionality can be added and will be the same as the function above)
    """
    #Initial section same as function above:
    
    if epi == True:
        
        if epi_str == 'high':
            sums_overall = sums_overall_h
            epi_mat = epi_mat_h
            
        if epi_str == 'mid':
            sums_overall = sums_overall_m
            epi_mat = epi_mat_m
            
        if epi_str == 'low':
            sums_overall = sums_overall_l
            epi_mat = epi_mat_l
                    
    animalHistory, pheno_0_B = gen_0_mass(epi,epi_str)
    
    animalPheno = np.zeros((gen,N)) #Matrix to store measured behavioural phenotype for each animal at every generation

    animalPheno[0] = pheno_0_B

    ani_range = np.arange(0,scr,1)
    
    varia_beh = []

    for i in range (1,gen):
        
        genoMeasured = np.zeros((scr,S,2))
        phenoMeasured = np.zeros(scr)
    
        animals_scr = np.random.choice(range(N), scr, replace = False)
        
        for j,k in enumerate(animals_scr):
            whichAnimal = k
            genoMeasured[j] = (animalHistory[(i-1),whichAnimal]) 
            phenoMeasured[j] = (animalPheno[(i-1),whichAnimal])
    
    
        stdev_mea = np.std(phenoMeasured)
        varia_beh.append(stdev_mea)

        #DISTINCT FROM FUNCTION ABOVE - Implementing alternative mass selection method:
        aniMeasuredSort = np.argsort(phenoMeasured) #Sort animals by measured phentoype
        picked_extreme_low = aniMeasuredSort[0:int(kept/2)] #Pick half the total number of animals to be selected for the next generation from the lower extreme of the distribution
        picked_extreme_high = aniMeasuredSort[-int(kept/2):] #Pick the other half from the upoper extreme of the distribution
        picked_ani = np.concatenate((picked_extreme_low,picked_extreme_high)) #Combine IDs of picked animals
        #END OF DISTINCT SELECTION REGIME

        #After selecting individuals - same steps as function above
        allInd = list(range(kept))
        femInd = np.random.choice(allInd, size = int(kept/2), replace = False)
        malInd = np.setdiff1d(allInd, femInd)

        fit_prob = np.random.choice(fitness_prob, size = len(malInd))
        fit_prob = fit_prob/np.sum(fit_prob)
                         
        for k in range (N):
            parent1ind = (np.random.choice(femInd))
            parent2ind = (np.random.choice(malInd, p = fit_prob))
            
            parent1 = picked_ani[parent1ind]
            parent2 = picked_ani[parent2ind]
            parent1Chr = genoMeasured[parent1].T
            parent2Chr = genoMeasured[parent2].T
        
            g1 = gen_gamete(parent1Chr)
            g2 = gen_gamete(parent2Chr)
            
            progenySites = np.vstack((g1,g2))
        
            if epi == True:
                
                epi_k = 0
                for b in progenySites:
                    indices = np.where(b!=0,1,0)
                    indices = indices.reshape((S,1))
                    epi_b = np.sum(np.multiply(np.multiply(epi_mat,indices),indices.T))
                    epi_k += epi_b
            
                progenySites = progenySites.T
            
                animalHistory[i,k] = progenySites
        
                pheno_ss_k = np.sum(progenySites)
            
                pheno_s_overall = (pheno_ss_k) + (epi_k)

                pheno_V_k = SD_map(pheno_s_overall, sums_overall, map_1)
        
            else:
                progenySites = progenySites.T
            
                animalHistory[i,k] = progenySites
        
                pheno_ss_k = np.sum(progenySites)
            
                pheno_V_k = SD_map(pheno_ss_k, sums_map, map_1)
        
            
            pheno_k = np.random.normal(mu_new,pheno_V_k)
            
            if pheno_k<0:
                pheno_k = 0
            if pheno_k >1:
                pheno_k = 1
        
            pheno_mea_k = (np.random.binomial(no_trials,pheno_k))/no_trials
            animalPheno[i,k] = pheno_mea_k
   
    return varia_beh


def fam_sel(sib, epi: bool, epi_str='No'):
   
    """
    This function runs family selection. 
        'sib' should be 'half' for half-sib selection, 'full' for full-sib selection.
        'epi' should be Boolean, with epi_str being 'high', 'low' or 'mid'.
    Returns population variability at each generation and matrix with genotype values of sites for each animal after selection
    """
    
    if epi == True:  # If epistasis is included in the model

        if epi_str == 'high':  # for high epistasis
            sums_overall = sums_overall_h
            epi_mat = epi_mat_h

        if epi_str == 'mid': # for medium epistasis
            sums_overall = sums_overall_m
            epi_mat = epi_mat_m

        if epi_str == 'low': # for low epistasis
            sums_overall = sums_overall_l
            epi_mat = epi_mat_l

    animalHistory, pheno_beh_0 = gen_0_fam(epi, epi_str) #Intitialize starting generation for family selection with gen_0_fam function defined in setup.py

    FemPheno = np.zeros((gen, fam, fem)) # Matrix to store measured behavioural phenotype for each female animal in every family at every generation
    MalPheno = np.zeros((gen, fam, mal)) # Matrix to store measured behavioural phenotype for each male animal in every family at every generation
    animalPheno = np.zeros((gen, fam, ani))  # Matrix to store measured behavioural phenotype for each animal at every generation

    FemPhenoScr = np.zeros(((gen - 1), fam, scr_f))  #Behavioral phenotypes for females screened from population 
    MalPhenoScr = np.zeros(((gen - 1), fam, scr_m)) #Behavioral phenotypes for males screened from population 
    animalPhenoScr = np.zeros(((gen - 1), fam, scr_ani)) #Behavioral phenotypes for all animals screened from population (Here, all animals are screened)

    animalPheno[0] = pheno_beh_0 # Initialize phenotypic matrix with gen 0 behavior data
    FemPheno[0] = pheno_beh_0[:, 0:fem] # Same as above for females
    MalPheno[0] = pheno_beh_0[:, fem:] # Same as above for males

    means_varia = []  # In each generation, store the mean variability across all families
    pheno_kept = np.zeros(((gen - 1), kept_fam, ani)) #Behavioral phenotypes for animals from the families that are kept after selection

    for i in range(1, gen): #For every generation after the initial (Gen 0)
        
        varFam = np.zeros((2, fam)) #To store index and family-wise variability for each family
        varFam[0] = range(fam) #Assign index to each family

        for f in range(fam): #Iterating over families

            phenof_beh = np.zeros(scr_f) # Initialize matrix to store phentypes of females screened from this family in this generation
            phenom_beh = np.zeros(scr_m) #Same as above for males

            Fem_scr = np.random.choice(range(fem), scr_f, replace=False) #Choose females to be screened from the total pool of females in the family (Here, all females are screened)
            Mal_scr = np.random.choice(range(mal), scr_m, replace=False) #Same as above for males

            for j, k in enumerate(Fem_scr): #For each female screened
                whichAnimal = k # Index of the animal
                phenof_beh[j] = FemPheno[(i - 1), f, whichAnimal] #Extract the phenotype of this animal from the FemPheno matrix using its index

            #Repeat process for males
            for l, m in enumerate(Mal_scr):
                whichAnimal = m
                phenom_beh[l] = MalPheno[(i - 1), f, whichAnimal]

            FemPhenoScr[(i - 1), f] = phenof_beh #Store screened female phenotypes
            MalPhenoScr[(i - 1), f] = phenom_beh #Store screened male phenotypes
            animalPhenoScr[(i - 1), f] = np.hstack((phenof_beh, phenom_beh)) #Store screened animal phenotypes

            varFam[1][f] = np.std(animalPhenoScr[(i - 1), f]) #Calculate variability for family f = standard deviation of turn bias phenotypes for all individuals in the family

        means_varia.append(np.mean(varFam[1])) #Calculate and store population variability = mean of all family-wise variabiltiies

        ind_sort = varFam[1, :].argsort() #Get indices to sort families by their variability 

        varFamSort = varFam[:, ind_sort] #Families sorted by variability

        FamKept = varFamSort[:, -kept_fam:] #Keep kept_fam number of families with the highest variability

        Fem_breed = [] #To store genotypic values of females used to create families for the next generation
        Mal_breed = [] #Same as above for males

        for p, q in enumerate(FamKept[0]): #For each family 'q' kept after selection

            pheno_kept[(i - 1), p] = animalPheno[(i - 1), int(q)] #Extract phenotype of individuals from family 'q' from the animalPheno matrix

            if sib == 'half': #If half-sibling selection is to be implemented
                Fem_breed.append(animalHistory[(i - 1), int(q), 0:fem]) #Extract genotypes for all females from selected family 'q' 
                Mal_breed.append(animalHistory[(i - 1), int(q), fem:]) #Extract genotypes for all males from selected family 'q' 

            if sib == 'full': #If full-sibling selection is to be implemented
                Fem_q = np.random.choice(range(fem), size=kept_fem, replace=False) #Randomly pick kept_fem unique females from this family to generate 'fam' families in the next generation
                Mal_q = np.random.choice(range(mal, mal + fem), size=kept_mal, replace=False) #Randomly pick kept_mal unique males from this family to generate 'fam' families in the next generation

                Fem_breed.append(animalHistory[(i - 1), int(q)][Fem_q]) #Extract genotypes for females kept from family 'q'
                Mal_breed.append(animalHistory[(i - 1), int(q)][Mal_q]) #Extract genotypes for males kept from family 'q'

        Fem_breed = np.array(Fem_breed) #Convert list to array
        Mal_breed = np.array(Mal_breed) #Convert list to array

        #End of loop for families

        #Reshape breeding individual's genotype arrays according to the method of selection
        if sib == 'half': 
            Fem_breed = Fem_breed.reshape((kept_fam * fem), S, 2)
            Mal_breed = Mal_breed.reshape((kept_fam * mal), S, 2)

        if sib == 'full':
            Fem_breed = Fem_breed.reshape((kept_fam * kept_fem), S, 2)
            Mal_breed = Mal_breed.reshape((kept_fam * kept_mal), S, 2)

        fit_prob = np.random.choice(fitness_prob, size=len(Mal_breed)) #Fitness probabilities for all male individuals from empirical data
        fit_prob = fit_prob / np.sum(fit_prob) #Normalize fitness probabilities

        femPar = np.random.choice(range(len(Fem_breed)), size=fam, replace=False) #Choose 'fam' number of females randomly from the breeding pool

        if sib == 'full': #If full-sib selection is to be implemented
            malPar = np.random.choice(range(len(Mal_breed)), size=fam, replace=False) #Choose 'fam' number of males randomly from the breeding pool

        for a in range(fam): #for each family 'a' to be generated for the next generation
            
            parent1ind = femPar[a] #The female parent for the family 

            if sib == 'full': 
                parent2ind = malPar[a] #The male parent for the family in full-sib selection

            for b in range(fem): #For each female offspring in family 'a'
                parent1Sites = Fem_breed[parent1ind] #Genotype of female parent

                if sib == 'half': #If half-sib selection is to be implemented
                    parent2ind = np.random.choice(range(len(Mal_breed)), p=fit_prob) #Choose a male parent for this individual from the breeding pool weighted by their fitness

                parent2Sites = Mal_breed[parent2ind] #Genotype of male parent

                parent1Chr = parent1Sites.T
                parent2Chr = parent2Sites.T

                #Generate gametes for both parents:
                g1 = gen_gamete(parent1Chr)
                g2 = gen_gamete(parent2Chr)

                progenySites = np.vstack((g1, g2)) #Stack gametes to create animal b's gentoype  
                progenySites = progenySites.T

                animalHistory[i, a, b] = progenySites  #Store the animal's genotype

                #End of loop for female progeny

            #Repeat steps above for male progeny
            for c in range(mal):

                parent1Sites = Fem_breed[parent1ind]

                if sib == 'half':
                    parent2ind = np.random.choice(range(len(Mal_breed)), p=fit_prob)

                parent2Sites = Mal_breed[parent2ind]

                parent1Chr = parent1Sites.T
                parent2Chr = parent2Sites.T

                g1 = gen_gamete(parent1Chr)
                g2 = gen_gamete(parent2Chr)

                progenySites = np.vstack((g1, g2))
                progenySites = progenySites.T

                animalHistory[i, a, (fem + c)] = progenySites

                #End of loop for male progeny
        
        #End of loop for families to be generated

        if epi == True: #If epistasis is included in the model
            animalSites_i = animalHistory[i] #Extract genotype of animals generated above
            Sites_rs_i = np.reshape(animalSites_i, (fam, ani, 2, S))  

            epi_all_i = np.zeros((fam, ani)) #To store epistatic effect sizes for all animals in each family
            
            for n, f in enumerate(Sites_rs_i): #For each newly generated family
                for m, a in enumerate(f): #For each animal in the family
                    epi_a = 0 #To store epistatic effect size for animal 'a'
                    for b in a: #For each chromosome in animal 'a'
                        indices = np.where(b != 0, 1, 0) # returns an array of 1s and 0s, where 1 indicates that the site has a non-zero effect size
                        indices = indices.reshape((S, 1)) # reshapes the indices array into a two-dimensional array

                        #Extract the epistatic value for animal 'a' for all interactions on chromosome 'b'
                        epi_b = np.sum(np.multiply(np.multiply(epi_mat, indices), indices.T)) # if site has a zero additive effect size, its epistatic value is also 0
                        epi_a += epi_b # add value for chromosome 'b' to cumulative epistatic value for animal a
                    epi_all_i[n, m] = epi_a #Store epistatic values for all animals in the family

            pheno_ss_i = np.sum(np.sum(animalHistory[i], axis=2), axis=2) #Calculate additive genotypic effect for all individuals

            s_overall_i = np.add(pheno_ss_i, epi_all_i) #Calculate overall (additive + epistatic) genotype for all individuals

            pheno_V_i = np.zeros((fam, ani)) #To store mapped genotypes

            #Map measured genotype to empirical standard deviation value and store for each animal in each family 
            for n, k in enumerate(s_overall_i): 
                for m, j in enumerate(k):
                    v_i = SD_map(j, sums_overall, map_1)
                    pheno_V_i[n, m] = v_i

        else: #If epistasis is not included in the model
            pheno_ss_i = np.sum(np.sum(animalHistory[i], axis=2), axis=2) #Calculate additive genotypic effect for all individuals

            pheno_V_i = np.zeros((fam, ani)) #To store mapped genotypes

            #Map measured genotype to empirical standard deviation value and store for each animal in each family
            for d, e in enumerate(pheno_ss_i):
                for g, h in enumerate(e):
                    v_i = SD_map(h, sums_map, map_1)
                    pheno_V_i[d, g] = v_i

        pheno_beh_i = np.zeros((fam, ani)) #To store measured turn bias behaviors for all animals in all newly generated families

        for n, o in enumerate(pheno_V_i): #For each family 'o'
            for r, s in enumerate(o): #For each mapped genotype value for animals within the family
                pheno_j = np.random.normal(mu, s) #Estimate the animal's true turn bias by drawing from a normal distribution with mean mu and standard deviation = mapped genotypic value 
                # Constrain the generated turn bias value to be between 0 and 1 (inclusive)
                if pheno_j < 0:
                    pheno_j = 0
                if pheno_j > 1:
                    pheno_j = 1

                # Estimate the measured behavior for this individual over several trials via binomial sampling
                pheno_mes = (np.random.binomial(no_trials, pheno_j)) / no_trials
                pheno_beh_i[n, r] = pheno_mes  #Store the animal's turn bias phentoype 
        
        #End of loop for generations

        animalPheno[i] = pheno_beh_i #Store phenotypes measured for all animals generated in this generated
        FemPheno[i] = pheno_beh_i[:, 0:fem]  #Store phenotypes measured for all females generated in this generated
        MalPheno[i] = pheno_beh_i[:, fem:] #Store phenotypes measured for all males generated in this generated

    
    #Estimating site genotype values after selection

    reshaped = animalHistory.reshape(gen, N, S, 2) #reshape the genotype matrix for finding end of selection site genotypes

    gt_sites = [] #To store genotype of sites after selection
    hom = np.zeros((N, S))
    het = np.zeros((N, S))

    for m, i in enumerate(reshaped[-1]): #For each individual i in the last generation
        for n, j in enumerate(i): #For each site j in individual i
            if np.any(j != 0): #If the site has any non-zero value (i.e. has a variability allele) on either chromosome
                gt_sites.append([m, n, j])  #Append that individual,site,genotype to the list of genotypes

    for f in gt_sites: #For each site combination with a variability allele 
        if np.all(f[2] != 0): #If both chromosomes have the variability allele
            hom[f[0], f[1]] = 1  #Assign that animal to be homozygous (GT = 1) for that site
        else: #If only one chromosome has the variability allele
            het[f[0], f[1]] = 0.5 #Assign that animal to be heterozygous (GT = 0.5) for that site

    gt_all = hom + het #Combine all information across animals and sites 
    gt_all = gt_all.T #Transpose combined matrix

    #gt_all will be an S x N matrix where for each site s and each animal i, 
        #gt_all[s,i] = 1 if animal i is homozygous for the variability allele at s, 
        #gt_all[s,i] = 0.5 if i is heterozygous and 
        #gt_all[s,i] = 0 if i is homozygous for the non-variability allele

    return means_varia, gt_all #Return the measured variability phenotype at every generation and the final genotype values after selection