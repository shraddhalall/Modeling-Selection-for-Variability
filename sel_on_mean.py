"""
This file defines functions to carry out mass or family selection on the mean value of a trait, with heritable mean and no heritable variability, over multiple generations 

It builds upon selection.py from the default model of selection on variability, but has no epistasis and site fixation estimation included

"""
#Import packages and variables
from var import *
from setup_mean import *
import pickle
import scipy.stats as stats
import numpy as np
rng = np.random.default_rng()

sites = load_pickle('sites.pickle') #The mean determining sites have the same effect sizes as the variability determining sites in the default model
sums_map = load_pickle('sums_map.pickle')

SD_tb = load_pickle('SD_tb.pickle') #Load empirical locomotor handedness variability data - this is used as mean value mapping data for the trait of interest
map_1 = np.sort(SD_tb) #sort SD_tb in ascending order for mapping

fitness_prob = load_pickle('fitness_prob.pickle') #Load empirical male Drosophila fitness data

sd_fix = 0.04 #Fixed value of standard deviation of trait for drawing measured phenotype with no heritable variability
    
def mass_sel_mean():
    
    """This function implements mass selection for increasing the mean value of the trait. Returns population mean values at each generation"""
    
    animalHistory, pheno_V_0 = gen_0_mass_mean()
    
    animalPheno = np.zeros((gen,N)) #Matrix to store measured phenotype for each animal at every generation

    animalPheno[0] = pheno_V_0

    ani_range = np.arange(0,scr,1)
    
    mean_beh = []

    for i in range (1,gen):
        
        genoMeasured = np.zeros((scr,S,2))
        phenoMeasured = np.zeros(scr)
    
        animals_scr = np.random.choice(range(N), scr, replace = False)
        
        for j,k in enumerate(animals_scr):
            whichAnimal = k
            genoMeasured[j] = (animalHistory[(i-1),whichAnimal]) 
            phenoMeasured[j] = (animalPheno[(i-1),whichAnimal])
    
        mean_mea = np.mean(phenoMeasured) #Average value of the trait measured across all animals in the population
        
        mean_beh.append(mean_mea)

        #Selection on increasing mean - pick animals with the highest meausured phentoypes
        picked_ani = []
        pheno_sorted = np.sort(phenoMeasured)[::-1] #Sort animals by the measured phenotype
        pheno_ind_sorted = np.argsort(phenoMeasured)[::-1] 
        picked_ani = pheno_ind_sorted[0:(kept+1)]
        
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
        
            progenySites = progenySites.T
            
            animalHistory[i,k] = progenySites
        
            pheno_ss_k = np.sum(progenySites)
                
            pheno_V_k = SD_map(pheno_ss_k, sums_map, map_1)
                
            pheno_k = np.random.normal(pheno_V_k,sd_fix)  #Draw phenotype with mean determined by genetic sum of sites and fixed standard deviation
            
            if pheno_k<0:
                pheno_k = 0
            if pheno_k >1:
                pheno_k = 1
        
            pheno_mea_k = (np.random.binomial(no_trials,pheno_k))/no_trials
            animalPheno[i,k] = pheno_mea_k

    return mean_beh
    
def fam_sel_mean(sib):
    
    """This function implements family selection for increasing the mean value of the trait. Returns population mean values at each generation"""
    
    animalHistory, pheno_V_0 = gen_0_fam_mean()
    FemPheno = np.zeros((gen,fam,fem))
    MalPheno = np.zeros((gen,fam,mal))
    animalPheno = np.zeros((gen,fam, ani)) #Matrix to store measured behavioural phenotype for each animal at every generation

    FemPhenoScr = np.zeros(((gen-1),fam,scr_f))
    MalPhenoScr = np.zeros(((gen-1),fam,scr_m))
    animalPhenoScr = np.zeros(((gen-1),fam, scr_ani))


    animalPheno[0] = pheno_V_0
    FemPheno[0] = pheno_V_0[:,0:fem]
    MalPheno[0] = pheno_V_0[:,fem:]

    means_trait = [] #In each generation, stores the mean trait across 35 families
    pheno_kept = np.zeros(((gen-1),kept_fam,ani))

    for i in range (1,gen):
        meanFam = np.zeros((2,fam))
        meanFam[0] = range(fam)
    
        for f in range(fam):
                
            phenof_beh = np.zeros(scr_f)
            phenom_beh = np.zeros(scr_m)

            Fem_scr = np.random.choice(range(fem), scr_f, replace = False)
            Mal_scr = np.random.choice(range(mal), scr_m, replace = False)
        
            for j,k in enumerate(Fem_scr):
        
                whichAnimal = k
            
                phenof_beh[j] = FemPheno[(i-1),f,whichAnimal]
            
            for l,m in enumerate(Mal_scr):
        
                whichAnimal = m
          
                phenom_beh[l] = MalPheno[(i-1),f,whichAnimal]
        
            FemPhenoScr[(i-1),f] = phenof_beh
            MalPhenoScr[(i-1),f] = phenom_beh
            animalPhenoScr[(i-1),f] = np.hstack((phenof_beh,phenom_beh))
        
            meanFam[1][f] = np.mean(animalPhenoScr[(i-1),f])
    
        means_trait.append(np.mean(meanFam[1]))

        
        #Implement selection for increasing mean value of trait - keep families with highest average value of trait   
        ind_sort = meanFam[1,:].argsort()
        meanFamSort = meanFam[:,ind_sort]
        FamKept = meanFamSort[:,-kept_fam:]
    
        Fem_breed = []
        Mal_breed = []
    
        for p,q in enumerate(FamKept[0]):
        
            pheno_kept[(i-1),p] = animalPheno[(i-1),int(q)]
        
            if sib == 'half':
                Fem_breed.append(animalHistory[(i-1),int(q),0:fem])
                Mal_breed.append(animalHistory[(i-1),int(q),fem:])
                
            if sib == 'full':
                Fem_q = np.random.choice(range(fem), size = kept_fem, replace = False)
                Mal_q = np.random.choice(range(mal, mal+fem), size = kept_mal, replace = False)
            
                Fem_breed.append(animalHistory[(i-1),int(q)][Fem_q])
                Mal_breed.append(animalHistory[(i-1),int(q)][Mal_q])
                
        Fem_breed = np.array(Fem_breed)
        Mal_breed = np.array(Mal_breed)
        
        if sib == 'half':
            Fem_breed = Fem_breed.reshape((kept_fam*fem),S,2)
            Mal_breed = Mal_breed.reshape((kept_fam*mal),S,2)
    
        if sib == 'full':
            Fem_breed = Fem_breed.reshape((kept_fam*kept_fem),S,2)
            Mal_breed = Mal_breed.reshape((kept_fam*kept_mal),S,2)
            
        fit_prob = np.random.choice(fitness_prob, size = len(Mal_breed))
        fit_prob = fit_prob/np.sum(fit_prob) 
           
        femPar = np.random.choice(range(len(Fem_breed)), size = fam, replace = False)
        
        if sib == 'full':
            malPar = np.random.choice(range(len(Mal_breed)), size = fam, replace = False)
        
        for a in range (fam):
            progeny_beh_list = []
            parent1ind = femPar[a]
            
            if sib == 'full':
                parent2ind = malPar[a]
                       
            for b in range(fem):
                parent1Sites = Fem_breed[parent1ind]
                
                if sib == 'half':
                    parent2ind = np.random.choice(range(len(Mal_breed)), p = fit_prob)
               
                parent2Sites = Mal_breed[parent2ind]
            
                parent1Chr = parent1Sites.T
                parent2Chr = parent2Sites.T
        
                g1 = gen_gamete(parent1Chr)
                g2 = gen_gamete(parent2Chr)
            
                progenySites = np.vstack((g1,g2))

                progenySites = progenySites.T
                progeny_SS = np.sum(progenySites)

                progeny_V = SD_map(progeny_SS,sums_map, map_1)
                    
                pheno_beh_A = np.random.normal(progeny_V,sd_fix)  #Draw phenotype with mean determined by genetic sum of sites and fixed standard deviation
                if pheno_beh_A<0:
                    pheno_beh_A = 0
                if pheno_beh_A >1:
                    pheno_beh_A = 1 
    
                progeny_beh = (np.random.binomial(no_trials,pheno_beh_A))/no_trials
        
                progeny_beh_list.append(progeny_beh)
        
                animalHistory[i,a,b] = progenySites
        
            for c in range(mal):
            
                parent1Sites = Fem_breed[parent1ind]
                
                if sib == 'half':
                    parent2ind = np.random.choice(range(len(Mal_breed)), p = fit_prob)
                
                parent2Sites = Mal_breed[parent2ind]
            
                parent1Chr = parent1Sites.T
                parent2Chr = parent2Sites.T
        
                g1 = gen_gamete(parent1Chr)
                g2 = gen_gamete(parent2Chr)
            
                progenySites = np.vstack((g1,g2))
               
                progenySites = progenySites.T
                progeny_SS = np.sum(progenySites)

                progeny_V = SD_map(progeny_SS,sums_map, map_1)
                    
                pheno_beh_A = np.random.normal(progeny_V,sd_fix)  #Draw phenotype with mean determined by genetic sum of sites and fixed standard deviation
                if pheno_beh_A<0:
                    pheno_beh_A = 0
                if pheno_beh_A >1:
                    pheno_beh_A = 1 
    
                progeny_beh = (np.random.binomial(no_trials,pheno_beh_A))/no_trials
        
                progeny_beh_list.append(progeny_beh)
        
                animalHistory[i,a,(fem+c)] = progenySites

        
        pheno_ss_i = np.sum(np.sum(animalHistory[i], axis = 2), axis = 2)

        pheno_V_i = np.zeros((fam, ani))
        for d,e in enumerate(pheno_ss_i):
            for g,h in enumerate(e):
                v_i = SD_map(h, sums_map, map_1)
                pheno_V_i[d,g] = v_i
        
        pheno_beh_i = np.zeros((fam,ani))

        for n,o in enumerate(pheno_V_i):
            for r,s in enumerate(o):
                pheno_j = np.random.normal(s,sd_fix)  #Draw phenotype with mean determined by genetic sum of sites and fixed standard deviation
                if pheno_j<0:
                    pheno_j = 0
                if pheno_j >1:
                    pheno_j = 1 
    
                pheno_mes = (np.random.binomial(no_trials,pheno_j))/no_trials
                pheno_beh_i[n,r] = pheno_mes

        animalPheno[i] = pheno_beh_i
        FemPheno[i] = pheno_beh_i[:,0:fem]
        MalPheno[i] = pheno_beh_i[:,fem:]
          
    return means_trait
