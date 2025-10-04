"""
This file defines functions to carry out mass or family selection on variability, for a trait with heritable mean and variability, over multiple generations 

It builds upon selection.py from the default model for a trait with only heritable variability
"""

from var_mv import *
from setup_mv import *
from site_funcs import *
import pickle
import scipy.stats as stats
import numpy as np
rng = np.random.default_rng()

# load all the pickles that will be used across selection regimes
# (will be imported into run script with 'from selection import *')

fitness_prob = load_pickle('fitness_prob.pickle') #Load empirical male Drosophila fitness data

#Load site effect sizes for mean and variability effects drawn from the same betaprime distribution
var_effs_bp = load_pickle('var_effs_bp.pickle')
mean_effs_bp = load_pickle('mean_effs_bp.pickle')

#Load site effect sizes for mean and variability effects drawn from different exponential distribution, with mean effects 15x larger than variability effects
var_effs_exp = load_pickle('var_effs_exp.pickle')
mean_effs_exp = load_pickle('mean_effs_exp.pickle')

def mass_sel_basic_mv(eff):

    """
    Implments the default mass selection method, with uniform sampling acorss the population for selection on variability
    Builds upon mass_sel() from selection.py in the default model
    
    eff: Set to 'bp' for identical mean and variability effects, 'exp' for mean effects 15x larger than variability effects
    """

    if eff == 'exp':
        v_effs = var_effs_exp
        m_effs = mean_effs_exp
    elif eff == 'bp':
        v_effs = var_effs_bp
        m_effs = mean_effs_bp
    else:
        print ('Please select exp or bp for effect size distributions')
		

    # set up chromosome layout 
    site_names, var_dict, mean_dict = chr_layout(v_effs, m_effs)
    
    animalHistory, pheno_0_mat = gen_0_mass_mv(site_names,var_dict,mean_dict) #set up generation 0
    
    animalPheno = np.zeros((gen, N, 3))  # Matrix to store indivar, indimean & measured phenotype for each animal at every generation
    animalPheno[0] = pheno_0_mat

    ani_range = np.arange(0, scr, 1)
    
    varia_beh = []  # variability across animals in behavior every generation
    mean_beh = []   # mean of measured behavior every generation

    for i in range(1, gen):
        beh_values = animalPheno[i - 1][:, 2]
        max_beh = np.max(beh_values)
        min_beh = np.min(beh_values)
        
        genoMeasured = np.full((scr, S, 2), 'null')
        phenoMeasured = np.zeros((scr, 3))
    
        animals_scr = np.random.choice(range(N), scr, replace=False)
        
        for j, k in enumerate(animals_scr):
            whichAnimal = k
            genoMeasured[j] = animalHistory[i - 1, whichAnimal]
            phenoMeasured[j] = animalPheno[i - 1, whichAnimal]
    
        stdev_mea = np.std(phenoMeasured[:, 2])
        mean_mea = np.mean(phenoMeasured[:, 2])
        varia_beh.append(stdev_mea)
        mean_beh.append(mean_mea)
    
        screened_fit = stats.norm.fit(phenoMeasured[:, 2])
    
        PDF_ratio = []
        for p in phenoMeasured[:, 2]:
            PDF_rat_i = (stats.uniform.pdf(p, loc=min_beh, scale=max_beh - min_beh)) / (stats.norm.pdf(p, *screened_fit))
            PDF_ratio.append(PDF_rat_i)
    
        PDF_ratio = PDF_ratio / np.sum(PDF_ratio)
    
        CDF_scr = np.cumsum(PDF_ratio)

        picked_ani = []
        
        for q in range(kept):
            r = np.random.uniform(0, 1)
            pick_ind = np.argwhere(CDF_scr == min(CDF_scr[(CDF_scr - r) >= 0]))
            pick_ind = pick_ind[0][0]
            picked_r = ani_range[int(pick_ind)]
            
            while picked_r in picked_ani:
                r = np.random.uniform(0, 1)
                pick_ind = np.argwhere(CDF_scr == min(CDF_scr[(CDF_scr - r) >= 0]))
                pick_ind = pick_ind[0][0]
                picked_r = ani_range[int(pick_ind)]
                
            picked_ani.append(picked_r)
        
        allInd = list(range(kept))
        femInd = np.random.choice(allInd, size=int(kept / 2), replace=False)
        malInd = np.setdiff1d(allInd, femInd)

        fit_prob = np.random.choice(fitness_prob, size=len(malInd))
        fit_prob = fit_prob / np.sum(fit_prob)
                         
        for k in range(N):
            parent1ind = np.random.choice(femInd)
            parent2ind = np.random.choice(malInd, p=fit_prob)
            
            parent1 = picked_ani[parent1ind]
            parent2 = picked_ani[parent2ind]

            parent1Chr = genoMeasured[parent1].T
            parent2Chr = genoMeasured[parent2].T
        
            g1 = gen_gamete(parent1Chr)
            g2 = gen_gamete(parent2Chr)
            
            progenySites = np.vstack((g1, g2))
            progenySites = progenySites.T
            
            animalHistory[i, k] = progenySites

            indiVar, indiMean = sum_effects(progenySites, var_dict, mean_dict)

            indiPheno = np.random.normal(indiMean, indiVar)
            animalPheno[i, k] = [indiVar, indiMean, indiPheno]
   
    return varia_beh, mean_beh, animalPheno

def mass_sel_alt_mv(eff):
    """
    Implments the alternative mass selection method, keeping individuals with extreme phenotypes for selection on variability
    Builds upon mass_sel_alt() from selection.py in the default model
    
    eff: Set to 'bp' for identical mean and variability effects, 'exp' for mean effects 15x larger than variability effects
    """
    if eff == 'exp':
        v_effs = var_effs_exp
        m_effs = mean_effs_exp
    elif eff == 'bp':
        v_effs = var_effs_bp
        m_effs = mean_effs_bp
    else:
        print ('Please select exp or bp for effect size distributions')

    # set up chromosome layout 
    site_names, var_dict, mean_dict = chr_layout(v_effs, m_effs)

    animalHistory, pheno_0_mat = gen_0_mass_mv(site_names,var_dict,mean_dict) #set up Gen 0
    
    animalPheno = np.zeros((gen, N, 3))  # Matrix to store indivar, indimean & measured phenotype for each animal at every generation
    animalPheno[0] = pheno_0_mat

    ani_range = np.arange(0, scr, 1)
    
    varia_beh = []  # variability across animals in behavior every generation
    mean_beh = []   # mean of measured behavior every generation

    for i in range(1, gen):
        genoMeasured = np.full((scr, S, 2), 'null')
        phenoMeasured = np.zeros((scr, 3))
    
        animals_scr = np.random.choice(range(N), scr, replace=False)
        
        for j, k in enumerate(animals_scr):
            whichAnimal = k
            genoMeasured[j] = animalHistory[i - 1, whichAnimal]
            phenoMeasured[j] = animalPheno[i - 1, whichAnimal]
    
        stdev_mea = np.std(phenoMeasured[:, 2])
        mean_mea = np.mean(phenoMeasured[:, 2])
        varia_beh.append(stdev_mea)
        mean_beh.append(mean_mea)
    
        aniMeasuredSort = np.argsort(phenoMeasured[:, 2])
        picked_extreme_low = aniMeasuredSort[0:int(kept / 2)]
        picked_extreme_high = aniMeasuredSort[-int(kept / 2):]
        picked_ani = np.concatenate((picked_extreme_low, picked_extreme_high))
        
        allInd = list(range(kept))
        femInd = np.random.choice(allInd, size=int(kept / 2), replace=False)
        malInd = np.setdiff1d(allInd, femInd)

        fit_prob = np.random.choice(fitness_prob, size=len(malInd))
        fit_prob = fit_prob / np.sum(fit_prob)
                         
        for k in range(N):
            parent1ind = np.random.choice(femInd)
            parent2ind = np.random.choice(malInd, p=fit_prob)
            
            parent1 = picked_ani[parent1ind]
            parent2 = picked_ani[parent2ind]

            parent1Chr = genoMeasured[parent1].T
            parent2Chr = genoMeasured[parent2].T
        
            g1 = gen_gamete(parent1Chr)
            g2 = gen_gamete(parent2Chr)
            
            progenySites = np.vstack((g1, g2))
            progenySites = progenySites.T
            
            animalHistory[i, k] = progenySites

            indiVar, indiMean = sum_effects(progenySites, var_dict, mean_dict)

            indiPheno = np.random.normal(indiMean, indiVar)
            animalPheno[i, k] = [indiVar, indiMean, indiPheno]
   
    return varia_beh, mean_beh, animalPheno

def fam_sel_mv(sib,eff):
    """
    Implments family selection, keeping individuals from families with the highest variability, for selection on variability
    Builds upon fam_sel() from selection.py in the default model
    
    sib: Set to 'half' for half-sib selection, 'full' for full-sib selection
    eff: Set to 'bp' for identical mean and variability effects, 'exp' for mean effects 15x larger than variability effects
    """
    
    if eff == 'exp':
        v_effs = var_effs_exp
        m_effs = mean_effs_exp
    elif eff == 'bp':
        v_effs = var_effs_bp
        m_effs = mean_effs_bp
    else:
        print ('Please select exp or bp for effect size distributions')

        
    # set up chromosome layout 
    site_names, var_dict, mean_dict = chr_layout(v_effs, m_effs)
        
    animalHistory, pheno_0_mat = gen_0_fam_mv(site_names,var_dict,mean_dict) #set up Gen 0
    
    FemPheno = np.zeros((gen, fam, fem, 3))
    MalPheno = np.zeros((gen, fam, mal, 3))
    animalPheno = np.zeros((gen, fam, ani, 3))  # Matrix to store indivar, indimean & measured phenotype for each animal at every generation
    
    animalPheno[0] = pheno_0_mat
    FemPheno[0] = pheno_0_mat[:, 0:fem]
    MalPheno[0] = pheno_0_mat[:, fem:]

    means_varia = []  # In each generation, what is the average phenotypic variability across 35 families
    means_mean = []   # In each generation, what is the average phenotypic mean across 35 families

    for i in range(1, gen):
        varFam = np.zeros((2, fam))
        varFam[0] = range(fam)

        meanFam = []
    
        for f in range(fam):
            phenof_beh = np.zeros(scr_f)
            phenom_beh = np.zeros(scr_m)

            Fem_scr = np.random.choice(range(fem), scr_f, replace=False)
            Mal_scr = np.random.choice(range(mal), scr_m, replace=False)
        
            for j, k in enumerate(Fem_scr):
                whichAnimal = k
                phenof_beh[j] = FemPheno[i - 1, f, whichAnimal, 2]
            
            for l, m in enumerate(Mal_scr):
                whichAnimal = m
                phenom_beh[l] = MalPheno[i - 1, f, whichAnimal, 2]

            pheno_beh_all = np.concatenate((phenof_beh, phenom_beh))
        
            varFam[1][f] = np.std(pheno_beh_all)
            meanFam.append(np.mean(pheno_beh_all))
    
        means_varia.append(np.mean(varFam[1]))
        means_mean.append(np.mean(meanFam))
          
        ind_sort = varFam[1, :].argsort()
    
        varFamSort = varFam[:, ind_sort]
        
        FamKept = varFamSort[:, -kept_fam:]
    
        Fem_breed = []
        Mal_breed = []
    
        for p, q in enumerate(FamKept[0]):
            if sib == 'half':
                Fem_breed.append(animalHistory[i - 1, int(q), 0:fem])
                Mal_breed.append(animalHistory[i - 1, int(q), fem:])
            if sib == 'full':
                Fem_q = np.random.choice(range(fem), size=kept_fem, replace=False)
                Mal_q = np.random.choice(range(mal, mal + fem), size=kept_mal, replace=False)
                Fem_breed.append(animalHistory[i - 1, int(q)][Fem_q])
                Mal_breed.append(animalHistory[i - 1, int(q)][Mal_q])
                
        Fem_breed = np.array(Fem_breed)
        Mal_breed = np.array(Mal_breed)
        
        if sib == 'half':
            Fem_breed = Fem_breed.reshape((kept_fam * fem), S, 2)
            Mal_breed = Mal_breed.reshape((kept_fam * mal), S, 2)
    
        if sib == 'full':
            Fem_breed = Fem_breed.reshape((kept_fam * kept_fem), S, 2)
            Mal_breed = Mal_breed.reshape((kept_fam * kept_mal), S, 2)
            
        fit_prob = np.random.choice(fitness_prob, size=len(Mal_breed))
        fit_prob = fit_prob / np.sum(fit_prob)
           
        femPar = np.random.choice(range(len(Fem_breed)), size=fam, replace=False)
        
        if sib == 'full':
            malPar = np.random.choice(range(len(Mal_breed)), size=fam, replace=False)
        
        for a in range(fam):
            parent1ind = femPar[a]
            if sib == 'full':
                parent2ind = malPar[a]
                       
            for b in range(fem):
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
    
                sumVar, sumMean = sum_effects(progenySites, var_dict, mean_dict)
    
                if map:
                    indiVar = trait_map(sumVar, sums_range, indiVar_range)
                    indiMean = trait_map(sumMean, sums_range, indiMean_range)
                else:
                    indiVar = sumVar
                    indiMean = sumMean
    
                indiPheno = np.random.normal(indiMean, indiVar)
    
                FemPheno[i, a, b] = [indiVar, indiMean, indiPheno]
                animalHistory[i, a, b] = progenySites
            
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
        
                indiVar, indiMean = sum_effects(progenySites, var_dict, mean_dict)
    
                indiPheno = np.random.normal(indiMean, indiVar)
    
                MalPheno[i, a, c] = [indiVar, indiMean, indiPheno]
                animalHistory[i, a, (fem + c)] = progenySites
    
        animalPheno[i] = np.concatenate((FemPheno[i], MalPheno[i]), axis=1)
  
    return means_varia, means_mean, animalPheno
