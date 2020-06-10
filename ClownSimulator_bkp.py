#!/usr/bin/env python 
from numpy import *
import numpy as np
np.set_printoptions(suppress=True)
np.set_printoptions(precision=3)  
# import utility functions (e.g. HPD, hist)
from utilities import *
from utility_functions import *

verbose = 0 
colony_size = 5
n_colonies = 10
n_loci = 1000
n_larvae = 10
recomb_freq = 0.01 # avg fraction of alleles recombining
n_time_steps = 500

fraction_of_larvae_spB = 0.001
scenario = 0
replace_no_going_up=0
only_one_event = 0
stop_when_mtDNA_is_gone=0
	
# SPECIAL SCENARIOS
if scenario ==1:
	fraction_of_larvae_spB = 0.3
	only_one_event = 1 # set to one to make only one hybridization event happen
elif scenario ==2:
	fraction_of_larvae_spB = 0.1
	stop_when_mtDNA_is_gone = 1
elif scenario ==3: 
	replace_no_going_up=0.025
	fraction_of_larvae_spB = 0



"""
Increase n_larvae?

fraction_of_larvae_spB should be a fraction of the amount of larvae in the pool (not a probability)

Move recombination to within female/male

amount of admixture could affect death of adults (worst being 50/50 A and B) <- Beta distribution

prob of a larva from the pool to enter anemonis could be a function of admixture as well

make hybridation rate decrease through time?

plot n. pure A individuals, n. admixt and n. pure B nDNA
plot n. pure A individuals and n. pure B mtDNA

OR plot each individual as a proportion of A vs B
color circles based on mtDNA (A: blue, B: red)

	
"""










# probabilities of death 
p_death = 1/(np.linspace(1,10,colony_size))
p_death = p_death/sum(p_death)



# indexes identifying individuals across all colonies
individual_counter = np.arange(n_colonies*colony_size)

species_colony_ind = np.zeros(n_colonies) # species assignemnt per colony

colony_ind = [[i]*colony_size for i in range(n_colonies)]
colony_ind = np.array(colony_ind).flatten()

species_ind = [[species_colony_ind[i]]*colony_size for i in range(n_colonies)]
species_ind = np.array(species_ind).flatten()

# possible_alleles_per_species = species_ind*2
# ideally 0,1 Sp. 0; 2,3 Sp. 1, etc
nDNA1 = np.random.choice(range(2),size=(len(species_ind),n_loci)) # ONLY WORKS for 1 SPECIES
nDNA2 = np.random.choice(range(2),size=(len(species_ind),n_loci))
nDNA = np.array([nDNA1,nDNA2])
nDNA_bck = nDNA + 0
# nDNA[0][1] : All alleles 1 on 2nd individual 
# nDNA[1][1] : All alleles 2 on 2nd individual 
# nDNA[:,0,:] : all paris of alleles for 1st individual
# nDNA[:,colony_ind==0,:] : all data for 1st colony



mtDNA =  np.random.choice(range(2),size=len(species_ind))

# hierarchy
hierarchy = np.array(list(range(colony_size))*n_colonies)

# Ages
ages = np.array(list(range(colony_size))*n_colonies)[::-1]


spB_nucleus = []
spB_mitoc = []
ind_spA_nDNA = []
ind_spA_mtDNA = []



for time_step in range(n_time_steps):
	# run reproduction
	females = individual_counter[hierarchy==0]
	males   = individual_counter[hierarchy==1]
	# choose the random alleles passed on to each larva by f/m (only for nDNA)
	allele_ind_female = np.random.choice([0,1],(n_larvae,n_colonies))
	allele_ind_male   = np.random.choice([0,1],(n_larvae,n_colonies))
	# allele_ind[0,:] : which allele in each larvae in 1st colony
	larvae_pool_nDNA = np.zeros((n_colonies,2,n_larvae,n_loci))
	larvae_pool_mtDNA = np.zeros((n_colonies,n_larvae))
	for j in range(n_colonies):
		female_allele_larvae = nDNA[allele_ind_female[:,j],females[j],:]
		male_allele_larvae   = nDNA[allele_ind_male[:,j],males[j],:]
		larvae_pool_nDNA[j] = np.array([female_allele_larvae,male_allele_larvae])
		larvae_pool_mtDNA[j,:] = mtDNA[females[j]] # 
	
	# print(shape(larvae_pool_nDNA))
	# (3, 2, 10, 100)
	# larvae_pool_nDNA[0] : larvae of 1st colony
	if verbose:		
		print(np.unique(np.sum(larvae_pool_nDNA[0,0],1))) # print unique female alleles
		print(np.unique(np.sum(larvae_pool_nDNA[0,1],1))) # unique male alleles	
	
	# recombination
	recomb_larvae_pool_nDNA  = larvae_pool_nDNA + 0
	female_inherited_alleles = recomb_larvae_pool_nDNA[:,0].flatten() # across all colonies
	female_inherited_alleles_bck = 	female_inherited_alleles + 0
	male_inherited_alleles   = recomb_larvae_pool_nDNA[:,1].flatten()
	
	n_all_larvae = n_loci*n_larvae*n_colonies
	n_recombining_alleles = np.sum(np.random.binomial(1,recomb_freq,n_all_larvae))
	allele_recomb_indx = np.random.choice(np.arange(n_all_larvae), n_recombining_alleles,replace=False)

	female_inherited_alleles[allele_recomb_indx] = male_inherited_alleles[allele_recomb_indx]
	male_inherited_alleles[allele_recomb_indx] = female_inherited_alleles_bck[allele_recomb_indx]
				
	# replace female alleles with male
	recomb_larvae_pool_nDNA[:,0] = female_inherited_alleles.reshape(n_colonies,n_larvae,n_loci)
	# replace male allele with female	
	recomb_larvae_pool_nDNA[:,1] = male_inherited_alleles.reshape(n_colonies,n_larvae,n_loci)	
	if verbose:		
		print("fraction changed by recombination",n_recombining_alleles/n_all_larvae)
		print("fraction changed by recombination (male)",sum(abs(larvae_pool_nDNA[:,1] - recomb_larvae_pool_nDNA[:,1])/np.size(recomb_larvae_pool_nDNA[:,1])))
	# check that the number of ones and zeros didnt change (they were only swapped)
	#print(np.sum(larvae_pool_nDNA)/np.size(larvae_pool_nDNA))
	#print(np.sum(recomb_larvae_pool_nDNA)/np.size(larvae_pool_nDNA))
	
	
	
	if np.random.random() < 1:
		# kill adult individuals
		rnd_colony = np.random.choice(range(n_colonies))
		rnd_individual = np.random.choice(range(colony_size), p=p_death)
		
		affected_colony_nDNA = nDNA[:,colony_ind==rnd_colony,:]+0
		affected_colony_mtDNA = mtDNA[colony_ind==rnd_colony]+0
				
		if np.random.random() < replace_no_going_up:
			nDNA_spB = np.random.choice(np.arange(2,4),(2,n_loci))
			mtDNA_spB = np.random.choice(np.arange(2,4))
			affected_colony_nDNA[:,rnd_individual,:] = nDNA_spB
			affected_colony_mtDNA[rnd_individual] = mtDNA_spB
			#replace_no_going_up = 0
		else:
			ind_going_up = np.arange(rnd_individual+1, colony_size)
			affected_colony_nDNA[:,ind_going_up-1,:] = affected_colony_nDNA[:,ind_going_up,:]	
			affected_colony_mtDNA[ind_going_up-1] = affected_colony_mtDNA[ind_going_up]
		
			# pick a larva
			ind_larva = [np.random.choice(range(n_colonies)), np.random.choice(range(n_larvae))]
			
			if np.random.random() < fraction_of_larvae_spB:
				nDNA_spB = np.random.choice(np.arange(2,4),(2,n_loci)) # values 2, 3
				mtDNA_spB = np.random.choice(np.arange(2,4)) # mt values 2, 3
				affected_colony_nDNA[:,colony_size-1,:] = nDNA_spB
				affected_colony_mtDNA[colony_size-1] = mtDNA_spB
				#print(affected_colony_mtDNA)
				if only_one_event: fraction_of_larvae_spB = 0
			else:
				affected_colony_nDNA[:,colony_size-1,:] = recomb_larvae_pool_nDNA[ind_larva[0],:,ind_larva[1],:] + 0
				affected_colony_mtDNA[colony_size-1] = larvae_pool_mtDNA[ind_larva[0],ind_larva[1]]+0
	
		# reset global DNA variables
		#print(affected_colony_mtDNA)
		nDNA[:,colony_ind==rnd_colony,:] = affected_colony_nDNA
		mtDNA[colony_ind==rnd_colony] = affected_colony_mtDNA

	spB_nucleus.append( size(nDNA[nDNA<=1])/float(np.size(nDNA)))
	spB_mitoc.append(size(mtDNA[mtDNA<=1])/float(np.size(mtDNA)))
	
	# individuals of pure spA nDNA
	a= np.max(nDNA, 2)
	
	np.max(nDNA, 2)
	
	
	
	ind_spA_nDNA.append()
	ind_spA_mtDNA.append()
	
	
	
	if stop_when_mtDNA_is_gone:
		if size(mtDNA[mtDNA<=1])/float(np.size(mtDNA))==0:
			fraction_of_larvae_spB = 0

	if time_step % 100 ==0:  
		if verbose:
			# print sum of full genome (both alleles) per individual : changes only with death/replacement
			print("\nvalue per indvidual (nDNA):", np.sum(np.sum(nDNA,2),0))
			#  print sum of femal genome per individual : changes only with death/replacement
			print("value per indvidual (nDNA,female):", np.sum(nDNA,2)[0])
			# print sum of female alleles across larvae: this changes also with just recombination 
			# without death/replacement
			print(np.sum(recomb_larvae_pool_nDNA[:,0],2).flatten()[0:15])
		
		str_summary = "%s nDNA species A: %s, mtDNA: %s" % (time_step,size(nDNA[nDNA<=1])/float(np.size(nDNA)),size(mtDNA[mtDNA<=1])/float(np.size(mtDNA)))
		print(str_summary)

fig = plt.figure()
fig.add_subplot(221)
plt.plot(range(n_time_steps),spB_nucleus)
fig.add_subplot(222)
plt.plot(range(n_time_steps),spB_mitoc)

plt.show()

"""
 plot n. pure A individuals, n. admixt and n. pure B nDNA
 plot n. pure A individuals and n. pure B mtDNA
 
"""












# import numpy as np
# import pandas as pd
# import matplotlib.pyplot as plt
#
#
# x = np.random.random(100)
# y = np.random.random(100)
#
#
# fig = plt.figure()
#
# fig.add_subplot(221)
# plt.plot(x,y)
# plt.plot(x+0.5,y+0.5,color='red')
# plt.plot(x+0.1,y+0.1,color='green')
# plt.plot(x+0.3,y+0.3,color='orange')
#
# # 2nd plot
# fig.add_subplot(222)
# plt.plot(x,y)
#
# # 2nd plot
# fig.add_subplot(223)
# plt.plot(x,y)
#
# # 2nd plot
# fig.add_subplot(224)
# plt.plot(x,y) (edited)