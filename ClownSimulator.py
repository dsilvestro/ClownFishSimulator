#!/usr/bin/env python 
from numpy import *
import numpy as np
np.set_printoptions(suppress=True)
np.set_printoptions(precision=3)  
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.axes

verbose = 0 
colony_size = 5
n_colonies = 20
n_loci = 1000
n_larvae = 25
recomb_freq = 0.01 # avg fraction of alleles recombining
n_time_steps = 3000

max_fraction_of_larvae_spB = 0.05 # max freq. of spB larvae at time 0
rate_fraction_of_larvae_spB_decline = 0.001 # set to small number for uniform prob of spB larvae
replace_no_going_up=0.005

w_dir = "/Users/danielesilvestro/Software/micro2macroEvolution"

rseed = np.random.randint(100)
random.seed(rseed)
np.random.seed(rseed)
print("random seed:", rseed)


"""
TO DO:

amount of admixture could affect death of adults (worst being 50/50 A and B) <- Beta distribution

prob of a larva from the pool to enter anemonis could be a function of admixture as well

DONE: plot individuals by colony: eg shape by hierachy, Y-axis by nDNA admixture, color by mtDNA

DONE: fraction_of_larvae_spB should be a fraction of the amount of larvae in the pool (not a probability)

DONE: Move recombination to within female/male

DONE: make hybridation rate decrease through time?

DONE plot n. pure A individuals, n. admixt and n. pure B nDNA
DONE plot n. pure A individuals and n. pure B mtDNA

DONE plot each individual as a proportion of A vs B
color circles based on mtDNA (A: blue, B: red)

	
"""





# probabilities of death 
# right now death is a function of hierarchy, not actual age
# at each generation one individual from one colony dies
p_death = (1/(np.linspace(1,10,colony_size)))**1 # increase exponent to skew more towards hierarchy =0
p_death = p_death/sum(p_death)



# indexes identifying individuals across all colonies
individual_counter = np.arange(n_colonies*colony_size)
total_pop_size = float(n_colonies*colony_size)

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
# nDNA.shape = 2 x 50 x 1000
# diploid x (n_colonies*n_individuals) x n_alleles


mtDNA =  np.random.choice(range(2),size=len(species_ind))

# hierarchy
hierarchy = np.array(list(range(colony_size))*n_colonies)

# Ages
ages = np.array(list(range(colony_size))*n_colonies)[::-1]


spB_nucleus = []
spB_mitoc = []
nDNA_indApure = []
nDNA_indBpure = []


for time_step in range(n_time_steps):
	# run reproduction
	females = individual_counter[hierarchy==0]
	males   = individual_counter[hierarchy==1]
	
	
	# recombination
	n_recombining =  np.random.poisson(recomb_freq*n_loci)
	
	# indexes of recombining alleles (females)
	xa = np.random.randint(0,n_loci,(n_recombining, n_colonies))
	# indexes of recombining alleles (males)
	xb = np.random.randint(0,n_loci,(n_recombining, n_colonies))
	
	nDNA_recomb = nDNA + 0
	nDNA_recomb[0,females,xa] = nDNA[1,females,xa] + 0
	nDNA_recomb[1,females,xa] = nDNA[0,females,xa] + 0
	nDNA_recomb[0,males,xb] = nDNA[1,males,xb] + 0
	nDNA_recomb[1,males,xb] = nDNA[0,males,xb] + 0
	
	# choose the random alleles passed on to each larva by f/m (only for nDNA)
	allele_ind_female = np.random.choice([0,1],(n_larvae,n_colonies))
	allele_ind_male   = np.random.choice([0,1],(n_larvae,n_colonies))
	# allele_ind[0,:] : which allele in each larvae in 1st colony
	larvae_pool_nDNA  = np.zeros((n_colonies,2,n_larvae,n_loci))
	larvae_pool_mtDNA = np.zeros((n_colonies,n_larvae))
	for j in range(n_colonies):
		female_allele_larvae   = nDNA_recomb[allele_ind_female[:,j],females[j],:]    + 0 # (n_larvae_per_colony x n_loci)
		male_allele_larvae     = nDNA_recomb[allele_ind_male[:,j],males[j],:]        + 0 
		larvae_pool_nDNA[j]    = np.array([female_allele_larvae,male_allele_larvae]) + 0 
		larvae_pool_mtDNA[j,:] = mtDNA[females[j]]                                   + 0 
	
	# make up a fraction of spB larvae
	r = rate_fraction_of_larvae_spB_decline * np.exp(-rate_fraction_of_larvae_spB_decline*time_step)
	fraction_of_larvae_spB = (max_fraction_of_larvae_spB*r) / rate_fraction_of_larvae_spB_decline
	
	
	n_spB_larvae = np.random.poisson(fraction_of_larvae_spB*(n_colonies*n_larvae)) # could replace with random Poisson
	if n_spB_larvae > 0:
		nDNA_spB  = np.random.choice(np.arange(2,4),(n_spB_larvae,2,n_loci))
		mtDNA_spB = np.random.choice(np.arange(2,4),(n_spB_larvae))
		x = np.random.randint(0,n_colonies,n_spB_larvae)
		y = np.random.randint(0,n_larvae,  n_spB_larvae)
		larvae_pool_nDNA[x,:,y,:]  = nDNA_spB
		larvae_pool_mtDNA[x,y]     = mtDNA_spB
	
	
	def returnLarvae(larvae_pool_nDNA,larvae_pool_mtDNA,n=1):
		x = np.random.randint(0,n_colonies,n)
		y = np.random.randint(0,n_larvae,  n)
		larvae_nDNA = larvae_pool_nDNA[x,:,y,:]
		larvae_mtDNA = larvae_pool_mtDNA[x,y]
		return( [larvae_nDNA, larvae_mtDNA] )
	
	
	if np.random.random() < 1:
		# kill adult individuals
		rnd_colony = np.random.choice(range(n_colonies))
		rnd_individual = np.random.choice(range(colony_size), p=p_death)
		
		affected_colony_nDNA = nDNA[:,colony_ind==rnd_colony,:]+0
		affected_colony_mtDNA = mtDNA[colony_ind==rnd_colony]+0
				
		if np.random.random() < replace_no_going_up:
			# fill the gap without scaling up individuals in the colony
			larva = returnLarvae(larvae_pool_nDNA,larvae_pool_mtDNA,n=1)
			affected_colony_nDNA[:,rnd_individual,:] = larva[0]
			affected_colony_mtDNA[rnd_individual]    = larva[1]
			#replace_no_going_up = 0
		else:
			# move up individuals in the colony
			ind_going_up = np.arange(rnd_individual+1, colony_size)
			affected_colony_nDNA[:,ind_going_up-1,:] = affected_colony_nDNA[:,ind_going_up,:]	
			affected_colony_mtDNA[ind_going_up-1] = affected_colony_mtDNA[ind_going_up]		
			# pick a larva
			larva = returnLarvae(larvae_pool_nDNA,larvae_pool_mtDNA,n=1)
			affected_colony_nDNA[:,colony_size-1,:] = larva[0] 
			affected_colony_mtDNA[colony_size-1]    = larva[1] 
	
		# reset global DNA variables
		#print(affected_colony_mtDNA)
		nDNA[:,colony_ind==rnd_colony,:] = affected_colony_nDNA
		mtDNA[colony_ind==rnd_colony] = affected_colony_mtDNA

	spB_nucleus.append( size(nDNA[nDNA<=1])/float(np.size(nDNA)))
	spB_mitoc.append(size(mtDNA[mtDNA<=1])/float(np.size(mtDNA)))
		
	# sp_indx1 = np.amax(nDNA,axis=(0,2)) # if sp_indx1 >= 2: spB is present; sp_indx1 <= 1: spB is absent
	# sp_indx2 = np.amin(nDNA,axis=(0,2)) # if sp_indx1 >= 2: spA is absent;  sp_indx1 <= 1: spA is present
	# 
	# sp_indx1[sp_indx1<=1] = 0
	# sp_indx1[sp_indx1> 1] = 1
	# sp_indx2[sp_indx2<=1] = 0
	# sp_indx2[sp_indx2> 1] = 1
	# nDNAadmixture =  sp_indx1+sp_indx2 # nDNAadmixture = 0 only spA, =1 admix, =2 only spB
	
	nDNAadmixture = nDNA+0
	nDNAadmixture[nDNAadmixture<=1] = 0
	nDNAadmixture[nDNAadmixture> 1] = 1
	# print(np.shape(nDNAadmixture))
	nDNAadmixture = np.mean(nDNAadmixture,axis=(0,2))
	# print(len(nDNAadmixture))
	# quit()
	# nDNAadmixture ==0 only spA, >0 admix, ==1 only spB
	#mtDNA
	
	nDNA_indApure.append(len(nDNAadmixture[nDNAadmixture==0])/total_pop_size)
	nDNA_indBpure.append(len(nDNAadmixture[nDNAadmixture==1])/total_pop_size)
	
	
	if time_step % 100 ==0:  
		if verbose:
			# print sum of full genome (both alleles) per individual : changes only with death/replacement
			print("\nvalue per indvidual (nDNA):", np.sum(np.sum(nDNA,2),0))
			#  print sum of femal genome per individual : changes only with death/replacement
			print("value per indvidual (nDNA,female):", np.sum(nDNA,2)[0])
			# print sum of female alleles across larvae: this changes also with just recombination 
			# without death/replacement
			print(np.sum(recomb_larvae_pool_nDNA[:,0],2).flatten()[0:15])
		
		str_summary = "%s nDNA species A: %s, mtDNA: %s (freq. B larvae: %s)" \
		% (time_step,size(nDNA[nDNA<=1])/float(np.size(nDNA)),size(mtDNA[mtDNA<=1])/float(np.size(mtDNA)),fraction_of_larvae_spB)
		print(str_summary)

fig = plt.figure(figsize=(12, 8))
Ylim = (-0.05,1.05)
fig.add_subplot(221)
plt.plot(range(n_time_steps),spB_nucleus)
plt.gca().set_ylim(Ylim)
plt.gca().set_title('Overall fraction of spA nDNA')

plt.plot(range(n_time_steps),spB_mitoc)
plt.gca().set_ylim(Ylim)
plt.gca().set_title('Overall fraction of spA mtDNA')

fig.add_subplot(222)
plt.plot(range(n_time_steps),nDNA_indApure,color= "blue")
plt.plot(range(n_time_steps),nDNA_indBpure,color= "orange")
plt.plot(range(n_time_steps),spB_mitoc,color= "green")
plt.gca().set_ylim(Ylim)
plt.gca().set_title('spA nDNA (blue), spA mtDNA (green), spB nDNA (orange)')


fig.add_subplot(223)
mtDNA[mtDNA<=1] = 0
mtDNA[mtDNA>1] = 1
colors = np.array(["blue","red"])
plt.scatter(np.arange(len(nDNAadmixture)),nDNAadmixture,marker="o",color=colors[mtDNA]) 
#plt.plot(nDNAadmixture[mtDNA==1],"bo",color="red") 
plt.gca().set_ylim(Ylim)
plt.gca().set_title('0: pure spA nDNA, 1: pure spB nDNA')
	
legend = """
A) overall fraction of spA nDNA (blue) and overall fraction of spA mtDNA (green)
B) fraction of individuals with pure spA nDNA (blue)
   fraction of individuals with pure spA mtDNA (green)
   fraction of individuals with pure spB nDNA (orange)
C) final individuals (0: pure spA nDNA, 1: pure spB nDNA)
   blue circles: spA mtDNA; red: spB mtDNA
D) final individuals by colony
   0: pure spA nDNA, 1: pure spB nDNA
   blue: spA mtDNA; red: spB mtDNA
   '+' female, '*' male, 'o' others 
"""

fig.add_subplot(224)
colors = np.array(["blue","red"])
shapes = np.array(["o","+","*"])

hierarchy[hierarchy>1] = 2

plt.scatter(colony_ind[hierarchy==0], nDNAadmixture[hierarchy==0],marker="+",color=colors[mtDNA[hierarchy==0]], alpha=1) 
plt.scatter(colony_ind[hierarchy==1], nDNAadmixture[hierarchy==1],marker="*",color=colors[mtDNA[hierarchy==1]], alpha=1) 
plt.scatter(colony_ind[hierarchy==2], nDNAadmixture[hierarchy==2],marker="o",color=colors[mtDNA[hierarchy==2]], alpha=0.3) 


#plt.scatter(colony_ind[mtDNA==0], nDNAadmixture[mtDNA==0],marker=(5, 0),color="blue") 
#plt.scatter(colony_ind[mtDNA==1], nDNAadmixture[mtDNA==1],marker="+",color="red") 
plt.gca().set_ylim(Ylim)
plt.gca().set_title('0: pure spA nDNA, 1: pure spB nDNA')


print(legend)
file_name = "%s/plot_r%s.pdf" % (w_dir,rseed)
plt.savefig(file_name)
#plt.show()



