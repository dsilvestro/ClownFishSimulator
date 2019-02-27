#!/usr/bin/env python 
from numpy import *
import numpy as np
np.set_printoptions(suppress=True)
np.set_printoptions(precision=3)  
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.axes
import matplotlib.backends.backend_pdf
import argparse, sys, os
import scipy.stats

p = argparse.ArgumentParser() #description='<input file>') 
p.add_argument('-out',    type=int, help='out name', default= 1, metavar= 1)
args = p.parse_args()


w_dir = max(os.path.dirname(sys.argv[0]), os.getcwd())
logfile_name = "%s/abc_fixSteepnessR%s.log" % (w_dir,args.out)
mcmc_logfile = open(logfile_name , "w",0) 
head = ["it","accepted","max_fraction_of_larvae_spB","mid_point","steepness","max_death_rate_hybrid","rseed"]
#  fraction of spB mtDNA, median admixture (>80% spA), admixt individuals (>80% spA) with spB mtDNA
for i in range(10):
	head.append("spBmtDNA_%s\tspAadmix_%s\tmixAB_%s" % (i,i,i) )
mcmc_logfile.write('\t'.join(head)+'\n')

plot = 0 # set to 1 to run the PDF plots
abc_iteration=0
verbose = 0 
target_samples = 1000
fix_seed = 0
abc_iter = 0

while abc_iteration < target_samples:
	if fix_seed:
		rseed = np.random.randint(0,100000)
		random.seed(rseed)
		np.random.seed(rseed)
	else: rseed = abc_iter
	
	colony_size = 5
	n_colonies = 20
	n_loci = 1000
	n_larvae = 25
	recomb_freq = 0.03 # avg fraction of alleles recombining
	n_time_steps = 3000

	max_fraction_of_larvae_spB = np.random.uniform(0,1) # max influx of spB larvae (fraction)
	mid_point = np.random.uniform(0,n_time_steps)  # mid point of logistic decrease in influx of spB larvae
	steepness = 1 #np.random.uniform(0,1)  # steepness of logistic decrease in influx of spB larvae
	
	replace_no_going_up= 0 #np.random.uniform(0,0.1) # frequency of larvae filling empty space in the colony ladder

	fraction_of_individuals_dying = 0.01 # avg fraction of individuals dying at each generation

	print "\nrandom seed:", rseed, mid_point, steepness
	
	max_death_rate_hybrid = np.random.uniform(0,1)
	
	"""
	TO DO:
	kill larvae based on Beta distribution based on nDNA
	
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
	
	
	
	
	
	# relative probabilities of death 
	# right now death is a function of hierarchy, not actual age
	# at each generation one individual from one colony dies
	p_death = (1/(np.linspace(1,10,colony_size)))**1 # increase exponent to skew more towards hierarchy =0
	p_death = p_death/sum(p_death)
	
	prm_logi = [mid_point,steepness] # mid point and steepness logistic decline in hybridization
	
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

	def get_plot(args):	
		[spB_nucleus, spB_mitoc, nDNA_indApure, nDNA_indBpure, hybridization_list, mtDNA, nDNAadmixture] = args
	
		x_axis = range(len(spB_nucleus))
	
		fig = plt.figure(figsize=(12, 8))
		Ylim = (-0.05,1.05)
		fig.add_subplot(221)
		plt.plot(x_axis,spB_nucleus)
		plt.plot(np.array(hybridization_list),np.zeros(len(hybridization_list)),".",color="black", alpha=0.3)
		plt.gca().set_ylim(Ylim)
		plt.gca().set_title('Overall fraction of spA nDNA')


		plt.plot(x_axis,spB_mitoc)
		plt.gca().set_ylim(Ylim)
		plt.gca().set_title('Overall fraction of spA mtDNA')

		fig.add_subplot(222)
		plt.plot(x_axis,nDNA_indApure,color= "blue")
		plt.plot(x_axis,nDNA_indBpure,color= "orange")
		plt.plot(x_axis,spB_mitoc,color= "green")
		plt.gca().set_ylim(Ylim)
		plt.gca().set_title('spA nDNA (blue), spA mtDNA (green), spB nDNA (orange)')


		fig.add_subplot(223)
		mtDNA_temp = mtDNA+0
		mtDNA_temp[mtDNA_temp<=1] = 0
		mtDNA_temp[mtDNA_temp>1] = 1
		colors = np.array(["blue","red"])
		plt.scatter(np.arange(len(nDNAadmixture)),nDNAadmixture,marker="o",color=colors[mtDNA_temp]) 
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
	
		hierarchy_temp = hierarchy+0
		hierarchy_temp[hierarchy_temp>1] = 2

		plt.scatter(colony_ind[hierarchy_temp==0], nDNAadmixture[hierarchy_temp==0],marker="+",color=colors[mtDNA_temp[hierarchy_temp==0]], alpha=1) 
		plt.scatter(colony_ind[hierarchy_temp==1], nDNAadmixture[hierarchy_temp==1],marker="*",color=colors[mtDNA_temp[hierarchy_temp==1]], alpha=1) 
		plt.scatter(colony_ind[hierarchy_temp==2], nDNAadmixture[hierarchy_temp==2],marker="o",color=colors[mtDNA_temp[hierarchy_temp==2]], alpha=0.3) 

		#plt.scatter(colony_ind[mtDNA==0], nDNAadmixture[mtDNA==0],marker=(5, 0),color="blue") 
		#plt.scatter(colony_ind[mtDNA==1], nDNAadmixture[mtDNA==1],marker="+",color="red") 
		plt.gca().set_ylim(Ylim)
		plt.gca().set_title('0: pure spA nDNA, 1: pure spB nDNA')
	
		return( fig, legend )

	def get_logistic(r0,prm,x):
		# r0 is the max value
		x0,k = prm # mid point and steepness
		kx = min(k * (x-x0), 50)
		rate_at_trait = r0 / ( 1. + exp( kx )    )
		return rate_at_trait

	def calc_death_prob(admix,multi,min_d=0.001):
		# if admix ==0.5 (max admxture) -> delta = multi*min_d
		# if admix ==0   (min admxture) -> delta = min_d
		delta = (admix*2)*multi + min_d
		return(delta)

	spB_nucleus = []
	spB_mitoc = []
	nDNA_indApure = []
	nDNA_indBpure = []
	hybridization_list = []
	figures_list = []
	summary_stats = []
	accepted = 0 

	for time_step in range(n_time_steps+1):
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
	
		# Logistic decline in hybridization rate
		fraction_of_larvae_spB = get_logistic(max_fraction_of_larvae_spB,prm_logi,time_step)
	
		# Exponential decline in hybridization rate
		#r = rate_fraction_of_larvae_spB_decline * np.exp(-rate_fraction_of_larvae_spB_decline*time_step)
		#fraction_of_larvae_spB = (max_fraction_of_larvae_spB*r) / rate_fraction_of_larvae_spB_decline
	
		# hybridization if n_spB_larvae > 0
		n_spB_larvae = np.random.poisson(fraction_of_larvae_spB*(n_colonies*n_larvae)) 
		if n_spB_larvae > 0:
			nDNA_spB  = np.random.choice(np.arange(2,4),(n_spB_larvae,2,n_loci))
			mtDNA_spB = np.random.choice(np.arange(2,4),(n_spB_larvae))
			x = np.random.randint(0,n_colonies,n_spB_larvae)
			y = np.random.randint(0,n_larvae,  n_spB_larvae)
			larvae_pool_nDNA[x,:,y,:]  = nDNA_spB
			larvae_pool_mtDNA[x,y]     = mtDNA_spB
			hybridization_list.append(time_step)
		
		def returnLarvae(larvae_pool_nDNA,larvae_pool_mtDNA, max_death_rate_hybrid,n=1):
			while True:
				x = np.random.randint(0,n_colonies,n)
				y = np.random.randint(0,n_larvae,  n)
				larvae_nDNA = larvae_pool_nDNA[x,:,y,:]
				larvae_mtDNA = larvae_pool_mtDNA[x,y]
				frac_a = np.size(larvae_nDNA[larvae_nDNA<=1])/np.float(np.size(larvae_nDNA))
				admix = min(frac_a, 1-frac_a)
				delta_prob = calc_death_prob(admix,max_death_rate_hybrid,min_d=0.001)
				if np.random.random() > delta_prob:
					break
			#print delta_prob, admix, frac_a
			return( [larvae_nDNA, larvae_mtDNA] )
	
	
		rr = np.random.random(int(total_pop_size))
		n_dead_ind = len(rr[rr<fraction_of_individuals_dying])
	
		if n_dead_ind >0:
			for i in range(n_dead_ind):
				# kill adult individuals
				rnd_colony = np.random.choice(range(n_colonies))
				rnd_individual = np.random.choice(range(colony_size), p=p_death)
	
				affected_colony_nDNA = nDNA[:,colony_ind==rnd_colony,:]+0
				affected_colony_mtDNA = mtDNA[colony_ind==rnd_colony]+0
			
				# pick a larva
				larva = returnLarvae(larvae_pool_nDNA,larvae_pool_mtDNA,max_death_rate_hybrid,n=1)
				#if larva[1]>1: hybridization_list.append(time_step)

				if np.random.random() < replace_no_going_up:
					# fill the gap without scaling up individuals in the colony
					affected_colony_nDNA[:,rnd_individual,:] = larva[0]
					affected_colony_mtDNA[rnd_individual]    = larva[1]
					#replace_no_going_up = 0
				else:
					# move up individuals in the colony
					ind_going_up = np.arange(rnd_individual+1, colony_size)
					affected_colony_nDNA[:,ind_going_up-1,:] = affected_colony_nDNA[:,ind_going_up,:]	
					affected_colony_mtDNA[ind_going_up-1] = affected_colony_mtDNA[ind_going_up]		
					affected_colony_nDNA[:,colony_size-1,:] = larva[0] 
					affected_colony_mtDNA[colony_size-1]    = larva[1] 
			
			# reset global DNA variables
			nDNA[:,colony_ind==rnd_colony,:] = affected_colony_nDNA
			mtDNA[colony_ind==rnd_colony] = affected_colony_mtDNA

		spB_nucleus.append( size(nDNA[nDNA<=1])/float(np.size(nDNA)))
		spB_mitoc.append(size(mtDNA[mtDNA<=1])/float(np.size(mtDNA)))
		
		nDNAadmixture = nDNA+0
		nDNAadmixture[nDNAadmixture<=1] = 0
		nDNAadmixture[nDNAadmixture> 1] = 1
		# print(np.shape(nDNAadmixture))
		nDNAadmixture = np.mean(nDNAadmixture,axis=(0,2))
	
		nDNA_indApure.append(len(nDNAadmixture[nDNAadmixture==0])/total_pop_size)
		nDNA_indBpure.append(len(nDNAadmixture[nDNAadmixture==1])/total_pop_size)
	
	
		# check acceptance
		if time_step % 1 ==0:  
			mtDNA_temp = mtDNA+0
			mtDNA_temp[mtDNA_temp<=1] = 0
			mtDNA_temp[mtDNA_temp>1] = 1
			if time_step % 300  ==0 and time_step>0: 
				try: summary_stats = summary_stats + [np.sum(mtDNA_temp)/total_pop_size,np.median(nDNAadmixture),np.min(nDNAadmixture[mtDNA_temp==1])]
				except: summary_stats = summary_stats + [np.sum(mtDNA_temp)/total_pop_size,np.median(nDNAadmixture),"NA"]
			try:
				#print(np.sum(mtDNA_temp),np.median(nDNAadmixture), np.min(nDNAadmixture[mtDNA_temp==0]) )
				#  fraction of spB mtDNA                        median admixture (>80% spA)       admix individuals (>80% spA) with spB mtDNA
				if np.sum(mtDNA_temp)/total_pop_size > 0.05 and np.median(nDNAadmixture)< 0.2 and np.min(nDNAadmixture[mtDNA_temp==1])<0.2:
					accepted+=1.
			except:
				pass
	
		if time_step % 1000 ==0:  
			if verbose:
				# print sum of full genome (both alleles) per individual : changes only with death/replacement
				print("\nvalue per indvidual (nDNA):", np.sum(np.sum(nDNA,2),0))
				#  print sum of femal genome per individual : changes only with death/replacement
				print("value per indvidual (nDNA,female):", np.sum(nDNA,2)[0])
				# print sum of female alleles across larvae: this changes also with just recombination 
				# without death/replacement
				print(np.sum(recomb_larvae_pool_nDNA[:,0],2).flatten()[0:15])
		
			str_summary = "%s nDNA species A: %s, mtDNA: %s (freq. B larvae: %s) %s" \
			% (time_step,size(nDNA[nDNA<=1])/float(np.size(nDNA)),size(mtDNA[mtDNA<=1])/float(np.size(mtDNA)),fraction_of_larvae_spB,accepted)
			print(str_summary)
			#_ if time_step>1: 
			#_ 	try: print(np.array([np.sum(mtDNA_temp)/total_pop_size,np.median(nDNAadmixture),np.min(nDNAadmixture[mtDNA_temp==1])]))
			#_ 	except: print(np.array([np.sum(mtDNA_temp)/total_pop_size,np.median(nDNAadmixture),np.nan]))
	
		if plot:
			if time_step % 500 ==0 and time_step>0 and accepted >= 0:
				fig, legend = get_plot([spB_nucleus, spB_mitoc, nDNA_indApure, nDNA_indBpure, hybridization_list, mtDNA, nDNAadmixture])
				figures_list.append(fig)
	
	#print(summary_stats)
	
	if plot:
		file_name = "%s/plots/plot_r%s.pdf" % (w_dir,rseed)
		pdf = matplotlib.backends.backend_pdf.PdfPages(file_name)
		
		for fig in figures_list: ## will open an empty extra figure :(
		    pdf.savefig( fig )
		pdf.close()
		print(legend)
	

	if accepted>0: 
		log_state = map(str, (abc_iteration,accepted/n_time_steps, max_fraction_of_larvae_spB,mid_point,steepness,max_death_rate_hybrid,rseed )) + map(str,summary_stats)
		mcmc_logfile.write('\t'.join(log_state)+'\n')
		mcmc_logfile.flush()
		abc_iteration+=1
	
	abc_iter+=1


	