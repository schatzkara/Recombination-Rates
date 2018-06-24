#! python3

# script to simulate DNA sequences mutating and produce a matrix of the ID% between them

import random
import numpy as np

def sim_id_percent(L, generations, GC_prop, kappa, phi):
	alpha = ((1/L) * kappa)/(kappa + 1) # probability of transitions
	beta = (1/L)/(kappa + 1) # probability of transversions
	ancestor = L*[None] # ancestor DNA strand
	strains = 2*[None] # list of the original DNA sequences
	# cms = generations*[None] # list of the number of convergent mutations that have arisen up to this point with the index = (number of SNPs - 1) aka (number of generations - 1)
	id_cms = {}
	nucleotides = ['A', 'T', 'G', 'C']
	# the following are the probabilites of mutating to each of the nucleotides depending on the original nucleotide
	weights_A = [0, (beta*phi), alpha, (beta*(1-phi))]
	weights_T = [(beta*phi), 0, (beta*(1-phi)), alpha]
	weights_G = [alpha, (beta*(1-phi)), 0, (beta*phi)]
	weights_C = [(beta*(1-phi)), alpha, (beta*phi), 0]

	# creates the initial strains
	# time comlexity: O(n), where n is L
	# AT = 100
	AT = int(L*(1-GC_prop)) # the number of As and Ts that will be present in the ancestor
	GC = int(L*(GC_prop)) # the number of Gs and Cs that will be present in the ancestor

	# generates the ancestral strain according to the number of As, Ts, Gs, and Cs necessary
	for x in range(AT):
		ancestor[x] = random.choice(['A', 'T'])
	for y in range(GC): 
		ancestor[AT + y] = random.choice(['G', 'C'])
	random.shuffle(ancestor) # randomly organizes the nucleotides to produce a random ancestor with the correct GC%
	ancestor = ''.join(ancestor)
	for z in range(2): # duplicates the ancestor to generate 2 strains
		strains[z] = ancestor

	# adds the SNPs to each strain
	mutation_sites = [[],[]] # a list of lists where each list contains the sites that were mutated in the corresponding strand; the strand number is the first index and the generation number in the second index
	c = 0 # counter for the number of convergent mutations that have arisen
	for gen in range(generations): # adds one SNP for each generation
		for child in range(2): # adds one SNP to each child strand
			t = list(strains[child])
			site = random.randint(0,L-1) # the site that is to be mutated
			# while(site in mutation_sites[child]): # forces a new mutation site
			# 	site = random.randint(0,L-1)
			if(site not in mutation_sites[child]):
				mutation_sites[child].append(site)
			current = t[site]
			# mutates the nucleotide based on the appropriate probabilties
			if current == 'A':
				t[site] = (random.choices(nucleotides, weights=weights_A, k=1))[0]
			elif current == 'T':
				t[site] = (random.choices(nucleotides, weights=weights_T, k=1))[0]
			elif current == 'G':
				t[site] = (random.choices(nucleotides, weights=weights_G, k=1))[0]
			elif current == 'C':
				t[site] = (random.choices(nucleotides, weights=weights_C, k=1))[0]
			t = ''.join(t)
			strains[child] = t # replaces the old child with the new one
		for site in mutation_sites[0]: # counts up the number of convergent mutations
			if(strains[0][site] == strains[1][site]):
				c += 1
		identity = (L - len(mutation_sites[0]) - len(mutation_sites[1]) + 2*c)/L
		id_cms[gen] = [identity, c]
		c = 0
	# print(len(id_cms.keys()))
	# print(len(id_cms.keys()))
	return id_cms

def id_matrix_sim(n, L, generations, mu, kappa, phi):
	alpha = (mu * kappa)/(kappa + 1) # probability of transitions
	beta = mu/(kappa + 1) # probability of transversions
	ancestor = '' # ancestor DNA strand
	strains = n*[None] # list of the child DNA sequences
	nucleotides = ['A', 'T', 'G', 'C']
	id_matrix = np.empty([n,n], dtype = float, order='C')
	cms_list_matrix = np.empty([n,n], dtype = {}, order='C')
	# the following are the probabilites of mutating to each of the nucleotides or staying the same depending on the original nucleotide
	weights_A = [(1-mu), (beta*phi), alpha, (beta*(1-phi))]
	weights_T = [(beta*phi), (1-mu), (beta*(1-phi)), alpha]
	weights_G = [alpha, (beta*(1-phi)), (1-mu), (beta*phi)]
	weights_C = [(beta*(1-phi)), alpha, (beta*phi), (1-mu)]

	for x in range(L): # builds a random ancestor strand of length L
		nucleotide = random.choice(nucleotides)
		ancestor += nucleotide
	for y in range(n):
		strains[y] = ancestor

	mutation_sites = [n*[None]]
	for gen in range(generations): # adds one SNP for each generation
		for child in range(n): # adds one SNP to each child strand
			t = list(strains[child])
			site = random.randint(0,L-1) # the site that is to be mutated
			if(site not in mutation_sites[child]):
				mutation_sites[child].append(site)
			current = t[site]
			# mutates the nucleotide based on the appropriate probabilties
			if current == 'A':
				t[site] = (random.choices(nucleotides, weights=weights_A, k=1))[0]
			elif current == 'T':
				t[site] = (random.choices(nucleotides, weights=weights_T, k=1))[0]
			elif current == 'G':
				t[site] = (random.choices(nucleotides, weights=weights_G, k=1))[0]
			elif current == 'C':
				t[site] = (random.choices(nucleotides, weights=weights_C, k=1))[0]
			t = ''.join(t)
			strains[child] = t # replaces the old child with the new one
		c = 0
		for strain1 in range(n):
			for strain2 in range(1,n):
				for site in mutation_sites[strain1]:
					if(strains[strain1][site] == strains[strain2][site]):
						c += 1
						cms_list_matrix[strain1,strain2][site] = strains[strain1][site]
						cms_list_matrix[strain2,strain1][site] = strains[strain1][site]
					identity = (L - len(mutation_sites[strain1]) - len(mutation_sites[strain2]) + 2*c)/L
					id_matrix[strain1,strain2] = identity
					id_matrix[strain2,strain1] = identity
				c = 0
	return id_matrix

id_matrix_sim(4,10,4,1/10,1,1/2)