#! python3

# script to simulate 2 DNA sequences mutating and detect the number of convergent mutations between them as compared to the number of SNPs

import random
import numpy as np
import csv

# simulation to mutate 2 DNA strands a given number of times and count up the total number of convergent mutations between them
# starts with identical but random strands; factors in the probabilites of A,T,C,G; forces a single mutation with each 'generation,'' does not allow for multiple mutations at one site
# this simulation adds a certain number of SNPs to each DNA strain and counts up the number of convergent mutations that arise after each SNP is added
# params: 
# 	L (int) = length of the DNA strands
# 	generations (int) = the number DNA replications to simulate; for our purposes, this is the number of SNPs that we want to insert
# 	GC (float) = the proportion of the genome that is comprised of Gs and Cs
# 	kappa (float) = ratio of transitions to transversions
# 	phi (float) = probability of transversion to its complementary base pair 
# return: list that contains the number of convergent mutations that have arisen after introducing mutations equivalent to the index + 1
# time complexity: O(n^3), where n is L (really it's L * generations^2)
def sim_snps(L, generations, GC_prop, kappa, phi):
	alpha = ((1/L) * kappa)/(kappa + 1) # probability of transitions
	beta = (1/L)/(kappa + 1) # probability of transversions
	ancestor = L*[None] # ancestor DNA strand
	strains = 2*[None] # list of the original DNA sequences
	cms = generations*[None] # list of the number of convergent mutations that have arisen up to this point with the index = (number of SNPs - 1) aka (number of generations - 1)
	nucleotides = ['A', 'T', 'G', 'C']
	# the following are the probabilites of mutating to each of the nucleotides depending on the original nucleotide
	weights_A = [0, (beta*phi), alpha, (beta*(1-phi))]
	weights_T = [(beta*phi), 0, (beta*(1-phi)), alpha]
	weights_G = [alpha, (beta*(1-phi)), 0, (beta*phi)]
	weights_C = [(beta*(1-phi)), alpha, (beta*phi), 0]

	# creates the initial strains
	# time comlexity: O(n), where n is L
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
	mutation_sites = [generations*[None],generations*[None]] # a list of lists where each list contains the sites that were mutated in the corresponding strand; the strand number is the first index and the generation number in the second index
	c = 0 # counter for the number of convergent mutations that have arisen
	for gen in range(generations): # adds one SNP for each generation
		for child in range(2): # adds one SNP to each child strand
			t = list(strains[child])
			site = random.randint(0,L-1) # the site that is to be mutated
			while(site in mutation_sites[child]): # forces a new mutation site
				site = random.randint(0,L-1)
			mutation_sites[child][gen] = site
			current = t[mutation_sites[child][gen]]
			# mutates the nucleotide based on the appropriate probabilties
			if current == 'A':
				t[mutation_sites[child][gen]] = (random.choices(nucleotides, weights=weights_A, k=1))[0]
			elif current == 'T':
				t[mutation_sites[child][gen]] = (random.choices(nucleotides, weights=weights_T, k=1))[0]
			elif current == 'G':
				t[mutation_sites[child][gen]] = (random.choices(nucleotides, weights=weights_G, k=1))[0]
			elif current == 'C':
				t[mutation_sites[child][gen]] = (random.choices(nucleotides, weights=weights_C, k=1))[0]
			t = ''.join(t)
			strains[child] = t # replaces the old child with the new one
		for site in mutation_sites[0]: # counts up the number of convergent mutations
			if site != None:
				if(list(strains[0])[site] == list(strains[1])[site]):
					c += 1
		cms[gen] = c
		c = 0
		
	return cms