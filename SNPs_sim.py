#! python3

# script to simulate DNA sequences mutating and detect the number of convergent mutations between them as compared to the number of SNPs

import random
import numpy as np
import csv

# simulation to mutate 2 DNA strands a given number of times and count up the total number of convergent mutations between them
# starts with identical but random strands; factors in the probabilites of A,T,C,G; allows 'mutation' to itself; does not allow for multiple mutations in one spot
# params: 
# 	L (int) = length of the DNA strands
# 	generations (int) = the number DNA replications to simulate; for our purposes, this is the number of SNPs that we want to insert
# 	GC (float) = the proportion of the genome that is comprised of Gs and Cs
# 	kappa (float) = ratio of transitions to transversions
# 	phi (float) = probability of transversion to its complementary base pair 
# return: list that contains the number of convergent mutations that have arisen after introducing mutations equivalent to the index + 1
# time complexity: O()
def sim_snps(L, generations, GC_prop, kappa, phi):
	alpha = ((1/L) * kappa)/(kappa + 1) # probability of transitions
	beta = (1/L)/(kappa + 1) # probability of transversions
	ancestor = L*[None] # ancestor DNA strand
	strains = 2*[None] # list of the original DNA sequences
	new = 2*[None] # list of the mutated DNA sequences
	cms = generations*[None] # list of the number of convergent mutations that have arisen up to this point with the index = (number of SNPs - 1) aka (number of generations - 1)
	nucleotides = ['A', 'T', 'G', 'C']
	# the following are the probabilites of mutating to each of the nucleotides or staying the same depending on the original nucleotide
	weights_A = [0, (beta*phi), alpha, (beta*(1-phi))]
	weights_T = [(beta*phi), 0, (beta*(1-phi)), alpha]
	weights_G = [alpha, (beta*(1-phi)), 0, (beta*phi)]
	weights_C = [(beta*(1-phi)), alpha, (beta*phi), 0]

	# creates the initial strains
	# time comlexity: O(n), where n is the bigger of L and n; technically it's L+n
	AT = int(L*(1-GC_prop))
	GC = int(L*(GC_prop))
	for x in range(AT):
		ancestor[x] = random.choice(['A', 'T'])
	for y in range(GC): 
		ancestor[AT + y] = random.choice(['G', 'C'])
	# ancestor = ''.join(ancestor)
	random.shuffle(ancestor)
	ancestor = ''.join(ancestor)
	# print(ancestor)
	for z in range(2): # duplicates the ancestor to generate n strains
		strains[z] = ancestor

	# print('ancestor: ' + ancestor)
	# print('strains: ')
	# print(strains)

	# adds the SNPs to each strain
	mutation_sites = [generations*[None],generations*[None]] # a list of lists where each list contains the sites that were mutated in the corresponding strand; the strand number is the first index and the gen number in the second index
	c = 0
	for gen in range(generations):
		for child in range(2):
			t = list(strains[child])
			site = random.randint(0,L-1)
			while(site in mutation_sites[child]):
				site = random.randint(0,L-1)
			mutation_sites[child][gen] = site
			current = t[mutation_sites[child][gen]]
			if current == 'A':
				t[mutation_sites[child][gen]] = (random.choices(nucleotides, weights=weights_A, k=1))[0]
			elif current == 'T':
				t[mutation_sites[child][gen]] = (random.choices(nucleotides, weights=weights_T, k=1))[0]
			elif current == 'G':
				t[mutation_sites[child][gen]] = (random.choices(nucleotides, weights=weights_G, k=1))[0]
			elif current == 'C':
				t[mutation_sites[child][gen]] = (random.choices(nucleotides, weights=weights_C, k=1))[0]
			t = ''.join(t)
			strains[child] = t
		for site in mutation_sites[0]:
			if site != None:
				if(list(strains[0])[site] == list(strains[1])[site]):
					c += 1
		cms[gen] = c
		c = 0
	# print('new: ')
	# print(strains)
	# print(cms)

	return cms

Ls = [1000]
generations = 300
GCs = [0.50]
kappas = [1.0]
phi = 1/2
iterations = 1000

for l in Ls: # iterates over every length desired
	for gc in GCs: # iterates over every mu desired
		for k in kappas:
			# data = {} # dictionary for data over all iterations (key: iteration number, value: average convergent mutations of all iterations up to that point)
			for z in range(iterations): # number of iterations to run
				with open(('SNPs_sim_data_' + str(l) + '_' + str(gc) + '_' + str(k) + '_' + str(z+1) + '.csv'), 'w', newline = '') as f: # writes the data to a .csv file
				    writer = csv.writer(f)
				    writer.writerow(['SNPs on Each Strain', 'CMs']) # column headers
				    data = [(list(range(1,generations+1))), (sim_snps(l,generations,gc,k,phi))]
				    data = zip(*data)
				    writer.writerows(data)