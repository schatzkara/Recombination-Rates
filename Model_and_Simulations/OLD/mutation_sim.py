#! python3

# script to simulate DNA sequences mutating and detect the number of convergent mutations between them

import random
import numpy as np

# this sim starts with all As
# assumes equal probabilities of C,T,G, forces a mutation, and does not allow for multiple mutations in one site
# params: 
# 	n (int) = number of DNA strands
# 	L (int) = length of the DNA strands
# 	mu (float) = mutation rate (units: mutations per base pair per generation
# return: int that equals the total number of convergent mutations between all pairs of DNA strands
# time comlexity: O(n^3), where n is L
def sim_equal_mutations(n, L, mu):
	totals = {}
	convergent_mutations = 0
	strains = []
	new = []
	nucleotides = ['A', 'T', 'G', 'C']

	# creates the initial strains
	# time complexity: O(n), where n is the bigger of L and n; technically it's L+n
	s = ''
	for y in range(L):
		s+='A'
	for x in range(n):
		strains.append(s)

	# mutates the strains
	# time complexity: O(n^2), where n is the bigger of n adn mu*L; technically it's n*mu*L
	mutations = int(mu * L)
	for s in strains:
		t = list(s)
		for m in range(mutations):
			f = random.randint(1,(L-1))
			while(t[f] != 'A'):
				f = random.randint(1,(L-1))
			while(t[f] == 'A'):
				t[f] = random.choice(nucleotides)
		t = ''.join(t)
		new.append(t)

	# counts up actual o and c values
	# time complexity: O(n^3), where n is the bigger of n and L; technically it's n^2 * L
	o = 0
	c = 0
	for a in range(len(new)):
		strain1 = list(new[a])
		for b in range((a+1),len(new)):
			strain2 = list(new[b])
			for d in range(len(strain1)):
				# counts up the number of overlapping mutation sites
				if (strain1[d] != 'A' and strain2[d] != 'A'):
					o += 1
				# counts up the number of overlapping and matching mutation sites
				if (strain1[d] != 'A' and strain1[d] == strain2[d]):
					c += 1
	totals['overlaps'] = o
	totals['matches'] = c
	convergent_mutations = c

	return convergent_mutations

# simulation to mutate DNA strands and count up the total number of convergent mutations between them 
# starts with identical but random strands; factors in the probabilites of A,T,C,G; allows 'mutation' to itself; does not allow for multiple mutations in one spot
# params: 
# 	n (int) = number of DNA strands
# 	L (int) = length of the DNA strands
# 	mu (float) = mutation rate (units: mutations per base pair per generation)
# 	kappa (float) = ratio of transitions to transversions
# 	phi (float) = probability of transversion to its complementary base pair 
# return: int that equals the total number of convergent mutations between all pairs of DNA strands
# time complexity: O(n^3), where n is L
def sim_unequal_mutations(n, L, mu, kappa, phi):
	alpha = (mu * kappa)/(kappa + 1) # probability of transitions
	beta = mu/(kappa + 1) # probability of transversions
	strains = [] # list of the original DNA sequences
	new = [] # list of the mutated DNA sequences
	nucleotides = ['A', 'T', 'G', 'C']
	# the following are the probabilites of mutating to each of the nucleotides or staying the same depending on the original nucleotide
	weights_A = [(1-mu), (beta*phi), alpha, (beta*(1-phi))]
	weights_T = [(beta*phi), (1-mu), (beta*(1-phi)), alpha]
	weights_G = [alpha, (beta*(1-phi)), (1-mu), (beta*phi)]
	weights_C = [(beta*(1-phi)), alpha, (beta*phi), (1-mu)]

	# creates the initial strains
	# time comlexity: O(n), where n is the bigger of L and n; technically it's L+n
	ancestor = ''
	for y in range(L): # builds a random ancestor strand of length L
		ancestor+=random.choice(nucleotides)
	for x in range(n): # duplicates the ancestor to generate n strains
		strains.append(ancestor)

	# mutates the strains
	# time complexity: O(n^2), where n is the bigger of n and L; technically it's n*L
	for child in strains:
		t = list(child)
		for s in range(len(t)): # goes through each nucleotide and mutates it or keeps it the same
			current = t[s]
			if current == 'A':
				t[s] = (random.choices(nucleotides, weights=weights_A, k=1))[0]
			elif current == 'T':
				t[s] = (random.choices(nucleotides, weights=weights_T, k=1))[0]
			elif current == 'G':
				t[s] = (random.choices(nucleotides, weights=weights_G, k=1))[0]
			elif current == 'C':
				t[s] = (random.choices(nucleotides, weights=weights_C, k=1))[0]
		t = ''.join(t)
		new.append(t)

	# counts up actual o and c values
	# time complexity: O(n^3), where n is the bigger of n and L; technically it's n^2 * L
	o = 0 # counter for number of overlapping sites
	c = 0 # counter for number of overlapping and matching sites aka convergent sites
	totals = {} # dictionary to hold the total overlaps and matches (key: site type, value: number)
	for a in range(len(new)): # uses each strain as 'strain1'
		strain1 = list(new[a]) 
		for b in range((a+1),len(new)): # uses each other strain to which it has not yet been compared as 'strain2'
			strain2 = list(new[b])
			for d in range(len(strain1)):
				# counts up the number of overlapping mutation sites
				if (strain1[d] != ancestor[d] and strain2[d] != ancestor[d]):
					o += 1
					# counts up the number of overlapping and matching mutation sites
					if (strain1[d] == strain2[d]):
						c += 1
	totals['overlaps'] = o
	totals['matches'] = c
	return c

# simulation to mutate 2 DNA strands and count up the total number of convergent mutations between them 
# starts with identical but random strands; factors in the probabilites of A,T,C,G; allows 'mutation' to itself; does not allow for multiple mutations in one spot
# params: 
# 	n (int) = number of DNA strands --- must be 2 for this version of the simulation
# 	L (int) = length of the DNA strands
# 	mu (float) = mutation rate (units: mutations per base pair per generation)
# 	kappa (float) = ratio of transitions to transversions
# 	phi (float) = probability of transversion to its complementary base pair 
# return: int that equals the total number of convergent mutations between all pairs of DNA strands
# time complexity: O(n), where n is L
def efficient_sim_unequal_mutations(n, L, mu, kappa, phi):
	alpha = (mu * kappa)/(kappa + 1) # probability of transitions
	beta = mu/(kappa + 1) # probability of transversions
	ancestor = '' # ancestor DNA strand
	strains = ['',''] # list of the child DNA sequences
	nucleotides = ['A', 'T', 'G', 'C']
	# the following are the probabilites of mutating to each of the nucleotides or staying the same depending on the original nucleotide
	weights_A = [(1-mu), (beta*phi), alpha, (beta*(1-phi))]
	weights_T = [(beta*phi), (1-mu), (beta*(1-phi)), alpha]
	weights_G = [alpha, (beta*(1-phi)), (1-mu), (beta*phi)]
	weights_C = [(beta*(1-phi)), alpha, (beta*phi), (1-mu)]

	c = 0
	for y in range(L): # builds a random ancestor strand of length L
		nucleotide = random.choice(nucleotides)
		ancestor += nucleotide
		for child in range(2):
			if nucleotide == 'A':
				strains[child] += (random.choices(nucleotides, weights=weights_A, k=1))[0]
			elif nucleotide == 'T':
				strains[child] += (random.choices(nucleotides, weights=weights_T, k=1))[0]
			elif nucleotide == 'G':
				strains[child] += (random.choices(nucleotides, weights=weights_G, k=1))[0]
			elif nucleotide == 'C':
				strains[child] += (random.choices(nucleotides, weights=weights_C, k=1))[0]
		if(strains[0][y] == strains[1][y] and strains[0][y] != ancestor[y]): 
			c += 1
	return c
