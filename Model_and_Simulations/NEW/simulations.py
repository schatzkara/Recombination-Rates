#! python 3

import random
import numpy as np

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
def mutations_sim(n, L, mu, kappa, phi):
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
def efficient_mutations_sim(n, L, mu, kappa, phi):
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

# simulation to mutate 2 DNA strands a given number of times and count up the total number of convergent mutations between them
# starts with identical but random strands; factors in the probabilites of A,T,C,G; forces a single mutation with each 'generation,' allows for multiple mutations at one site
# this simulation adds a certain number of SNPs to each DNA strain and counts up the number of convergent mutations that arise after each SNP is added
# params: 
# 	L (int) = length of the DNA strands
# 	generations (int) = the number DNA replications to simulate; for our purposes, this is the number of SNPs that we want to insert
# 	GC (float) = the proportion of the genome that is comprised of Gs and Cs
# 	kappa (float) = ratio of transitions to transversions
# 	phi (float) = probability of transversion to its complementary base pair 
# return: list that contains the number of convergent mutations that have arisen after introducing mutations equivalent to the index + 1
# time complexity: O(n^3), where n is L (really it's L * generations^2)
def mutations_over_generations_sim(L, generations, GC_prop, kappa, phi):
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
			if(strains[0][site] == strains[1][site] != ancestor[site]):
				c += 1
		cms[gen] = c
		c = 0

	return cms

def id_percent_sim(L, generations, GC_prop, kappa, phi):
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
	o = 0
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
			if t[site] == ancestor[site]:
				mutation_sites[child].remove(site)
			t = ''.join(t)
			strains[child] = t # replaces the old child with the new one
		for site in mutation_sites[0]: # counts up the number of convergent mutations
			if site in mutation_sites[1]:
				o += 1
				if(strains[0][site] == strains[1][site] != ancestor[site]):
					c += 1
		identity = (L - len(mutation_sites[0]) - len(mutation_sites[1]) + o + c)/L
		id_cms[gen] = [identity, c]
		o = 0
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
	id_matrix = np.empty([n,n], dtype = np.float, order='C')
	# cms_list_matrix = np.empty([n,n], order='C')
	cms_list_matrix = [[{} for x in range(n)] for y in range(n)]
	# print(cms_list_matrix)
	# the following are the probabilites of mutating to each of the nucleotides or staying the same depending on the original nucleotide
	weights_A = [0, (beta*phi), alpha, (beta*(1-phi))]
	weights_T = [(beta*phi), 0, (beta*(1-phi)), alpha]
	weights_G = [alpha, (beta*(1-phi)), 0, (beta*phi)]
	weights_C = [(beta*(1-phi)), alpha, (beta*phi), 0]

	for x in range(L): # builds a random ancestor strand of length L
		nucleotide = random.choice(nucleotides)
		ancestor += nucleotide
	for y in range(n):
		strains[y] = ancestor

	mutation_sites = n*[None]
	for x in range(n):
		mutation_sites[x] = []
	# print(mutation_sites)
	for gen in range(generations): # adds one SNP for each generation
		for child in range(n): # adds one SNP to each child strand
			# print(child)
			t = list(strains[child])
			site = random.randint(0,L-1) # the site that is to be mutated
			# print(mutation_sites[child])
			if(site not in mutation_sites[child]):
				mutation_sites[child].append(site)
				# print(mutation_sites)
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
			if t[site] == ancestor[site]:
				mutation_sites[child].remove(site)
			t = ''.join(t)
			strains[child] = t # replaces the old child with the new one
	o = 0
	c = 0
	# print(ancestor)
	# print(strains)
	# print(mutation_sites)
	for strain1 in range(n):
		id_matrix[strain1][strain1] = 1
		for strain2 in range(strain1+1,n):
			# print(str(strain1) + ',' + str(strain2))
			# cms_list_matrix[strain1,strain2] = {}
			# cms_list_matrix[strain2,strain1] = {}
			for site in mutation_sites[strain1]:
				if site in mutation_sites[strain2]:
					o += 1
					if(strains[strain1][site] == strains[strain2][site] != ancestor[site]):
						c += 1
						cms_list_matrix[strain1][strain2][site] = strains[strain1][site]
						cms_list_matrix[strain2][strain1][site] = strains[strain1][site]
					# cms_list_matrix[strain1,strain2][site] = (int(site), strains[strain1][site])
					# cms_list_matrix[strain2,strain1][site] = (int(site), strains[strain1][site])
			# print(len(mutation_sites[strain1]))
			# print(len(mutation_sites[strain2]))
			identity = (L - len(mutation_sites[strain1]) - len(mutation_sites[strain2]) + o + c)/L
			id_matrix[strain1,strain2] = identity
			id_matrix[strain2,strain1] = identity
			print('1: ' + str(strain1) + ', 2: ' + str(strain2) + ', m1: ' + str(len(mutation_sites[strain1])) + ', m2: ' + str(len(mutation_sites[strain2])) + ', o: ' + str(o) + ', c: ' + str(c))
			o = 0
			c = 0
	# print(cms_list_matrix)
	# print(id_matrix)
	# return id_matrix

	cms = L*[None]
	strains_with_site = L*[None]
	# print(strains_with_site)
	for x in range(L):
		cms[x] = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
		strains_with_site[x] = []
	for strain1 in range(n):
		for strain2 in range(strain1+1, n):
			# print(str(strain1) + ',' + str(strain2))
			keys = list(cms_list_matrix[strain1][strain2].keys())
			for z in range(len(keys)):
				site = keys[z]
				# print(cms)
				nucleotide = cms_list_matrix[strain1][strain2][site]
				site_count = cms[site]
				# print(site_count)
				# print(nucleotide)
				# print(site_count[nucleotide])
				if strain1 not in strains_with_site[site]:
					strains_with_site[site].append(strain1)
					site_count[nucleotide] += 1
				if strain2 not in strains_with_site[site]:
					strains_with_site[site].append(strain2)
					site_count[nucleotide] += 1

	# print(cms)

	multiple_cms = 0
	for site in cms:
		for base in nucleotides:
			if site[base] > 2:
				multiple_cms += 1
	print(multiple_cms)
	return multiple_cms

# id_matrix_sim(10,10,6,1/10,1,1/2)