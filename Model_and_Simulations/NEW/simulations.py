
import random
import numpy as np

# keyword parameters: GC_prop
def generate_ancestor(L, **keyword_parameters):
	nucleotides = ['A', 'T', 'G', 'C']
	if 'GC_prop' in keyword_parameters:
		GC_prop = keyword_parameters['GC_prop']
		ancestor = L*[None] # ancestor DNA strand
		GC = int(L*(GC_prop)) # the number of Gs and Cs that will be present in the ancestor
		AT = int(L*(1-GC_prop)) # the number of As and Ts that will be present in the ancestor
		while AT + GC != L:
			if AT % 2 != 0:
				AT += 1
			elif GC % 2 != 0:
				GC += 1
		# generates the ancestral strand according to the number of As, Ts, Gs, and Cs necessary
		for x in range(AT):
			ancestor[x] = random.choice(['A', 'T'])
		for y in range(GC):
			ancestor[AT + y] = random.choice(['G', 'C'])
		random.shuffle(ancestor) # randomly organizes the nucleotides to produce a random ancestor with the correct GC%
		ancestor = ''.join(ancestor)
	else:
		ancestor = ''
		for y in range(L): # builds a random ancestor strand of length L
			ancestor+=random.choice(nucleotides)
	return ancestor
	

def mutate(current_nucleotide, mu, kappa, phi, force_mutation):
	nucleotides = ['A', 'T', 'G', 'C']
	alpha = (mu * kappa)/(kappa + 1) # probability of transitions
	beta = mu/(kappa + 1) # probability of transversions
	# the following are the probabilites of mutating to each of the nucleotides or staying the same depending on the original nucleotide
	if force_mutation:
		weights_A = [0, (beta*phi)/(alpha+beta), alpha/(alpha+beta), (beta*(1-phi))/(alpha+beta)]
		weights_T = [(beta*phi)/(alpha+beta), 0, (beta*(1-phi))/(alpha+beta), alpha/(alpha+beta)]
		weights_G = [alpha/(alpha+beta), (beta*(1-phi))/(alpha+beta), 0, (beta*phi)/(alpha+beta)]
		weights_C = [(beta*(1-phi))/(alpha+beta), alpha/(alpha+beta), (beta*phi)/(alpha+beta), 0]
	else:
		weights_A = [(1-mu), (beta*phi), alpha, (beta*(1-phi))]
		weights_T = [(beta*phi), (1-mu), (beta*(1-phi)), alpha]
		weights_G = [alpha, (beta*(1-phi)), (1-mu), (beta*phi)]
		weights_C = [(beta*(1-phi)), alpha, (beta*phi), (1-mu)]
	if current_nucleotide == 'A':
		new_nucleotide = (np.random.choice(nucleotides, p=weights_A))#, k=1))[0]
	elif current_nucleotide == 'T':
		new_nucleotide = (np.random.choice(nucleotides, p=weights_T))#, k=1))[0]
	elif current_nucleotide == 'G':
		new_nucleotide = (np.random.choice(nucleotides, p=weights_G))#, k=1))[0]
	elif current_nucleotide == 'C':
		new_nucleotide = (np.random.choice(nucleotides, p=weights_C))#, k=1))[0]
	return str(new_nucleotide)

# keyword parameters: mutation_sites1 and mutation_sites2, return_sites
def detect_o_and_c(ancestor, strain1, strain2, **keyword_parameters):
	o = 0
	c = 0
	sites = {}
	if 'mutation_sites1' and 'mutation_sites2' in keyword_parameters:
		mutation_sites1 = keyword_parameters['mutation_sites1']
		mutation_sites2 = keyword_parameters['mutation_sites2']
		for site in mutation_sites1: # counts up the number of convergent mutations
			if site in mutation_sites2:
				o += 1
				if(strain1[site] == strain2[site] and strain1[site] != ancestor[site]):
					c += 1
					sites[site] = strain1[site]
		if 'return_sites' in keyword_parameters and keyword_parameters['return_sites']:
			return_value = {'o': o, 'c': c, 'sites': sites}
		else:                                        
			return_value = {'o': o, 'c': c}
	
	else:
		for site in range(len(strain1)):
			if (strain1[site] != ancestor[site] and strain2[site] != ancestor[site]):
				o += 1
				if (strain1[site] == strain2[site]):
					c += 1
					sites[site] = strain1[site]
		if 'return_sites' in keyword_parameters and keyword_parameters['return_sites']:
			return_value = {'o': o, 'c': c, 'sites': sites}
		else:
			return_value = {'o': o, 'c': c}
	return return_value
        
# simulation to mutate DNA strands and count up the total number of convergent mutations between them 
# starts with identical but random strands; factors in the probabilites of A,T,C,G; allows 'mutation' to itself; does not allow for multiple mutations in one spot
# params: 
# 	n (int) = number of DNA strands
# 	L (int) = length of the DNA strands
# 	mu (float) = mutation rate (units: mutations per base pair per generation)
#       generations (int) = number of generations to elapse
# 	kappa (float) = ratio of transitions to transversions
# 	phi (float) = probability of transversion to its complementary base pair 
# return: int that equals the total number of convergent mutations between all pairs of DNA strands
# time complexity: O(n^3), where n is L
def mutation_sim(n, L, mu, generations, kappa, phi):
	strains = n*[None] # list of the daughter DNA sequences
	nucleotides = ['A', 'T', 'G', 'C']

	# creates the initial strains
	# time comlexity: O(n), where n is the bigger of L and n; technically it's L+n
	ancestor = generate_ancestor(L)
	for x in range(n): # duplicates the ancestor to generate n strains
		strains[x] = ancestor

	# mutates the strains
	# time complexity: O(n^2), where n is the bigger of n and L; technically it's n*L
	for child in range(n):
		t = list(child)
		for site in range(len(t)): # goes through each nucleotide and mutates it or keeps it the same
			for g in generations:
				current_nucleotide = t[site]
				t[site] = mutate(current_nucleotide, mu, kappa, phi, False)
		t = ''.join(t)
		strains[child] = t # replaces the old strand with the mutated one

	# counts up actual o and c values
	# time complexity: O(n^3), where n is the bigger of n and L; technically it's n^2 * L
	for a in range(len(strains)): # uses each strain as 'strain1'
		strain1 = strains[a] 
		for b in range((a+1),len(new)): # uses each other strain to which it has not yet been compared as 'strain2'
			strain2 = strains[b]
			totals = detect_o_and_c(ancestor, strain1, strain2)
	return totals['c']

# simulation to mutate 2 DNA strands and count up the total number of convergent mutations between them 
# starts with identical but random strands; factors in the probabilites of A,T,C,G; allows 'mutation' to itself; does not allow for multiple mutations in one spot
# params: 
# 	L (int) = length of the DNA strands
# 	mu (float) = mutation rate (units: mutations per base pair per generation)
#       generations (int) = number of generations to elapse
# 	kappa (float) = ratio of transitions to transversions
# 	phi (float) = probability of transversion to its complementary base pair 
# return: int that equals the total number of convergent mutations between all pairs of DNA strands
# time complexity: O(n^2), where n is L
def efficient_mutation_sim(L, mu, generations, kappa, phi):
	ancestor = '' # ancestor DNA strand
	strains = ['',''] # list of the child DNA sequences
	nucleotides = ['A', 'T', 'G', 'C']
	c = 0
	for y in range(L): # builds a random ancestor strand of length L
		nucleotide = random.choice(nucleotides)
		ancestor += nucleotide
		strains[0] += nucleotide
		strains[1] += nucleotide
		for child in range(2):
			t = list(strains[child])
			for g in range(generations):
				current_nucleotide = t[y]
				t[y] = mutate(current_nucleotide, mu, kappa, phi, False)
			t = ''.join(t)
			strains[child] = t
		if(strains[0][y] == strains[1][y] and strains[0][y] != ancestor[y]): 
			c += 1
	return c

# simulation to mutate 2 DNA strands a given number of times and count up the total number of convergent mutations between them
# starts with identical but random strands; factors in the probabilites of A,T,C,G; forces an exact number of mutations with each 'generation,' allows for multiple mutations at one site
# params: 
# 	L (int) = length of the DNA strands
#       mu (float) = mutation rate (units: mutations per base pair per generation)
# 	generations (int) = the number DNA replications to simulate; for our purposes, this is the number of SNPs that we want to insert
# 	GC (float) = the proportion of the genome that is comprised of Gs and Cs
# 	kappa (float) = ratio of transitions to transversions
# 	phi (float) = probability of transversion to its complementary base pair 
# return: list that contains the number of convergent mutations that have arisen after eacg generation
def generation_sim(L, mu, generations, GC_prop, kappa, phi):
	cms = generations*[None] # list of the number of convergent mutations that have arisen up to this point; index = (generation - 1)

	# this section creates the initial strains
	# time comlexity: O(n), where n is L
	ancestor = generate_ancestor(L, GC_prop=GC_prop)
	# duplicates the ancestor strand to generate 2 daughter strands
	strains = 2*[None] # list of the 2 daughter strands
	for z in range(2):
		strains[z] = ancestor

	# this section adds the mutations to each daughter strain
	m = int(mu*L) # the number of mutations to add per generation
	mutation_sites = [[],[]] # a list of lists where each list contains the sites that were mutated in the corresponding strand; the strand number is the first index and the generation number in the second index
	for g in range(generations): # adds m mutations for each generation
		for child in range(2): # adds m mutations to each child strand
			t = list(strains[child])
			for x in range(m): # mutates the appropriate number of sites for a generation
				site = random.randint(0,L-1) # the site that is to be mutated
				if(site not in mutation_sites[child]):
					mutation_sites[child].append(site)
				current_nucleotide = t[site]
				t[site] = mutate(current_nucleotide, mu, kappa, phi, True)
				# mutates the nucleotide based on the appropriate probabilties
			t = ''.join(t)
			strains[child] = t # replaces the old child with the new one

		totals = detect_o_and_c(ancestor, strains[0], strains[1], mutation_sites1 = mutation_sites[0], mutation_sites2 = mutation_sites[1])

		cms[g] = totals['c']
	# print(cms)
	return cms

# simulation to mutate 2 DNA strands a given number of times and calculate the identity percentage and the total number of convergent mutations between them
# starts with identical but random strands; factors in the probabilites of A,T,C,G; forces an exact number of mutations with each 'generation,' allows for multiple mutations at one site
# params: 
# 	L (int) = length of the DNA strands
#       mu (float) = mutation rate (units: mutations per base pair per generation)
# 	generations (int) = the number DNA replications to simulate; for our purposes, this is the number of SNPs that we want to insert
# 	GC (float) = the proportion of the genome that is comprised of Gs and Cs
# 	kappa (float) = ratio of transitions to transversions
# 	phi (float) = probability of transversion to its complementary base pair 
# return: dictionary that contains the identity percentage and number of convergent mutations after each generation
def identity_sim(L, mu, generations, GC_prop, kappa, phi):
	id_cms = {} # a dictionary to hold the ID% and convergent mutations after each generation; key: generation number, value: (ID%, convergent mutations)

	# this section creates the initial strains
	# time comlexity: O(n), where n is L
	ancestor = generate_ancestor(L, GC_prop=GC_prop)
	# duplicates the ancestor to generate 2 daughter strands
	strains = 2*[None] # list of the 2 daughter strands
	for z in range(2):
		strains[z] = ancestor
			
	# this section adds the mutations to each daughter strain
	nucleotides = ['A', 'T', 'G', 'C']
	m = int(mu*L) # the number of mutations to add per generation
	mutation_sites = [[],[]] # a list of lists where each list contains the sites that were mutated in the corresponding strand; the strand number is the first index and the generation number in the second index
	for g in range(generations): # adds m mutations for each generation
		for child in range(2): # adds m mutations to each daughter strand
			t = list(strains[child])
			for x in range(m): # mutates the appropriate number of sites for a generation
				site = random.randint(0,L-1) # the site that is to be mutated
				if(site not in mutation_sites[child]):
					mutation_sites[child].append(site)
				current_nucleotide = t[site]
				# mutates the nucleotide based on the appropriate probabilties
				t[site] = mutate(current_nucleotide, mu, kappa, phi, True)
				if t[site] == ancestor[site]:
					mutation_sites[child].remove(site)
			t = ''.join(t)
			strains[child] = t # replaces the old child with the new one
		totals = detect_o_and_c(ancestor, strains[0], strains[1], mutation_sites1 = mutation_sites[0], mutation_sites2 = mutation_sites[1])
		identity = (L - len(mutation_sites[0]) - len(mutation_sites[1]) + totals['o'] + totals['c'])/L
		id_cms[g+1] = [identity, totals['c']]
	return id_cms
        
def id_matrix_sim(n, L, mu, generations, kappa, phi):
	ancestor = '' # ancestor DNA strand
	strains = n*[None] # list of the child DNA sequences
	id_matrix = np.empty([n,n], dtype = np.float, order='C')
	c_matrix = np.empty([n,n], dtype = np.int64, order='C')
	ancestor = generate_ancestor(L)
	for y in range(n):
		strains[y] = ancestor
	
	# this section adds the mutations to each daughter strain
	nucleotides = ['A', 'T', 'G', 'C']
	m = int(mu*L) # the number of mutations to add per generation
	mutation_sites = n*[None] # a list of lists where each list contains the sites that were mutated in the corresponding strand; the strand number is the first index and the generation number in the second index
	for x in range(n):
		mutation_sites[x] = []
	for g in range(int(generations)): # adds m mutations for each generation
		for child in range(n): # adds m mutations to each daughter strand
			t = list(strains[child])
			for x in range(m): # mutates the appropriate number of sites for a generation
				site = random.randint(0,L-1) # the site that is to be mutated
				if(site not in mutation_sites[child]):
					mutation_sites[child].append(site)
				current_nucleotide = t[site]
				# mutates the nucleotide based on the appropriate probabilties
				t[site] = mutate(current_nucleotide, mu, kappa, phi, True)
				if t[site] == ancestor[site]:
					mutation_sites[child].remove(site)
			t = ''.join(t)
			strains[child] = t # replaces the old child with the new one
				
	for strain1 in range(n):
		id_matrix[strain1][strain1] = 1
		c_matrix[strain1][strain1] = 0
		for strain2 in range(strain1+1,n):
			totals = detect_o_and_c(ancestor, strains[strain1], strains[strain2], mutation_sites1 = mutation_sites[strain1], mutation_sites2 = mutation_sites[strain2])
			identity = (L - len(mutation_sites[strain1]) - len(mutation_sites[strain2]) + totals['o'] + totals['c'])/L
			id_matrix[strain1,strain2] = identity
			id_matrix[strain2,strain1] = identity
			c_matrix[strain1,strain2] = totals['c']
			c_matrix[strain2,strain1] = totals['c']
	return {'id_matrix': id_matrix, 'c_matrix': c_matrix}

# simulation to mutate 2 DNA strands a given number of times and count the number of convergent mutations that are present between 2,3,...,n of them
# starts with identical but random strands; factors in the probabilites of A,T,C,G; forces an exact number of mutations with each 'generation,' allows for multiple mutations at one site
# params: 
#       n (int) = the number of strains in the population
# 	L (int) = length of the DNA strands
#       mu (float) = mutation rate (units: mutations per base pair per generation)
# 	generations (int) = the number DNA replications to simulate; for our purposes, this is the number of SNPs that we want to insert
# 	kappa (float) = ratio of transitions to transversions
# 	phi (float) = probability of transversion to its complementary base pair 
# return: list that contains the number of convergent mutations between each number of them; index = q - 2
def c_q_sim(n, L, mu, generations, kappa, phi):
	strains = n*[None] # list of the child DNA sequences
	id_matrix = np.empty([n,n], dtype = np.float, order='C') # matrix of identity percentages; the (i,j) entry is the ID% between strain i and strain j
	c_list_matrix = [[{} for x in range(n)] for y in range(n)] # matrix of the convergent mutation sites; the (i,j) entry is a dictionary of the convergent mutation sites between strain i and strain j; key = site, value = nucleotide

	ancestor = generate_ancestor(L)
	for y in range(n): # duplicates the ancestor to make n identical strains
		strains[y] = ancestor

	nucleotides = ['A', 'T', 'G', 'C']
	m = int(mu*L) # the number of mutations to add per generation
	mutation_sites = n*[None] # a list of lists where each list contains the sites that were mutated in the corresponding strand; the strand number is the first index and the generation number in the second index
	for x in range(n):
		mutation_sites[x] = []
	for g in range(generations): # adds m mutations for each generation
		for child in range(n): # adds m mutations to each daughter strand
			t = list(strains[child])
			for x in range(m): # mutates the appropriate number of sites for a generation
				site = random.randint(0,L-1) # the site that is to be mutated
				if(site not in mutation_sites[child]):
					mutation_sites[child].append(site)
				current_nucleotide = t[site]
				# mutates the nucleotide based on the appropriate probabilties
				t[site] = mutate(current_nucleotide, mu, kappa, phi, True)
				if t[site] == ancestor[site]: # removes the site from the list of mutation sites if it mutates back to match the ancestral strand
					mutation_sites[child].remove(site)
			t = ''.join(t)
			strains[child] = t # replaces the old child with the new one

	for strain1 in range(n): # allows each strain to be strain1
		id_matrix[strain1][strain1] = 1 
		for strain2 in range(strain1+1,n): # allows each strain after strain1 to be strain2; avoids doing any combination of strains twice
			totals = detect_o_and_c(ancestor, strains[strain1], strains[strain2], mutation_sites1 = mutation_sites[strain1], mutation_sites2 = mutation_sites[strain2], return_sites = True)
			c_list_matrix[strain1][strain2] = totals['sites']
			c_list_matrix[strain2][strain1] = totals['sites']
	# print(c_list_matrix)
	site_counts = L*[None] # list of dictionaries to keep track of which nucleotides are at each convergent site; index = site; key = nucleotide, value = number of strains with that nucleotide
	strains_with_site = L*[None] # list of the strains that have a convergent mutation at each site; index = site
	for x in range(L):
		site_counts[x] = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
		strains_with_site[x] = []
	for strain1 in range(n): # allows each strain to be strain1
		for strain2 in range(strain1+1, n): # allows each strain after strain1 to be strain2; avoids doing any combination of strains twice
			for site in c_list_matrix[strain1][strain2].keys(): # goes through every sites that is convergent between strains 1 and 2
				nucleotide = c_list_matrix[strain1][strain2][site] # extracts the nucleotide that strain1 and strain2 have at the site
				if strain1 not in strains_with_site[site]: # avoids double counting strain1 as convergent at that site
					strains_with_site[site].append(strain1)
					site_counts[site][nucleotide] += 1
				if strain2 not in strains_with_site[site]: # avoids double counting strain2 as convergent at that site
					strains_with_site[site].append(strain2)
					site_counts[site][nucleotide] += 1

	c_q = (n-1)*[None] # list of the number of convergent mutations between q strains; index = q - 2
	for x in range(n-1):
		c_q[x] = 0
	for site in site_counts:
		for base in nucleotides:
			for q in range(2,n+1):
				if site[base] == (q):
					c_q[q-2] += 1
	print(c_q)
	return c_q

def mutation_sites_sim(L, mu, generations, kappa, phi):
	cms = generations*[None] # list of the number of convergent mutations that have arisen up to this point; index = (generation - 1)
	# mutation_sites = generations*[None]
	# this section creates the initial strains
	# time comlexity: O(n), where n is L
	ancestor = generate_ancestor(L)
	daughter = ancestor
	# duplicates the ancestor strand to generate 2 daughter strands
	# strains = 2*[None] # list of the 2 daughter strands
	# for z in range(2):
	# 	strains[z] = ancestor

	# this section adds the mutations to each daughter strain
	m = int(mu*L) # the number of mutations to add per generation
	mutation_sites = [] # a list of lists where each list contains the sites that were mutated in the corresponding strand; the strand number is the first index and the generation number in the second index
	for g in range(generations): # adds m mutations for each generation
		# for child in range(2): # adds m mutations to each child strand
		t = list(daughter) # list(strains[child])
		for x in range(m): # mutates the appropriate number of sites for a generation
			site = random.randint(0,L-1) # the site that is to be mutated
			if(site not in mutation_sites):
				mutation_sites.append(site)
			current_nucleotide = t[site]
			t[site] = mutate(current_nucleotide, mu, kappa, phi, True)
			# mutates the nucleotide based on the appropriate probabilties
		t = ''.join(t)
		daughter = t # replaces the old child with the new one

		# totals = detect_o_and_c(ancestor, strains[0], strains[1], mutation_sites1 = mutation_sites[0], mutation_sites2 = mutation_sites[1])

		cms[g] = len(mutation_sites) # totals['c']
	# print(cms)
	return cms