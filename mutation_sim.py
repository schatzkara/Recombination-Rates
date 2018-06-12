import random
import numpy as np

# this sim starts with all As
# assumes equal probabilities of C,T,G
# forces a mutation
# and does not allow for multiple mutations in one site
# time comlexity: O(n^3), where n is the bigger of n and L; technically it's n^3 + n^2 + n
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

# this sim starts 2 identical but random strands
# factors in the probabilites of A,T,C,G
# allows 'mutation' to itself******************
# does not allow for multiple mutations in one spot
# time complexity: O(n^3), where n is the bigger of n and L; technically it's n^3 + n^2 + n
def sim_unequal_mutations(n, L, mu, kappa, phi):
	alpha = (mu * kappa)/(kappa + 1) # probability of transitions
	beta = mu/(kappa + 1) # probability of transversions
	totals = {}
	convergent_mutations = 0
	strains = []
	new = []
	nucleotides = ['A', 'T', 'G', 'C']
	weights_A = [(1-mu), (beta*phi), alpha, (beta*(1-phi))]
	weights_T = [(beta*phi), (1-mu), (beta*(1-phi)), alpha]
	weights_G = [alpha, (beta*(1-phi)), (1-mu), (beta*phi)]
	weights_C = [(beta*(1-phi)), alpha, (beta*phi), (1-mu)]

	# creates the initial strains
	# time comlexity: O(n), where n is the bigger of L and n; technically it's L+n
	ancestor = ''
	for y in range(L):
		ancestor+=random.choice(nucleotides)
	for x in range(n):
		strains.append(ancestor)

	# mutates the strains
	# time complexity: O(n^2), where n is the bigger of n and mu*L; technically it's n*mu*L
	for child in strains:
		t = list(child)
		for s in range(len(t)):
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
	o = 0
	c = 0
	for a in range(len(new)):
		strain1 = list(new[a])
		for b in range((a+1),len(new)):
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
	convergent_mutations = c
	# print(ancestor)
	# print(new)
	return convergent_mutations
