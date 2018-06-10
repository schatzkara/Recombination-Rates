import random
import numpy as np

# this sim starts with all As
# assumes equal probabilities of C,T,G
# forces a mutation
# and does not allow for multiple mutations in one site
def sim_equal_mutations(n, L, mu):
	totals = {}
	convergent_mutations = 0
	strains = []
	new = []
	nucleotides = ['A', 'T', 'C', 'G']

	# creates the initial strains
	for x in range(n):
		s = ''
		for y in range(L):
			s+='A'
		strains.append(s)

	# mutates the strains
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
def sim_unequal_mutations(n, L, mu):
	totals = {}
	convergent_mutations = 0
	strains = []
	new = []
	nucleotides = ['A', 'T', 'C', 'G']
	weights_A = [0.2, 0.5, 0.2, 0.1]
	weights_T = [0.2, 0.5, 0.2, 0.1]
	weights_G = [0.2, 0.5, 0.2, 0.1]
	weights_C = [0.2, 0.5, 0.2, 0.1]

	# creates the initial strains
	ancestor = ''
	for y in range(L):
		ancestor+=random.choice(nucleotides)
	for x in range(n):
		strains.append(ancestor)

	# mutates the strains
	mutations = int(mu * L)
	for child in strains:
		t = list(child)
		f = random.randint((L-1), mutations, replace=False)
		for m in f:
			current = t[m]
			if current == 'A':
				t[m] = random.choice(nucleotides, p=weights_A)
			elif current == 'T':
				t[m] = random.choice(nucleotides, p=weights_T)
			elif current == 'C':
				t[m] = random.choice(nucleotides, p=weights_C)
			elif current == 'G':
				t[m] = random.choice(nucleotides, p=weights_G)
		t = ''.join(t)
		new.append(t)

	# counts up actual o and c values
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

	return convergent_mutations
