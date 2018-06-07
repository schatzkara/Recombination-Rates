import random

nucleotides = ['A', 'T', 'C', 'G']

def sim_mutations(n, L, mu):
	totals = {}
	strains = []
	new = []
	nucleotides = ['A', 'T', 'C', 'G']

	# creates the initial strains
	for x in range(n):
		s = ''
		for y in range(L):
			s+='A'
			# s+=random.choice(nucleotides)
		strains.append(s)
	# print(strains)

	# mutates the strains
	mutations = int(mu * L)
	for s in strains:
		t = list(s)
		for m in range(mutations):
			f = random.randint(1,(L-1))
			while(t[f] == 'A'):
				t[f] = random.choice(nucleotides)
		t = ''.join(t)
		new.append(t)
	# print(new)

	# counts up actual o and c values
	o = 0
	c = 0
	for a in range(len(new)):
		strain1 = list(new[a])
		for b in range((a+1),len(new)):
			strain2 = list(new[b])
			for d in range(len(strain1)):
				# print(strain1[d] + strain2[d])
				# counts up the number of overlapping mutation sites
				if (strain1[d] != 'A' and strain2[d] != 'A'):
					o += 1
				# counts up the number of overlapping and matching mutation sites
				if (strain1[d] != 'A' and strain1[d] == strain2[d]):
					c += 1
	totals['overlaps'] = o
	totals['matches'] = c

	return totals




# print(sim_mutations(2,10,.50))