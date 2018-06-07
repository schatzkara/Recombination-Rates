import random

nucleotides = ['A', 'T', 'C', 'G']

def sim_overlaps(n, L, mu):
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
	# return strains
	mutations = int(mu * L)
	# mutates the strains
	for s in strains:
		t = list(s)
		for m in range(mutations):
			c = random.randint(1,(L-1))
			while(t[c] == 'A'):
				t[c] = random.choice(nucleotides)
		t = ''.join(t)
		new.append(t)
	print(new)

	# counts up the number of overlapping mutation sites
	o = 0
	for a in range(len(new)):
		strain1 = list(new[a])
		for b in range((a+1),len(new)):
			strain2 = list(new[b])
			for d in range(len(strain1)):
				print(strain1[d] + strain2[d])
				if (strain1[d] != 'A' and strain2[d] != 'A'):
					o += 1
	return o

def sim_matches(n, L, mu):
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
	# return strains
	mutations = int(mu * L)
	# mutates the strains
	for s in strains:
		t = list(s)
		for m in range(mutations):
			c = random.randint(1,(L-1))
			while(t[c] == 'A'):
				t[c] = random.choice(nucleotides)
		t = ''.join(t)
		new.append(t)
	print(new)

	# counts up the number of overlapping and matching mutation sites
	c = 0
	for a in range(len(new)):
		strain1 = list(new[a])
		for b in range((a+1),len(new)):
			strain2 = list(new[b])
			for d in range(len(strain1)):
				print(strain1[d] + strain2[d])
				if (strain1[d] != 'A' and strain1[d] == strain2[d]):
					# print(strain1[d])
					c += 1
	return c




print(sim_overlaps(2,10,.50))
print(sim_matches(2,10,.50))