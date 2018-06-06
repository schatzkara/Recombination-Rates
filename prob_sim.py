import random



mu = .20
L = 10
nucleotides = ['A', 'T', 'C', 'G']

def sim_overlaps(n, L, mu):
	strains = []
	new = []
	nucleotides = ['A', 'T', 'C', 'G']
	for x in range(n):
		s = ''
		for y in range(L):
			s+='A'
			# s+=random.choice(nucleotides)
		strains.append(s)
	# return strains
	mutations = int(mu * L)
	for s in strains:
		t = list(s)
		for m in range(mutations):
			c = random.randint(1,(L-1))
			while(t[c] == 'A'):
				t[c] = random.choice(nucleotides)
			t = ''.join(t)
		new.append(t)
	return new




print(sim_overlaps(10,10,.20))