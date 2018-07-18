# (O(L,g))/(L^g)
def prob_no_d(L, g):
	num = np.arange(L, L-g+1,-1, dtype = np.float)
	denom = np.full((1,g), L)
	pairs = num/denom
	prob = np.prod(pairs)
	return prob ** 2
	
L = [10]
gen = 300
	
for l in L:
	for g in range(gen):
		print(prob_no_d(l,g))