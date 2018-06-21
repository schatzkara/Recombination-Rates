#! python3

# script for our model

import numpy as np
from scipy import stats
from scipy import special

# function to calculate the expected number of convergent mutations for a given mutation rate and DNA sequence length
# accounts for the unequal probabilities of switching to each other nucleotide
# params:
# 	mu (float) = mutation rate (units: mutations per base pair per generation)
# 	L (int) = length of DNA strand
# 	kappa (float) = proportion of transistions to transversions
# 	phi (float) = probability of a transition to one's complementary base
# return: float that equals the expected number of convergent mutations between 2 DNA strands of length L with mutation rate mu
# time complexity: O(n^3)
def expected_cms(L,mu,kappa,phi):
	sum1 = 0 # counter for sum of all o and c combinations
	sum2 = 0 # counter for sum of all o, c, and m2 combinations
	total = 0 # counter for total sum

	m_probs = prob_m(L,mu) # ordered list of the Poisson probabilties of each number of mutations with length L

	c_probs = prob_c(L,kappa,phi) # ordered list of the expected values of c for each possible value of o

	mutation_combos = combos(L) # ordered list of all the possible 'L choose m' values

	for m1 in range(L+1): # allows for all possible values of m1
		x = m_probs[m1]
		for m2 in range(L+1): # allows for all possible values of m2
			y = x * m_probs[m2]
			for o in range(min(m1,m2)+1): # allows for all possible values of o (note that o cannot be greater m1 OR m2 because then there can be no overlaps)
				# print('prob overlapping: ' + str(prob_overlapping(L,o,m1,m2,mutation_combos)) + ' c_prob: ' + str(c_probs[o]))
				z = prob_overlapping(L,o,m1,m2,mutation_combos) * c_probs[o]
				sum2 += z
			sum1 += y * sum2
			sum2 = 0
		total += sum1
		sum1 = 0

	return total

# function to calculate P(o;m1,m2,L): the probability of o overlapping sites given m1 and m2 mutations on strain 1 and strain 2 respectively
# params:
# 	o (int) = overlapping sites
# 	m1 (int) = mutations in strain 1
# 	m2 (int) = mutations in strain 2
# 	L (int) = length of DNA strand
# 	mutation_combos (list of ints) = ordered list of values of 'L choose m' for m = 0 to L
# return: float that equals the probability of o given m1 and m2
# time complexity: O(1)
def prob_overlapping(L,o,m1,m2,mutation_combos):
	prob = 0
	# mutation_combos = (L+1)*[None] # will be populated as an ordered list of 'L choose m' for m = 0 to L
	# for m in range(L+1): 
	# 	mutation_combos[m] = special.comb(L,m,exact=False,repetition=False)
	if(L-m1-m2-o < 1):
		prob = 0
	else:
		num = np.arange(L,L-m1-m2+o,-1, dtype = np.float) # array with each integer that would be multiplied to calculate L!/(L-m1-m2+o)!
		denom = np.concatenate((np.arange(o,0,-1, dtype = np.float),np.arange(m1-o,0,-1, dtype = np.float),np.arange(m2-o,0,-1, dtype = np.float))) # array with each integer that would be multiplied to calculate o!, m1!, and m2!
		# pairs = num/denom # elementalwise division
		# pairs = (np.arange(L,L-m1-m2+o,-1, dtype = np.float))/(np.concatenate((np.arange(o,0,-1, dtype = np.float),np.arange(m1-o,0,-1, dtype = np.float),np.arange(m2-o,0,-1, dtype = np.float))))
		# product = np.prod(pairs) # value of 'L choose o, m1, m2, L-m1-m2+o' (number of successful ways to get o overlapping sites)
		# combos = mutation_combos[m1] * mutation_combos[m2]

		# works all the computation in log form and then exponentiates at the end; keeps the magnitude in check
		combos = np.log(mutation_combos[m1]) + np.log(mutation_combos[m2]) # total number of ways to arrange m1 and m2 mutations on strain 1 and strain 2 respectively
		num = np.log(num)
		denom = np.log(denom)
		pairs = num-denom
		product = np.sum(pairs)
		prob = np.exp(product-combos)
	return prob

# function to calculate the Poisson probabilities of all possible values of m with expected value = mu*L
# params:
# 	L (int) = length of DNA strand
# 	mu (int) = mutation rate (in mutations per base pair per generation)
# return: ordered list of the Poisson probabilities for m = 0 to L (the ith element in the list corresponds to the probability of m=i)
# time complexity: O(n)
def prob_m(L,mu):
	m_probs = (L+1)*[None] # will be populated as an ordered list of the Poisson probabilities for m = 0 to L
	
	for m in range(L+1): 
		m_probs[m] = stats.poisson.pmf(m,(mu*L))

	return m_probs

# function to calculate the expected number of convergent mutations for all possible values of o
# params:
# 	L (int) = length of DNA strand
# 	kappa (float) = proportion of transistions to transversions
# 	phi (float) = probability of a transition to one's complementary base
# return: ordered list of the expected values of c for o = 0 to L (the ith element of the list corresponds to the expected value with o=i)
# time complexity: O(n^2)
def prob_c(L,kappa,phi):
	c_probs = (L+1)*[None] # will be populated as an ordered list of the expected values of c for o = 0 to L

	summation = 0 # counter for the total expected value

	for o in range(L+1): # allows for all possible values of o
		for c in range(o+1): # allows for all possible values of c
			summation += c * pi_bar(c,o,kappa,phi) # calculates the particular contribution to the expected value
		c_probs[o] = summation
		summation = 0

	return c_probs

# function to calculate the value of 'L choose m' for all possible values of m
# params:
# 	L (int) = length of DNA strand
# return: ordered list of the values of 'L choose m' for m = 0 to L (the ith element of the list corresponds to 'L choose i')
# time complexity: O(n)
def combos(L):
	mutation_combos = (L+1)*[None] # will be populated as an ordered list of 'L choose m' for m = 0 to L

	for m in range(L+1): # allows for all possible values of m
		mutation_combos[m] = special.comb(L,m,exact=False,repetition=False)

	return mutation_combos

# function to calculate pi_bar(c;o): the probability of c convergent mutations given o overlapping sites
# accounts for different probabilties of transitions versus transversions
# params: 
# 	c (int) = convergent sites
# 	o (int) = overlapping sites
# 	kappa (float) = ratio of transitions to transversions
# 	phi (float) = probability of transversion to its complementary base pair
# return: float that equals the probability of c given o
# time complexity: 0(1)
def pi_bar(c,o,kappa,phi):
	prob = (kappa**2 + 1 - 2*phi + 2*(phi)**2)/((kappa+1)**2) # the probability of some convergent mutation aka a 'success' in the binomial probability
	return stats.binom.pmf(c,o,prob) # pi_bar is equivalent to the binomial probability density function

def mutation_matrix(mu, kappa, phi, generations):
	alpha = (mu * kappa)/(kappa + 1) # probability of transitions
	beta = mu/(kappa + 1) # probability of transversions
	m = np.matrix([[(1-mu), beta*phi, alpha, beta*(1-phi)], [beta*phi, (1-mu), beta*(1-phi), alpha], [alpha, beta*(1-phi), (1-mu), beta*phi], [beta*(1-phi), alpha, beta*phi, (1-mu)]])
	mg = np.linalg.matrix_power(m, generations)
	# for g in range(generations-1):
	# 	m = np.dot(m,m)
	return mg 

def expected_idp(mu, kappa, phi, generations):
	print('entered function')
	print(generations)
	mg = mutation_matrix(mu, kappa, phi, generations)
	print('got m^g')
	expected_idp = 0
	for x in range(4):
		expected_idp += ((mg.item((0,x)))**2)
		print(x)
		print(expected_idp)
	# x.item((0, 1))

def expected_cms_given_m(L,mutations,kappa,phi):
	# sum1 = 0 # counter for sum of all o and c combinations
	# sum2 = 0 # counter for sum of all o, c, and m2 combinations
	total = 0 # counter for total sum

	# m_probs = prob_m(L,mu) # ordered list of the Poisson probabilties of each number of mutations with length L

	c_probs = prob_c(L,kappa,phi) # ordered list of the expected values of c for each possible value of o

	mutation_combos = combos(L) # ordered list of all the possible 'L choose m' values

	# for m1 in range(L+1): # allows for all possible values of m1
	# 	x = m_probs[m1]
	# 	for m2 in range(L+1): # allows for all possible values of m2
	# 		y = x * m_probs[m2]
	for o in range(mutations+1): # allows for all possible values of o (note that o cannot be greater m1 OR m2 because then there can be no overlaps)
		# print('prob overlapping: ' + str(prob_overlapping(L,o,m1,m2,mutation_combos)) + ' c_prob: ' + str(c_probs[o]))
		total += prob_overlapping(L,o,mutations,mutations,mutation_combos) * c_probs[o]
		# sum2 += z
	# sum1 += y * sum2
	# sum2 = 0
	# total += sum1
	# sum1 = 0

	return total

def expected_cms_with_mg(L, mu, kappa, phi):
	sum1 = 0 # counter for sum of all o and c combinations
	sum2 = 0 # counter for sum of all o, c, and m2 combinations
	total = 0 # counter for total sum

	m_probs = prob_m(L,mu) # ordered list of the Poisson probabilties of each number of mutations with length L

	c_probs = prob_c_with_mg(L,kappa,phi) # ordered list of the expected values of c for each possible value of o

	mutation_combos = combos(L) # ordered list of all the possible 'L choose m' values

	for m1 in range(L+1): # allows for all possible values of m1
		x = m_probs[m1]
		for m2 in range(L+1): # allows for all possible values of m2
			y = x * m_probs[m2]
			for o in range(min(m1,m2)+1): # allows for all possible values of o (note that o cannot be greater m1 OR m2 because then there can be no overlaps)
				# print('prob overlapping: ' + str(prob_overlapping(L,o,m1,m2,mutation_combos)) + ' c_prob: ' + str(c_probs[o]))
				z = prob_overlapping(L,o,m1,m2,mutation_combos) * c_probs[o]
				sum2 += z
			sum1 += y * sum2
			sum2 = 0
		total += sum1
		sum1 = 0

	return total

def prob_c_with_mg(L,kappa,phi):
	c_probs = (L+1)*[None] # will be populated as an ordered list of the expected values of c for o = 0 to L

	summation = 0 # counter for the total expected value

	for o in range(L+1): # allows for all possible values of o
		for c in range(o+1): # allows for all possible values of c
			summation += c * pi_bar(c,o,kappa,phi) # calculates the particular contribution to the expected value
		c_probs[o] = summation
		summation = 0

	return c_probs

def pi_bar_with_mg(c,o,kappa,phi):
	prob = (kappa**2 + 1 - 2*phi + 2*(phi)**2)/((kappa+1)**2) # the probability of some convergent mutation aka a 'success' in the binomial probability
	return stats.binom.pmf(c,o,prob) # pi_bar is equivalent to the binomial probability density function