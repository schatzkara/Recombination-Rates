#! python3

# script for our model

import math
import numpy as np
from scipy import stats
from scipy import special

# function to calculate the multinomial coefficient for our purposes for the P(o;m1,m2,L) function: the probability of o overlapping sites given m1 and m2 mutations on strain 1 and strain 2 respectively
# params: 
# 	L (int) = length of DNA strand
# 	o (int) = overlapping sites
# 	m1 (int) = mutations in strain 1
# 	m2 (int) = mutations in strain 2
# return: int that equals the number of ways of 'successfully' arranging m1 and m2 along L to get o
# time complexity: O(1)
def multiple(L, o, m1, m2):
	# the multinomial coefficient calculated via multiplying 3 binomial coefficients
	return ((special.comb(L,o,exact=False,repetition=False)) * (special.comb((L-o),(m2-o),exact=False,repetition=False)) * (special.comb((L-m2),(m1-o),exact=False,repetition=False)))

# function to calculate pi_bar(c;o): the probability of c convergent mutations given o overlapping sites
# accounts for different probabilties of transitions versus transversions
# params: 
# 	c (int) = convergent sites
# 	o (int) = overlapping sites
# 	kappa (float) = ratio of transitions to transversions
# 	phi (float) = probability of transversion to its complementary base pair
# return: float that equals the probability of c given o
# time complexity: 0(1)
def pi_bar_formula(c, o, kappa, phi):
	prob = (kappa**2 + 1 - 2*phi + 2*(phi)**2)/((kappa+1)**2) # the probability of some convergent mutation aka a 'success' in the binomial probability
	return stats.binom.pmf(c,o,prob) # pi_bar is equivalent to the binomial probability density function

# function to calculate P(o;m1,m2,L): the probability of o overlapping sites given m1 and m2 mutations on strain 1 and strain 2 respectively
# params:
# 	o (int) = overlapping sites
# 	m1 (int) = mutations in strain 1
# 	m2 (int) = mutations in strain 2
# 	L (int) = length of DNA strand
# 	mutation_combos (list of ints) = ordered values of 'L choose m' for m = 0 to L
# return: float that equals the probability of o given m1 and m2
# time complexity: O(1)
def overlapping_prob(o, m1, m2, L, mutation_combos):
	if(m1 < m2): # forces m1 > m2 WLOG for computational ease
		temp = m1
		m1 = m2
		m2 = temp
	combos = mutation_combos[m1] * mutation_combos[m2]
	# combos = (special.comb(L,m1,exact=False,repetition=False) * special.comb(L,m2,exact=False,repetition=False)) # total number of ways to place m1 and m2 along L
	return (multiple(L, o, m1, m2) / combos) 

# function to calcluate P(c): the probability of c convergent mutations over all values of o,m1,m2
# assumes that the probability of switching to each other nucleotide is equivalent
# params: 
# 	c (int) = convergent sites
# 	mu (float) = mutation rate (units: mutations per base pair per generation)
# 	L (int) = length of DNA strand
# return: float that equals the total probability of c (over all values of o,m1,m2)
# time complexity: O(n^3), where n is L
def prob_cm(c, mu, L):
	expected = mu*L
	mu_probs = []
	for m in range(L+1):
		mu_probs.append(stats.poisson.pmf(m,expected))

	innersum = 0
	outersum = 0
	prob = []
	for o in range(L+1):
		for m1 in range(L+1):
			for m2 in range(L+1):
				x = mu_probs[m1]
				y = mu_probs[m2]
				z = overlapping_prob(o, m1, m2, L)
				innersum += (x * y *z)
		outersum = (innersum * stats.binom.pmf(c,o,1/3)) # pi formula
		prob.append(outersum)
		innersum = 0
	return sum(prob)

# function to calcluate modified P(c): the probability of c convergent mutations over all values of o,m1,m2
# accounts for the different probabilities of transitions versus transversions
# params:
# 	c (int) = convergent sites
# 	mu (float) = mutation rate (units: mutations per base pair per generation)
# 	L (int) = length of DNA strand
# 	kappa (float) = ratio of transitions to transversions
# 	phi (float) = probability of transversion to its complementary base pair
# 	mu_probs (list of floats) = ordered list of the Poisson probabilites for m = 0 to L
# 	mutation_combos (list of ints) = ordered values of 'L choose m' for m = 0 to L
# return: float that equals the total probability of c (over all values of o,m1,m2)
# time complexity: O(n^3), where n is L
def better_prob_cm(c, mu, L, kappa, phi, mu_probs, mutation_combos):
	innersum = 0 # counter for sum of the inner 2 summations
	outersum = 0 # counter for the total summation
	# prob = [] # list to hold the proabilities of c for each combination of m1,m2,L

	for o in range(L+1):
		for m1 in range(L+1):
			for m2 in range(L+1):
				x = mu_probs[m1]
				y = mu_probs[m2]
				z = overlapping_prob(o, m1, m2, L, mutation_combos)
				innersum += (x * y * z)
		outersum += (innersum * pi_bar_formula(c, o, kappa, phi))
		# prob.append(outersum)
		innersum = 0
	return outersum
	# return sum(prob)

# function to calculate the expected number of convergent mutations for a given mutation rate and DNA sequence length
# assumes that the probability of switching to each other nucleotide is equivalent
# params:
# 	mu (float) = mutation rate (units: mutations per base pair per generation)
# 	L (int) = length of DNA strand
# return: float that equals the expected number of convergent mutations between 2 DNA strands of length L with mutation rate mu
# time complexity: O(n^4), where n is L
def expected_cms(mu, L):
	value = 0
	# total = 0
	for c in range(L+1):
		value += (c * prob_cm(c, mu, L))
		# total += prob_cm(c, mu, L)
	return value

# function to calculate the expected number of convergent mutations for a given mutation rate and DNA sequence length
# accounts for the different probabilities of transitions versus transversions
# params:
# 	mu (float) = mutation rate (units: mutations per base pair per generation)
# 	L (int) = length of DNA strand
# 	kappa (float) = ratio of transitions to transversions
# 	phi (float) = probability of transversion to its complementary base pair
# return: float that equals the expected number of convergent mutations between 2 DNA strands of length L with mutation rate mu
# time complexity: o(n^4), where n is L
def better_expected_cms(mu, L, kappa, phi):
	value = 0
	# expected = mu*L
	mu_probs = (L+1)*[None] # will be populated as an ordered list of the Poisson probabilities for m = 0 to L
	mutation_combos = (L+1)*[None] # will be populated as an ordered list of 'L choose m' for m = 0 to L
	for m in range(L+1): 
		mu_probs[m] = stats.poisson.pmf(m,(mu*L))
		mutation_combos[m] = special.comb(L,m,exact=False,repetition=False)

	for c in range(L+1): # calculates the expected value
		value += (c * better_prob_cm(c, mu, L, kappa, phi, mu_probs, mutation_combos))
	return value