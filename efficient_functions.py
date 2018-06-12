import math
import numpy as np
from scipy import stats
from scipy import special

# time complexity: O(1)
def multiple(L, o, m1, m2):
	return ((special.comb(L,o,exact=False,repetition=False)) * (special.comb((L-o),(m2-o),exact=False,repetition=False)) * (special.comb((L-m2),(m1-o),exact=False,repetition=False)))

# time complexity: 0(1)
def pi_bar_formula(c, o, kappa, phi):
	prob = (kappa**2 + 1 - 2*phi + 2*(phi)**2)/((kappa+1)**2)
	return stats.binom.pmf(c,o,prob)

# time complexity: O(1)
def overlapping_prob(o, m1, m2, L):
	if(m1 < m2):
		temp = m1
		m1 = m2
		m2 = temp
	combos = (special.comb(L,m1,exact=False,repetition=False) * special.comb(L,m2,exact=False,repetition=False))
	return (multiple(L, o, m1, m2) * 1/(combos))

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

def better_prob_cm(c, mu, L, kappa, phi):
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
				innersum += (x * y * z)
		outersum = (innersum * pi_bar_formula(c, o, kappa, phi))
		prob.append(outersum)
		innersum = 0
	return sum(prob)

# time complexity: O(n^4), where n is L
def expected_cms(mu, L):
	value = 0
	# total = 0
	for c in range(L+1):
		value += (c * prob_cm(c, mu, L))
		# total += prob_cm(c, mu, L)
	return value

# time complexity: o(n^4), where n is L
def better_expected_cms(mu, L, kappa, phi):
	value = 0
	# total = 0
	for c in range(L+1):
		value += (c * better_prob_cm(c, mu, L, kappa, phi))
		# total += better_prob_cm(c, mu, L, kappa)
	# print(total)
	return value