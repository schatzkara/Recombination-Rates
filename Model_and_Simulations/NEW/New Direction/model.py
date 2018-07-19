# # script for our model

import numpy as np
from scipy import stats
from scipy import special
# from simulations import c_q_sim

# # function to calculate the expected number of convergent mutations for a given mutation rate and DNA sequence length
# # accounts for the unequal probabilities of switching to each other nucleotide
# # params:
# # 	mu_star (float) = effective mutation rate (units: mutations per base pair)
# # 	L (int) = length of DNA strand
# # 	kappa (float) = proportion of transistions to transversions
# # 	phi (float) = probability of a transition to one's complementary base
# # return: float that equals the expected number of convergent mutations between 2 DNA strands of length L with mutation rate mu
# # time complexity: O(n^3)
# def expected_cms(L,mu_star,kappa,phi):
# 	sum1 = 0 # counter for sum of all o and c combinations
# 	sum2 = 0 # counter for sum of all o, c, and m2 combinations
# 	total = 0 # counter for total sum
	
# 	cutoff = int(stats.poisson.ppf(.99, mu_star*L))+1 # this is the point at which m1 and m2 become negligible with 99.99% confidence

# 	m_probs = prob_m(L,mu,cutoff) # ordered list of the Poisson probabilties of each number of mutations with length L

# 	c_probs = prob_c(L,kappa,phi,cutoff) # ordered list of the expected values of c for each possible value of o

# 	mutation_combos = combos(L,cutoff) # ordered list of all the possible 'L choose m' values

# 	for m1 in range(cutoff+1): # (L+1): # allows for all possible values of m1
# 		x = m_probs[m1]
# 		for m2 in range(cutoff+1): # (L+1): # allows for all possible values of m2
# 			y = x * m_probs[m2]
# 			for o in range(min(m1,m2)+1): # allows for all possible values of o (note that o cannot be greater m1 OR m2 because then there can be no overlaps)
# 				z = prob_overlapping(L,o,m1,m2,mutation_combos) * c_probs[o]
# 				sum2 += z
# 			sum1 += y * sum2
# 			sum2 = 0
# 		total += sum1
# 		sum1 = 0
# 	# print('expected cms')
# 	return total

# function to calculate P(o;m1,m2,L): the probability of o overlapping sites given m1 and m2 mutations on strain 1 and strain 2 respectively
# params:
# 	o (int) = overlapping sites
# 	m1 (int) = mutations in strain 1
# 	m2 (int) = mutations in strain 2
# 	L (int) = length of DNA strand
# 	mutation_combos (list of ints) = ordered list of values of 'L choose m' for m = 0 to L
# return: float that equals the probability of o given m1 and m2
# time complexity: O(1)
def prob_overlapping(L,o,m1,m2,mutation_combos_m1, mutation_combos_m2):
	prob = 0
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
		combos = np.log(mutation_combos_m1) + np.log(mutation_combos_m2) # total number of ways to arrange m1 and m2 mutations on strain 1 and strain 2 respectively
		num = np.log(num)
		denom = np.log(denom)
		pairs = num-denom
		product = np.sum(pairs)
		prob = np.exp(product-combos)
	# print('prob overlapping done')
	return prob

# # function to calculate the Poisson probabilities of all possible values of m with expected value = mu*L
# # params:
# # 	L (int) = length of DNA strand
# # 	mu (int) = mutation rate (in mutations per base pair per generation)
# # return: ordered list of the Poisson probabilities for m = 0 to L (the ith element in the list corresponds to the probability of m=i)
# # time complexity: O(n)
# def prob_m(L,mu,cutoff):
# 	m_probs = (cutoff+1)*[None] # (L+1)*[None] # will be populated as an ordered list of the Poisson probabilities for m = 0 to L
	
# 	for m in range(cutoff+1): # (L+1): 
# 		m_probs[m] = stats.poisson.pmf(m,(mu*L))

# 	# print('prob m done')
# 	return m_probs

# # function to calculate the expected number of convergent mutations for all possible values of o
# # params:
# # 	L (int) = length of DNA strand
# # 	kappa (float) = proportion of transistions to transversions
# # 	phi (float) = probability of a transition to one's complementary base
# # return: ordered list of the expected values of c for o = 0 to L (the ith element of the list corresponds to the expected value with o=i)
# # time complexity: O(n^2)
# def prob_c(L,kappa,phi,cutoff):
# 	c_probs = (cutoff+1)*[None] # (L+1)*[None] # will be populated as an ordered list of the expected values of c for o = 0 to L

# 	summation = 0 # counter for the total expected value

# 	for o in range(cutoff+1): # (L+1): # allows for all possible values of o
# 		for c in range(o+1): # allows for all possible values of c
# 			summation += c * pi_bar(c,o,kappa,phi) # calculates the particular contribution to the expected value
# 		c_probs[o] = summation
# 		summation = 0

# 	# print('prob c done')
# 	return c_probs

# function to calculate the value of 'L choose m' for all possible values of m
# params:
# 	L (int) = length of DNA strand
# return: ordered list of the values of 'L choose m' for m = 0 to L (the ith element of the list corresponds to 'L choose i')
# time complexity: O(n)
def combos(L,cutoff):
	mutation_combos = (cutoff+1)*[None] # (L+1)*[None] # will be populated as an ordered list of 'L choose m' for m = 0 to cutoff value

	for m in range(cutoff+1): # (L+1): # allows for all possible values of m
		mutation_combos[m] = special.comb(L,m,exact=False,repetition=False)

	# print('mutation combos done')
	return mutation_combos

# # function to calculate pi_bar(c;o): the probability of c convergent mutations given o overlapping sites
# # accounts for different probabilties of transitions versus transversions
# # params: 
# # 	c (int) = convergent sites
# # 	o (int) = overlapping sites
# # 	kappa (float) = ratio of transitions to transversions
# # 	phi (float) = probability of transversion to its complementary base pair
# # return: float that equals the probability of c given o
# # time complexity: 0(1)
# def pi_bar(c,o,kappa,phi):
# 	prob = (kappa**2 + 1 - 2*phi + 2*(phi)**2)/((kappa+1)**2) # the probability of some convergent mutation aka a 'success' in the binomial probability
# 	# print('pi bar done')
# 	return stats.binom.pmf(c,o,prob) # pi_bar is equivalent to the binomial probability density function

	











# # script for our model using the mutation matrix to the power of the generation (M^g)

# def expected_cms_with_mg(L, mu, kappa, phi, generations):
# 	sum1 = 0 # counter for sum of all o and c combinations
# 	sum2 = 0 # counter for sum of all o, c, and m2 combinations
# 	total = 0 # counter for total sum

# 	cutoff = int(stats.poisson.ppf(.99, mu*L*generations))+1 # this is the point at which m1 and m2 become negligible with 99.99% confidence
# 	print('cutoff: ' + str(cutoff))
# 	m_probs = prob_m(L,mu,cutoff) # ordered list of the Poisson probabilties of each number of mutations with length L
# 	print('m probs done')
# 	mg = mutation_matrix(mu, kappa, phi, generations)
# 	print('m^g done')
# 	c_probs = prob_c_with_mg(L,kappa,phi,mg,cutoff) # ordered list of the expected values of c for each possible value of o
# 	print('c probs done')
# 	mutation_combos = combos(L, cutoff) # ordered list of all the possible 'L choose m' values
# 	print('mutation combos done')
# 	for m1 in range(cutoff+1): # allows for all possible values of m1
# 		x = m_probs[m1]
# 		for m2 in range(cutoff+1): # allows for all possible values of m2
# 			y = x * m_probs[m2]
# 			for o in range(min(m1,m2)+1): # allows for all possible values of o (note that o cannot be greater m1 OR m2 because then there can be no overlaps)
# 				z = prob_overlapping(L,o,m1,m2,mutation_combos) * c_probs[o]
# 				print('prob overlapping done')
# 				sum2 += z
# 			sum1 += y * sum2
# 			sum2 = 0
# 		print('getting there')
# 		total += sum1
# 		sum1 = 0

# 	return total

def prob_c_with_mg(L,kappa,phi,mg_1,mg_2,cutoff):
	c_probs = (cutoff+1)*[None] # will be populated as an ordered list of the expected values of c for o = 0 to cutoff

	summation = 0 # counter for the total expected value

	for o in range(1,cutoff+1): # allows for all possible values of o
		for c in range(1,o+1): # allows for all possible values of c
			summation += c * pi_bar_with_mg(c,o,kappa,phi,mg_1,mg_2) # calculates the particular contribution to the expected value
			print(c)
			print(pi_bar_with_mg(c,o,kappa,phi,mg_1,mg_2))
		c_probs[o] = summation
		summation = 0

	return c_probs

def pi_bar_with_mg(c,o,kappa,phi,mg_1,mg_2):
	# prob = (mg[0,1]**2 + mg[0,2]**2 + mg[0,3]**2)/((mg[0,1] + mg[0,2] + mg[0,3])**2)
	prob = (mg_1[0,1]*mg_2[0,1] + mg_1[0,2]*mg_2[0,2] + mg_1[0,3]*mg_2[0,3])/((mg_1[0,1] + mg_1[0,2] + mg_1[0,3])*(mg_2[0,1] + mg_2[0,2] + mg_2[0,3]))
	num = (mg_1[0,1]*mg_2[0,1] + mg_1[0,2]*mg_2[0,2] + mg_1[0,3]*mg_2[0,3])
	denom = ((mg_1[0,1] + mg_1[0,2] + mg_1[0,3])*(mg_2[0,1] + mg_2[0,2] + mg_2[0,3]))
	# print(num)
	# print(denom)
	return stats.binom.pmf(c,o,prob) # pi_bar is equivalent to the binomial probability density function

def mutation_matrix(mu, kappa, phi, generations):
	alpha = (mu * kappa)/(kappa + 1) # probability of transitions
	beta = mu/(kappa + 1) # probability of transversions
	m = np.matrix([[(1-mu), beta*phi, alpha, beta*(1-phi)], [beta*phi, (1-mu), beta*(1-phi), alpha], [alpha, beta*(1-phi), (1-mu), beta*phi], [beta*(1-phi), alpha, beta*phi, (1-mu)]], dtype = np.float)
	mg = np.linalg.matrix_power(m, generations)
	return mg

# print(mutation_matrix(.01, 2.5, 0.5, 0))
	
# # model using uvwxyz
	
# def mutation_matrix_uvwxyz(mu, u,v,w,x,y,z, generations):
# 	m = np.matrix([[(1-mu), w, u, x], [w, (1-mu), y, v], [u, y, (1-mu), z], [x, v, z, (1-mu)]], dtype = np.float)
# 	mg = np.linalg.matrix_power(m, generations)
# 	return mg
	
# def prob_c_with_mg_uvwxyz(L,kappa,phi,mg,cutoff):
# 	c_probs = (cutoff+1)*[None] # will be populated as an ordered list of the expected values of c for o = 0 to L

# 	summation = 0 # counter for the total expected value

# 	for o in range(cutoff+1): # allows for all possible values of o
# 		for c in range(o+1): # allows for all possible values of c
# 			summation += c * pi_bar_with_mg_uvwxyz(c,o,kappa,phi,mg) # calculates the particular contribution to the expected value
# 		c_probs[o] = summation
# 		summation = 0

# 	return c_probs
	
# def pi_bar_with_mg_uvwxyz(c,o,kappa,phi,mg):
# 	prob = (mg[0,1]**2 + mg[0,2]**2 + mg[0,3]**2)/((mg[0,1] + mg[0,2] + mg[0,3])**2)
# 	return stats.binom.pmf(c,o,prob) # pi_bar is equivalent to the binomial probability density function
	
# def expected_cms_with_mg_uvwxyz(L, mu, kappa, phi, generations):
# 	sum1 = 0 # counter for sum of all o and c combinations
# 	sum2 = 0 # counter for sum of all o, c, and m2 combinations
# 	total = 0 # counter for total sum

# 	cutoff = int(stats.poisson.ppf(.99, mu*L*generations))+1 # this is the point at which m1 and m2 become negligible with 99.99% confidence

# 	m_probs = prob_m(L,mu,cutoff) # ordered list of the Poisson probabilties of each number of mutations with length L

# 	mg = mutation_matrix(mu, kappa, phi, generations)

# 	c_probs = prob_c_with_mg_uvwxyz(L,kappa,phi,mg,cutoff) # ordered list of the expected values of c for each possible value of o

# 	mutation_combos = combos(L, cutoff) # ordered list of all the possible 'L choose m' values

# 	for m1 in range(cutoff+1): # allows for all possible values of m1
# 		x = m_probs[m1]
# 		for m2 in range(cutoff+1): # allows for all possible values of m2
# 			y = x * m_probs[m2]
# 			for o in range(min(m1,m2)+1): # allows for all possible values of o (note that o cannot be greater m1 OR m2 because then there can be no overlaps)
# 				z = prob_overlapping(L,o,m1,m2,mutation_combos) * c_probs[o]
# 				sum2 += z
# 			sum1 += y * sum2
# 			sum2 = 0
# 		total += sum1
# 		sum1 = 0

# 	return total












# def expected_cms_given_m(L,mu,generations,kappa,phi):
# 	total = 0 # counter for total sum
# 	mutations = int(mu*L*generations)

# 	# cutoff = int(stats.poisson.ppf(.9999, mu*L))+1 # this is the point at which m1 and m2 become negligible with 99.99% confidence

# 	mg = mutation_matrix(mu, kappa, phi, generations)

# 	c_probs = prob_c_with_mg(L,kappa,phi,mg,mutations) # ordered list of the expected values of c for each possible value of o

# 	mutation_combos = combos(L,mutations) # ordered list of all the possible 'L choose m' values

# 	# mutation_sites = int((mu*L)/expected_m_at_site(mu, generations))

# 	for o in range(mutations+1): # allows for all possible values of o (note that o cannot be greater than m1 OR m2 because then there can be no overlaps)
# 		total += prob_overlapping(L,o,mutations,mutations,mutation_combos) * c_probs[o]

# 	return total


# def expected_m_at_site(mu, g):
# 	return (g*(mu**(g+1)))/(1-mu) + (mu - mu**(g+1))/((1-mu)**2)

# def expected_mutation_sites(L, mi):
# 	total = 0
# 	for k in range(1,mi):
# 		num = special.comb((mi-1),(k-1),exact=False,repetition=False) * special.comb(L,k)
# 		denom = L**mi
# 		total += (num/denom) * k
# 	return total


def expected_c_given_ms(L, m1, m2, mu, generations_1, generations_2, kappa, phi):
	print('Calculating the expected number of convergent mutations.\n')
	total = 0 # counter for total sum
	# mutations = int(mu*L*generations)

	# cutoff = int(stats.poisson.ppf(.9999, mu*L))+1 # this is the point at which m1 and m2 become negligible with 99.99% confidence
	# print('\tGenerating the mutation matrices.')
	mg_1 = mutation_matrix(mu, kappa, phi, generations_1)
	mg_2 = mutation_matrix(mu, kappa, phi, generations_2)
	# print(mg_1)
	# print(mg_2)

	# print('\tGetting the c probs.')
	# c_probs = prob_c_with_mg(L, kappa, phi, mg_1, mg_2, min(m1,m2)) # ordered list of the expected values of c for each possible value of o
	# c_probs_2 = prob_c_with_mg(L, kappa, phi, mg_2, m2)

	o_probs = [] # index = o-1

	mutation_combos_1 = special.comb(L,m1,exact=False,repetition=False)
	mutation_combos_2 = special.comb(L,m2,exact=False,repetition=False)
	# print(mutation_combos_1)
	# print(mutation_combos_2)
	# mutation_combos = combos(L, max(m1,m2)) # ordered list of all the possible 'L choose m' values

	# mutation_sites = int((mu*L)/expected_m_at_site(mu, generations))
	# print('\tGetting o probs.')
	for o in range(1, min(m1,m2)+1): # allows for all possible values of o (note that o cannot be greater than m1 OR m2 because then there can be no overlaps)
		prob_o = prob_overlapping(L, o, m1, m2, mutation_combos_1, mutation_combos_2)
		o_probs.append(prob_o)
		if prob_o < .0000001:
			break
	max_o = len(o_probs) - 1

	# print(o_probs)

	# print('\tGetting c probs.')
	c_probs = prob_c_with_mg(L, kappa, phi, mg_1, mg_2, max_o) # ordered list of the expected values of c for each possible value of o
	# print(c_probs)

	# print('\tGetting expexcted value')
	for o in range(1,max_o):
		total += o_probs[o-1] * c_probs[o]

	# print(total)

	return total





# def expected_idp(mu, kappa, phi, generations):
# 	# print('entered function')
# 	# print(generations)
# 	mg = mutation_matrix(mu, kappa, phi, generations)
# 	# print('got m^g')
# 	expected_idp = 0
# 	for x in range(4):
# 		expected_idp += (mg[0,x])**2
# 	return expected_idp

# def id_equation(L,m1,m2,o,c):
#         return (L-m1-m2+o+c)/L






# def calc_correction_factor(n, c_q_list):
# 	correction_factor = 0
# 	for q in range(3,n+1):
# 		c_q = c_q_list[q-2]
# 		overcounts = ((q*(q-1))/2) - 1
# 		correction_factor += c_q * overcounts
# 	return correction_factor

# def c_q_list(n, L, mu, generations, kappa, phi):
# 	c_q_list = (n-1)*[None]
# 	all_c_q = (n-1)*[None]
# 	for x in range(n-1):
# 		all_c_q[x] = 0 # 1000*[None]
# 	for i in range(1000):
# 		c_qs = c_q_sim(n, L, mu, generations, kappa, phi)
# 		for q in range(2,n+1):
# 			all_c_q[q-2] += c_qs[q-2] # all_c_q[q-2][i] = c_qs[q-2]
# 	for y in range(len(all_c_q)):
# 		c_q_list[y] = all_c_q[y]/1000 # np.mean(all_c_q[y])

# def expected_h_c_mu_star(n, L, mu, generations, kappa, phi):
# 	mu_star = mu*generations
# 	num_pairs = ((n*(n-1))/2)
# 	# correction_factor = calc_correction_factor(n, c_q_list(n, L, mu, generations, kappa, phi))
# 	return expected_cms(L, mu_star, kappa, phi) * num_pairs


# def expected_h_c_mg(n, L, mu, generations, kappa, phi):
# 	num_pairs = ((n*(n-1))/2)
# 	# correction_factor = calc_correction_factor(n, c_q_list(n, L, mu, generations, kappa, phi))
# 	return expected_cms_with_mg(L, mu, kappa, phi, generations) * num_pairs

# for g in range(5000,5100):
# 	mg = mutation_matrix(.00004, 2.5, 0.5, g)
# 	print(mg)
# 	pi_bar_with_mg(2,4,2.5,0.5,mg,mg)
