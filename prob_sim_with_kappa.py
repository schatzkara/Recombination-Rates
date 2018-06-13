#! python3

# script to run the mutation simulation over multiple iterations and compare the actual number of convergent mutations to that from our model
# time complexity for each mu,L combination: O(n^4), where n is the biggest of L, n, iterations; technically it's 2n^4 + n^3 + 2n^2 + 2n ((L^4 + L) + iterations * (n^2 * L + n*L + n + L) + iterations)
# time complexity: O(n^6), where n is the biggest of len(L), len(mu), L, n, iterations; technically it's 2n^6 + n^5 + 2n^4 + 2n^3 (len(L) * len(mu) * [(L^4 + L) + iterations * (n^2 * L + n*L + n + L) + iterations])

from efficient_functions import better_expected_cms
from mutation_sim import sim_unequal_mutations
import csv


# params:
# 	n (int) = number of DNA strands
# 	L (list of ints) = list of DNA strand lengths (units: base pairs)
# 	mu (list of floats) = list of mutation rates (units: mutations per base pair per generation)
# 	kappa (float) = ratio of transitions to transversions
# 	phi (float) = probabiltiy of transversion to its complementary base pair
# 	iterations (int) = number of iterations to run
# return: a .csv file with the data from all iterations and the expected value

# enter parameters here
n = 2 # number of DNA strands
L = list(range(200, 1001, 100)) # format for range(): list(range(starting length, ending length+1, increment))
mu = [0.1]
kappa = 3 
phi = 1/2 
iterations = 1000

# actually runs the simulations
for l in L: # iterates over every length desired
	for m in mu: # iterates over every mu desired
		total_cms = 0 # counter for total convergent mutations over all iterations
		data = {} # dictionary for data over all iterations (key: iteration number, value: average convergent mutations of all iterations up to that point)
		better_expected = better_expected_cms(m, l, kappa, phi) # number of expected convergent mutations from our model
		for z in range(iterations): # number of iterations to run
			actual = sim_unequal_mutations(n, l, m, kappa, phi) # number of actual convergent mutations from simulation
			total_cms += actual
			data[z+1] = (total_cms)/(z+1)
		with open(('better_sim_data_' + str(l) + '_' + str(m) + '.csv'), 'w', newline = '') as f: # writes the data to a .csv file
		    writer = csv.writer(f)
		    writer.writerow(['sim_num',('average_cms_' + str(m))]) # column headers
		    writer.writerows(data.items()) 
		    writer.writerow([str(better_expected)])
		    writer.writerow([((m)**2 * ((kappa**2 + 1/2)/(kappa + 1)**2) * l)]) # writes the new 'mu^2 / 3'
		print('expected: ' + str(((m)**2 * ((kappa**2 + 1/2)/(kappa + 1)**2) * l)) + ' better_expected: ' + str(better_expected)) 
