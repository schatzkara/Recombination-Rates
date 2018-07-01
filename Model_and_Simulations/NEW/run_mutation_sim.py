#! python3
# export PATH=/opt/rh/rh-python36/root/usr/bin:$PATH


# script to run the mutation simulation over multiple iterations and compare the actual number of convergent mutations to that from our model
# time complexity for each mu,L combination: O(n^4), where n is L

from run_sims import run_mutation_sim

# params:
# 	n (int) = number of DNA strands
# 	L (list of ints) = list of DNA strand lengths (units: base pairs)
# 	mu (list of floats) = list of mutation rates (units: mutations per base pair per generation)
# 	kappa (float) = ratio of transitions to transversions
# 	phi (float) = probabiltiy of transversion to its complementary base pair
# 	iterations (int) = number of iterations to run
# return: a .csv file will be produced in the directory where this file is located with the data from all iterations and the expected value

# enter parameters here
n = 2
L = [100] # format for range(): list(range(starting length, ending length+1, increment))
mu = [0.01, 0.1]
generations = [1]
kappa = [0.25, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0] 
phi = [0.0, 0.25, 0.5, 0.75, 1.0] 
iterations = 1000

run_mutation_sim(n, L, mu, generations, kappa, phi, iterations)
