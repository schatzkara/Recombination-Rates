#! python 3

# script to run the generation simulation according to the users inputs
# time complexity for each L, GC, kappa combination: O(n^3)

from run_sims import run_generation_sim

# params:
# 	L (list of ints) = list of DNA strand lengths (units: base pairs)
# 	generations (list of ints) = list of the number of SNPs to introduce to each strain
# 	GC (list of floats) = list of the GC% composition of the ancestral strain
# 	kappa (list of floats) = list of the ratios of transitions to transversions
# 	phi (float) = probabiltiy of transversion to its complementary base pair
# 	iterations (int) = number of iterations to run
# 	one_file (boolean) = True if you just want one output file that contains the averages over all iterations; False if you want one file for each iteration
# return: a .csv file will be produced in the directory where this file is located with the data for a single interation (one .csv file for each iteration)

# enter parameters here
L = [1000]
mu = [1/1000]
generations = [300]
GC = [0.1, 0.2, 0.3, 0.4, 0.5]
kappa = [0.50, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
phi = [0.25, 0.50, 0.75, 1.00]
iterations = 1000
one_file = True

run_generation_sim(L, mu, generations, GC, kappa, phi, iterations, one_file)

