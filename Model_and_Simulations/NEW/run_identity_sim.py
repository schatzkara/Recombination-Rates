#! python 3

# script to run the ID% simulation according to the users inputs
# time complexity for each L, GC, kappa combination: O(n^3)

from run_sims import run_identity_sim


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
GC = [0.5]
kappa = [1.86836732388]
# kappa = [6.1,6.2,.3,6.4,6.5,6.6,6.7,6.8,6.9,7.0,7.1,7.2,7.3,7.4,7.5,7.6,7.7,7.8,7.9,8.0,8.1,8.2,8.3,8.4,8.5,8.6,8.7,8.8,8.9,9.0,9.1,9.2,9.3,9.4,9.5,9.6,9.7,9.8,9.9,10.0]
# kappa = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0]
phi = [0.50]
iterations = 1000
one_file = True

run_identity_sim(L, mu, generations, GC, kappa, phi, iterations, one_file)

# L = [1000]
# mu = [1/1000]
# generations = [300]
# GC = [0.8,0.9,1.0]
# kappa = [1.0,2.0,3.0]
# phi = [0.50]
# iterations = 1000
# one_file = True

# run_identity_sim(L, mu, generations, GC, kappa, phi, iterations, one_file)
