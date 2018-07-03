
# script to run the mutation simulation over multiple iterations and compare the actual number of convergent mutations to that from our model
# time complexity for each mu,L combination: O(n^4), where n is L

from run_sims import run_mutation_sim
# from simulations import efficient_mutation_sim
# from model import expected_cms

# params:
# 	n (int) = number of DNA strands
# 	L (list of ints) = list of DNA strand lengths (units: base pairs)
# 	mu (list of floats) = list of mutation rates (units: mutations per base pair per generation)
# 	kappa (float) = ratio of transitions to transversions
# 	phi (float) = probabiltiy of transversion to its complementary base pair
# 	iterations (int) = number of iterations to run
# return: a .csv file will be produced in the directory where this file is located with the data from all iterations and the expected value

# enter parameters here
L = 1000
mu = [1/1000]
generations = 300
kappa = [1.0,2.0,3.0,4.0,5.0,6.0] 
phi = 0.5
iterations = 10

run_mutation_sim(2, L, mu, generations, kappa, phi, iterations)

# for m in mu:
    # for k in kappa:
        # for i in range(iterations):
            # # we want to output the results of these two functions
            # efficient_mutation_sim(L, m, generations, k, phi)
            # expected_cms(L, m, k, phi)
