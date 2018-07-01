
# script to run the mutation simulation over multiple iterations and compare the actual number of convergent mutations to that from our model
# time complexity for each mu,L combination: O(n^4), where n is L

from simulations import efficient_mutation_sim
from model import expected_cms

# params:
# 	n (int) = number of DNA strands
# 	L (list of ints) = list of DNA strand lengths (units: base pairs)
# 	mu (list of floats) = list of mutation rates (units: mutations per base pair per generation)
# 	kappa (float) = ratio of transitions to transversions
# 	phi (float) = probabiltiy of transversion to its complementary base pair
# 	iterations (int) = number of iterations to run
# return: a .csv file will be produced in the directory where this file is located with the data from all iterations and the expected value

# enter parameters here
L = 10000
mu = [0.001, 0.01, 0.1]
generations = 1
kappa = [1.0,2.0,3.0,4.0,5.0,6.0] 
phi = 0.5
iterations = 10

for m in mu:
    for k in kappa:
        for i in range(iterations):
            # we want to output the results of these two functions
            efficient_mutation_sim(L, m, generations, k, phi)
            expected_cms(L, m, k, phi)
