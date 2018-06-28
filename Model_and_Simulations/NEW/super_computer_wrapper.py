from simulations import mutations_sim
from simulations import efficient_mutations_sim
from model import expected_cms
import csv

n = 2
L = [10000]
mu = [.0001, .001, .01, .1]
kappa = [1.0,2.0,3.0]
phi = [1/2]
iterations = 1000

for l in L:
	for m in mu:
		for k in kappa:
			for p in phi:
				data = [] # list to hold the data values from each iteration of the simulation
				model_expected = expected_cms(l, m, k, p) # expected number of convergent mutations from our model
				for x in range(iterations):
					if n == 2: # when there are only two DNA strands, the sim can be more efficient
						data.append(efficient_mutations_sim(n, l, m, k, p)) # number of actual convergent mutations from simulation
					else: 
						data.append(mutations_sim(n, l, m, k, p))
				with open(('mutation_sim_data_ ' + str(l) + '_' + str(m) + '_' + str(k) + '_' + str(p) + '.csv'), 'w', newline = '') as f: # writes all the data to a .csv file 
					writer = csv.writer(f)
					writer.writerow(['Iteration Number', 'CMs', 'L', 'mu', 'kappa', 'phi', 'Expected CMs'])
					data_columns = [list(range(1,iterations+1)), data, iterations*[l], iterations*[m], iterations*[k], iterations*[p], iterations*[model_expected]] # all the columns for the .csv file
					data_rows = zip(*data_columns) # formats the data so it can be written as rows, not columns
					writer.writerows(data_rows)
				print('expected for ' + str(l) + ',' + str(m) + ',' + str(k) + ',' + str(p) + ': ' + str(model_expected))
