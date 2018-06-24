#! python3

from run_sims import run_mutation_sim_with_equal_chances
# from functions import expected_cms 
# from mutation_sim import sim_equal_mutations
# import numpy as np
# import csv

# time complexity: O(n^6); technically it's 2n^6
n = 2
L = list(range(50,101))
mu = [.2, .4, .6, .8]

run_mutation_sim_with_equal_chances(n, L, mu, iterations)

def run_mutation_sim_with_equal_chances(n, L, mu, iterations):
	for l in L:
		for m in mu:
			total_cms = 0
			data = {}
			expected = expected_cms(m, l)
			for z in range(iterations):
				actual = sim_equal_mutations(n, l, m)
				total_cms += actual
				data[z+1] = (total_cms)/(z+1)
			with open(('sim_data_' + str(l) + '_' + str(m) + '.csv'), 'w', newline = '') as f:
			    writer = csv.writer(f)
			    writer.writerow(['sim_num',('average_cms_' + str(m))])
			    writer.writerows(data.items())
			    writer.writerow([str(expected)])
			    # writer.writerow([((m**2)/3) * l])
			print('expected: ' + str(expected)) 

# for l in L:
# 	for m in mu:
# 		total_cms = 0
# 		data = {}
# 		expected = expected_cms(m, l)
# 		for z in range(1000):
# 			actual = sim_equal_mutations(n, l, m)
# 			total_cms += actual
# 			data[z+1] = (total_cms)/(z+1)
# 		with open(('sim_data_' + str(l) + '_' + str(m) + '.csv'), 'w', newline = '') as f:
# 		    writer = csv.writer(f)
# 		    writer.writerow(['sim_num',('average_cms_' + str(m))])
# 		    writer.writerows(data.items())
# 		    writer.writerow([str(expected)])
# 		    writer.writerow([((m**2)/3) * l])
# 		print('expected: ' + str(expected)) 
# # print(data)
