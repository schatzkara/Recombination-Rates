#! python3

from functions import expected_cms 
from mutation_sim import sim_equal_mutations
import numpy as np
import csv

# time complexity: O(n^6); technically it's 2n^6
n=2
mu = [.2, .4, .6, .8]

for L in range(41,101):
	for m in mu:
		total_cms = 0
		data = {}
		expected = expected_cms(m, L)
		for z in range(1000):
			actual = sim_equal_mutations(n, L, m)
			total_cms += actual
			data[z+1] = (total_cms)/(z+1)
		with open(('sim_data_' + str(L) + '_' + str(m) + '.csv'), 'w', newline = '') as f:
		    writer = csv.writer(f)
		    writer.writerow(['sim_num',('average_cms_' + str(m))])
		    writer.writerows(data.items())
		    writer.writerow([str(expected)])
		    writer.writerow([((m**2)/3) * L])
		print('expected: ' + str(expected)) 
# print(data)
