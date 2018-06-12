from efficient_functions import better_expected_cms
from mutation_sim import sim_unequal_mutations
import csv

# time complexity: O(n^6); technically it's 2n^6
n=2
mu = [0.2,0.4]
kappa = 3
phi = 1/2

for L in range(4,10):
	for m in mu:
		total_cms = 0
		data = {}
		better_expected = better_expected_cms(m, L, kappa, phi)
		for z in range(1000):
			actual = sim_unequal_mutations(n, L, m, kappa)
			total_cms += actual
			data[z+1] = (total_cms)/(z+1)
		with open(('better_sim_data_' + str(L) + '_' + str(m) + '.csv'), 'w', newline = '') as f:
		    writer = csv.writer(f)
		    writer.writerow(['sim_num',('average_cms_' + str(m))])
		    writer.writerows(data.items())
		    writer.writerow([str(better_expected)])
		    writer.writerow([((m)**2 * ((kappa**2 + 1/2)/(kappa + 1)**2) * L)])
		print('expected: ' + str(((m)**2 * ((kappa**2 + 1/2)/(kappa + 1)**2) * L)) + ' better_expected: ' + str(better_expected)) 
