from functions import expected_cms 
from mutation_sim import sim_mutations
import numpy as np
import csv

n=2
mu = [.05, .1, .25, .5, .75]

for L in range(4,16):
	for m in mu:
		total_cms = 0
		data = {}
		expected = expected_cms(m, L)
		for z in range(10):
			# expected = expected_cms(m, L)
			actual = sim_mutations(n, L, m)
			total_cms += actual
			data[z+1] = (total_cms)/(z+1)
		with open(('sim_data_' + str(L) + '_' + str(m) + '.csv'), 'w', newline = '') as f:
		    writer = csv.writer(f)
		    writer.writerow(['sim_num',('average_cms_' + str(m))])
		    writer.writerows(data.items())
		    writer.writerow([str(expected)])
		print('expected: ' + str(expected)) 
# print(data)


# expected = {}
# actual = {}
# differences = {}

# for z in range(10):
# # for n in range(2):
# 	for L in range(3000, 3001):
# 		# while 0.01 <= mu <= 0.05:
# 		expected['overlaps'] = np.longdouble(expected_overlaps(int(mu*L), int(mu*L), L))
# 		expected['matches'] = np.longdouble(expected_matches(int(mu*L), int(mu*L), L))
# 		actual = sim_mutations(n, L, mu)
# 		print('expected: ')
# 		print(expected)
# 		print('actual: ')
# 		print(actual)
# 		differences['overlaps'] = expected['overlaps'] - actual['overlaps']
# 		differences['matches'] = expected['matches'] - actual['matches']
# 		# print('the error for ' + str(n) + ' strains of length ' + str(L) + ' with mu = ' + str(mu))
# 		# print(differences)
# 		# mu += 0.01
# # mu = 0.01