from functions import expected_cms 
from mutation_sim import sim_mutations
import numpy as np
# import xlwt

n=2
mu = .5
L = 40
data = {}
# book = xlwt.Workbook()
# sheet = book.add_sheet(sheet)
# variables = [x, y, z]
# x_desc = 'Simulation Number'
# y_desc = 'Average Number of CMs'
# desc = [x_desc, y_desc]
# col1_name = 'SimNum'
# col2_name = 'CMs'

total_cms = 0
# for z in range(1):
expected = expected_cms(mu, L)
actual = sim_mutations(n, L, mu)
total_cms += actual
data[1] = (total_cms)/(1)
# sheet.write()
	# print('trial ' + str(z+1))
	# print('the actual cms: ' + str(actual))

print('expected: ' + str(expected)) 
print(data)

# book.save(sim.csv)


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