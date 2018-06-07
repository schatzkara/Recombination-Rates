from functions import expected_overlaps
from functions import expected_matches
from mutation_sim import sim_mutations
import numpy as np

mu = 0.05
x = 3
n = 2
expected = {}
actual = {}
differences = {}

for z in range(10):
# for n in range(2):
	for L in range(3000, 3001):
		# while 0.01 <= mu <= 0.05:
		expected['overlaps'] = np.longdouble(expected_overlaps(int(mu*L), int(mu*L), L))
		expected['matches'] = np.longdouble(expected_matches(int(mu*L), int(mu*L), L))
		actual = sim_mutations(n, L, mu)
		print('expected: ')
		print(expected)
		print('actual: ')
		print(actual)
		differences['overlaps'] = expected['overlaps'] - actual['overlaps']
		differences['matches'] = expected['matches'] - actual['matches']
		# print('the error for ' + str(n) + ' strains of length ' + str(L) + ' with mu = ' + str(mu))
		# print(differences)
		# mu += 0.01
# mu = 0.01