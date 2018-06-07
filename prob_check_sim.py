from functions import expected_overlaps
from functions import expected_matches
from mutation_sim import sim_mutations

mu = 0.01
x = 3
expected = {}
actual = {}
differences = {}

# for z in range(10):
for n in range(100, 105):
	for L in range(100, 105):
		while 0.01 <= mu <= 0.05:
			expected['overlaps'] = expected_overlaps(int(mu*L), int(mu*L), L)
			expected['matches'] = expected_matches(int(mu*L), int(mu*L), L)
			actual = sim_mutations(n, L, mu)
			differences['overlaps'] = expected['overlaps'] - actual['overlaps']
			differences['matches'] = expected['matches'] - actual['matches']
			print('the error for ' + str(n) + ' strains of length ' + str(L) + ' with mu = ' + str(mu))
			print(differences)
			mu += 0.01
	mu = 0.01