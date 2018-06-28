from run_sims import run_mutation_sim_super_computer

L = [10000]
mu = [.0001, .001, .01, .1]
kappa = [1.0,2.0,3.0]
phi = [1/2]
iterations = 1000

for l in L:
	for m in mu:
		for k in kappa:
			for p in phi:
				data = []
				model_expected = expected_cms(l, m, k, p)
				for x in range(iterations):
					if n == 2:
						data.append(efficient_mutations_sim(n, l, m, k, p)) # number of actual convergent mutations from simulation
					else: 
						data.append(mutations_sim(n, l, m, k, p))
				with open((mutation_sim_data.csv), 'w', newline = '') as f:
					writer = csv.writer(f)
					writer.writerow(['Iteration Number', 'CMs', 'L', 'mu', 'kappa', 'phi', 'Expected CMs'])
					data_columns = [list(range(1,iterations=1)), data, iterations*[l], , iterations*[m], iterations*[k], iterations*[p], iterations*[model_expected]]
					data_rows = zip(*data_columns)
					writer.writerows(data_rows)
				print('expected for ' + str(l) + ',' + str(m) + ',' + str(k) + ',' + str(p) + ': ' + str(model_expected))
