
from model import expected_cms_given_m
import csv

L = [1000]
mu = [1/1000]
kappa = [0.50,1.0,2.0,3.0,4.0,5.0,6.0]
phi = [0.25,0.50,0.75,1.00]
generations = 300

for l in L:
	for m in mu:
		for k in kappa:
			for p in phi:
				with open(('model_c_vs_m_data_' + str(l) + '_' + str(m) + '_' + str(k) + '_' + str(p) + '.csv'), 'w', newline = '') as f:
					writer = csv.writer(f)
					writer.writerow(['Number of Mutations on Each Strand', 'Expected Number of Convergent Mutations', 'Length', 'mu', 'kappa', 'phi'])
					for g in range(generations):
						writer.writerow([g+1, expected_cms_given_m(l,m,g+1,k,p), l, m, k, p])
						print(g)