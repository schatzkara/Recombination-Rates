
# from model import expected_cms_given_m
from model_evaluations import c_vs_m_evaluation
import csv

L = [1000]
mu = [1/1000]
kappa = [0.50,1.0,2.0,3.0,4.0,5.0,6.0]
phi = [0.50]
generations = [300]

c_vs_m_evaluation(L, mu, generations, kappa, phi)


# for l in L:
# 	for m in mu:
# 		for k in kappa:
# 			for p in phi:
# 				for gen in generations:
# 					with open(('model_c_vs_m_data_' + str(l) + '_' + str(m) + '_' + str(gen) + '_' + str(k) + '_' + str(p) + '.csv'), 'w', newline = '') as f:
# 						writer = csv.writer(f)
# 						writer.writerow(['Generation', 'Number of Mutations on Each Strand', 'E[c]', 'L', 'mu', 'kappa', 'phi'])
# 						for g in range(gen):
# 							writer.writerow([g+1, int(m*l*g), expected_cms_given_m(l,m,g+1,k,p), l, m, k, p])
# 							print(g)