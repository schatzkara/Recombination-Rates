from model import expected_cms
from model import expected_cms_with_mg
import csv

L = [1000]
generations = [300]
mu = [1/1000]
# GC = [1.0]
kappa = [1.0, 2.0, 3.0]
phi = 1/2
iterations = 1000

for l in L:
	for g in generations:
		for m in mu:
			# for gc in GC:
			for k in kappa:
				with_mg = g*[None]
				with_mustar = g*[None] # index = g - 1
				for i in range(g):
					# print(1/l)
					# print(k)
					# print(phi)
					# print(i+1)
					cms_with_mg = expected_cms_with_mg(l, m, k, phi, i+1) # THIS DOES IT WITH M^G # params: L, mu, kappa, phi, generations
					cms_with_mustar = expected_cms(l,m*(i+1),k,phi) # THIS DOES IT WITH MU* # params: L,mu,kappa,phi
					with_mg[i] = cms_with_mg
					with_mustar[i] = cms_with_mustar
					print(i)
				with open(('model_mg_vs_mustar_data_' + str(l) + '_' + str(g) + '_' + str(m) + '_' + str(k) + '.csv'), 'w', newline = '') as f:
					writer = csv.writer(f)
					writer.writerow(['Generation', 'Expected CMs using M^g', 'Expected CMs using Mu*']) # column headers
					data = [list(range(1,g+1)), with_mg, with_mustar]
					data = zip(*data) # formats the data so it can be written row by row and produce two long columns
					writer.writerows(data)