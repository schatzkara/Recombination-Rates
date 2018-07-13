#! python 3

# from model import expected_cms_given_m
# from model import expected_idp
# import csv
from model_evaluations import c_vs_id_evaluation

L = [1000]
mu = [1/1000]
generations = [300]
# GC = [1.0]
kappa = [0.50,1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
phi = [0.25, 0.50, 0.75, 1.00]
# iterations = 1000

c_vs_id_evaluation(L, mu, generations, kappa, phi)

# for l in L:
# 	for g in generations:
# 		# for gc in GC:
# 		for k in kappa:
# 			idps = g*[None]
# 			cms = g*[None] # index = g - 1
# 			for i in range(g):
# 				print(1/l)
# 				print(k)
# 				print(phi)
# 				print(i+1)
# 				idp = expected_idp(1/l, k, phi, i+1)
# 				expected_cms = expected_cms_given_m(l,i+1,k,phi)
# 				idps[i] = idp
# 				cms[i] = expected_cms
# 				print(i)
# 			with open(('model_id_vs_cm_data_' + str(l) + '_' + str(k) + '.csv'), 'w', newline = '') as f:
# 				writer = csv.writer(f)
# 				writer.writerow(['Number of Mutations per Strand', 'Expected ID%', 'Expected CMs']) # column headers
# 				data = [list(range(1,g+1)), idps, cms]
# 				data = zip(*data) # formats the data so it can be written row by row and produce two long columns
# 				writer.writerows(data)