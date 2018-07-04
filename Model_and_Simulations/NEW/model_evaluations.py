#! python 3

from model import expected_cms_given_m
from model import expected_idp
from model import expected_cms
from model import expected_cms_with_mg
import csv

def c_vs_m_evaluation(L, mu, generations, kappa, phi):
	for l in L:
		for m in mu:
			for gen in generations:
				for k in kappa:
					for p in phi:
						with open(('model_c_vs_m_data_' + str(l) + '_' + str(m) + '_' + str(gen) + '_' + str(k) + '_' + str(p) + '.csv'), 'w', newline = '') as f:
							writer = csv.writer(f)
							writer.writerow(['Generation', 'Number of Mutations on Each Strand', 'E[c]', 'L', 'mu', 'kappa', 'phi'])
							for g in range(gen):
								writer.writerow([g+1, int(m*l*g), expected_cms_given_m(l,m,g+1,k,p), l, m, k, p])
								print(g)

def c_vs_id_evaluation(L, mu, generations, kappa, phi):
	for l in L: # iterates over all genome lengths desired
		for m in mu:
			for g in generations: # iterates over all numbers of generations desired
				# for gc in GC:
				for k in kappa: # iterates over all values of kappa desired
					for p in phi: # iterates over all values of phi desired
						idps = g*[None] # list to hold the identity percentage after each generation; index = g - 1
						cms = g*[None] # list to hold the number of convergent mutations after each generation; index = g - 1
						for i in range(g): # performs calcuations for each generation
							idp = expected_idp(m, k, p, i+1) # gets the expected identity percentage from the model for the current generation
							expected_cms = expected_cms_given_m(l, m, i+1, k, p) # gets the expected number of convergent mutations from the model for the current generation
							idps[i] = idp
							cms[i] = expected_cms
							print(i)
						with open(('model_c_vs_id_data_' + str(l) + '_' + str(g) + '_' + str(k) + '_' + str(p) + '.csv'), 'w', newline = '') as f: # opens the file that the data will be written to
							writer = csv.writer(f)
							writer.writerow(['Mutations per Strand', 'E[ID%]', 'E[c]', 'L', 'mu', 'generations', 'kappa', 'phi']) # column headers
							data = [list(range(int(m*l),int(m*l*g),int(m*l))), idps, cms, g*[l], g*[m], g*[g], g*[k], g*[p]] # the lists that will form the columns of the output file
							data = zip(*data) # formats the data so it can be written row by row and produce two long columns
							writer.writerows(data)

def mg_vs_mustar_evaluation(L, generations, mu, kappa, phi):
	for l in L: # iterates over all genome lengths desired
		for g in generations: # iterates over all numbers of generations desired
			for m in mu: # iterates over all mutation rates desired
				# for gc in GC:
				for k in kappa: # iterates over all values of kappa desired
					for p in phi: # iterates over all values of phi desired
						with_mg = g*[None] # list of the model values using M^g for each generation; index = g - 1
						with_mustar = g*[None] # list of the model values using mu* for each generation; index = g - 1
						for i in range(g): # performs calculations for each generation
							cms_with_mg = expected_cms_with_mg(l, m, k, p, i+1) # THIS DOES IT WITH M^G # params: L, mu, kappa, phi, generations
							cms_with_mustar = expected_cms(l, m*(i+1), k, p) # THIS DOES IT WITH MU* # params: L, mu, kappa, phi
							with_mg[i] = cms_with_mg
							with_mustar[i] = cms_with_mustar
							print(i)
						with open(('model_mg_vs_mustar_data_' + str(l) + '_' + str(g) + '_' + str(m) + '_' + str(k) + '_' + str(p) + '.csv'), 'wb') as f: # , newline = '') as f: # opens the file that the data will be written to
							writer = csv.writer(f)
							writer.writerow(['Generation', 'E[c] using M^g', 'E[c] using Mu*', 'L', 'mu', 'kappa', 'phi']) # column headers
							data = [list(range(1,g+1)), with_mg, with_mustar, g*[l], g*[m], g*[k], g*[p]] # the three lists that will form the columns of the output file
							data = zip(*data) # formats the data so it can be written row by row and produce two long columns
							writer.writerows(data)
