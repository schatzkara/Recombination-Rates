#! python 3

from model import expected_cms_given_m
from model import expected_idp
from model import expected_cms
from model import expected_cms_with_mg
import csv

def cm_vs_id_evaluation(L, generations, kappa, phi):

	# L = [1000]
	# generations = [300]
	# # GC = [1.0]
	# kappa = [2.0, 3.0]
	# phi = 1/2
	# # iterations = 1000

	for l in L: # iterates over all genome lengths desired
		for g in generations: # iterates over all numbers of generations desired
			# for gc in GC:
			for k in kappa: # iterates over all values of kappa desired
				for p in phi: # iterates over all values of phi desired
					idps = g*[None] # list to hold the identity percentage after each generation; index = g - 1
					cms = g*[None] # list to hold the number of convergent mutations after each generation; index = g - 1
					for i in range(g): # performs calcuations for each generation
						idp = expected_idp(1/l, k, p, i+1) # gets the expected identity percentage from the model for the current generation
						expected_cms = expected_cms_given_m(l, i+1, k, p) # gets the expected number of convergent mutations from the model for the current generation
						idps[i] = idp
						cms[i] = expected_cms
						print(i)
					with open(('model_id_vs_cm_data_' + str(l) + '_' + str(g) + '_' + str(k) + '_' + str(p) + '.csv'), 'w', newline = '') as f: # opens the file that the data will be written to
						writer = csv.writer(f)
						writer.writerow(['Number of Mutations per Strand', 'Expected ID%', 'Expected CMs']) # column headers
						data = [list(range(1,g+1)), idps, cms] # the three lists that will form the columns of the output file
						data = zip(*data) # formats the data so it can be written row by row and produce two long columns
						writer.writerows(data)

def mg_vs_mustar_evaluation(L, generations, mu, kappa, phi):

	# L = [1000]
	# generations = [300]
	# mu = [1/1000]
	# # GC = [1.0]
	# kappa = [1.0, 2.0, 3.0]
	# phi = [1/2]
	# iterations = 1000

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
						with open(('model_mg_vs_mustar_data_' + str(l) + '_' + str(g) + '_' + str(m) + '_' + str(k) + '_' + str(p) + '.csv'), 'w', newline = '') as f: # opens the file that the data will be written to
							writer = csv.writer(f)
							writer.writerow(['Generation', 'Expected CMs using M^g', 'Expected CMs using Mu*']) # column headers
							data = [list(range(1,g+1)), with_mg, with_mustar] # the three lists that will form the columns of the output file
							data = zip(*data) # formats the data so it can be written row by row and produce two long columns
							writer.writerows(data)
