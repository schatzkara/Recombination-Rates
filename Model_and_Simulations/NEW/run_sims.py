#! python 3

# script to run each of the simulations according to use input

from model import expected_cms 
from simulations import mutations_sim
from simulations import efficient_mutations_sim
from simulations import mutations_over_generations_sim
from simulations import id_percent_sim
from simulations import id_matrix_sim
# from mutation_sim import sim_equal_mutations
# from mutation_sim import sim_unequal_mutations
# from mutation_sim import efficient_sim_unequal_mutations
# from SNPs_sim import sim_snps
# from SNPs_sim import sim_id_percent
from data_averages import get_SNP_sim_averages
from data_averages import get_ID_sim_averages
import csv
import os

# function to run the mutation simulation over multiple iterations and compare the actual number of convergent mutations to that from our model
# time complexity for each mu,L,kappa,phi combination: O(n^4), where n is L
# params:
# 	n (int) = number of DNA strands
# 	L (list of ints) = list of DNA strand lengths (units: base pairs)
# 	mu (list of floats) = list of mutation rates (units: mutations per base pair per generation)
# 	kappa (list of floats) = list of ratios of transitions to transversions
# 	phi (list of floats) = list of probabilities of transversion to its complementary base pair
# 	iterations (int) = number of iterations to run
# return: a .csv file will be produced in the directory where this file is located with the data from all iterations and the expected value
def run_mutation_sim(n, L, mu, kappa, phi, iterations):
	for l in L: # iterates over every length desired
		for m in mu: # iterates over every mu desired
			for k in kappa: # iterates over every kappa desired
				for p in phi: # iterates over every phi desired
					total_cms = 0 # counter for total convergent mutations over all iterations
					data = iterations*[None] # dictionary for data over all iterations (key: iteration number, value: average convergent mutations of all iterations up to that point)
					model_expected = expected_cms(l, m, k, p) # better_expected_cms(m, l, kappa, phi) # number of expected convergent mutations from our model
					for z in range(iterations): # number of iterations to run
						if n == 2:
							actual = efficient_mutations_sim(n, l, m, k, p) # number of actual convergent mutations from simulation
						else: 
							actual = mutations_sim(n, l, m, k, p)
						total_cms += actual
						data[z] = (total_cms)/(z+1)
					with open(('mutation_sim_' + str(l) + '_' + str(m) + '_' + str(k) + '_' + str(p) + '.csv'), 'w', newline = '') as f: # writes the data to a .csv file
					    writer = csv.writer(f)
					    writer.writerow(['Iteration Number', 'Average CMs', 'L', 'mu', 'kappa', 'phi', 'Expected CMs']) # column headers
					    data_columns = [(list(range(1,iterations+1))), data, iterations*[l], iterations*[m], iterations*[k], iterations*[p], iterations*[model_expected]]
					    data_columns = zip(*data_columns)
					    writer.writerows(data_columns)
					    # writer.writerows(data.items()) 
					    # writer.writerow([str(model_expected)])
					    # writer.writerow([((m)**2 * ((k**2 + 1/2)/(k + 1)**2) * l)]) # writes the expected value
					print('expected: ' + str(model_expected))
					# print('expected: ' + str(((m)**2 * ((k**2 + 1/2)/(k + 1)**2) * l)) + ' model_expected: ' + str(model_expected)) 

def run_SNP_sim(L, generations, GC, kappa, phi, iterations, one_file):
	current_path = os.path.dirname(os.path.abspath(__file__))
	for l in L: # iterates over every length desired
		for g in generations: # iterates over every number of SNPs desired
			for gc in GC: # iterates over every GC% desired
				for k in kappa: # iterates over every kappa desired
					for p in phi: # iterates over every phi desired
						for z in range(iterations): # runs it for the desired iterations
							# current_path = os.path.dirname(os.path.abspath(__file__))
							save_path = current_path + '/SNP Sim/L = ' + str(l) + '/GC% = ' + str(gc) + '/kappa = ' + str(k) + '/phi = ' + str(p)
							if not os.path.exists(save_path):
								os.makedirs(save_path)
							# save_path = 'C:/Users/Owner/Documents/UNCG REU/Project/Recombination-Rates/Data/SNP Sim/L = ' + str(l) + '/GC% = ' + str(gc) + '/kappa = ' + str(int(k))
							file_name = 'SNPs_sim_data_' + str(l) + '_' + str(g) + '_' + str(gc) + '_' + str(k) + '_' + str(p) + '_' + '{0:04}'.format(z+1)
							full_name = os.path.join(save_path, file_name + '.csv')   
							with open((full_name), 'w', newline = '') as f: # writes the data to a .csv file
							    writer = csv.writer(f)
							    writer.writerow(['Mutations on Each Strain', 'CMs', 'L', 'GC%', 'kappa', 'phi']) # column headers
							    data = [list(range(1,g+1)), (mutations_over_generations_sim(l, g, gc, k, p)), g*[l], g*[gc], g*[k], g*[p]] # one variable is the number of SNPs and the other is the number of convergent mutations
							    data = zip(*data) # formats the data so it can be written row by row and produce two long columns
							    writer.writerows(data) # writes the data
	if(one_file):
		new_path = current_path + '/SNP Sim'
		get_SNP_sim_averages(L, generations, GC, kappa, phi, save_path, new_path)

def run_ID_sim(L, generations, GC, kappa, phi, iterations, one_file):
	current_path = os.path.dirname(os.path.abspath(__file__))
	for l in L: # iterates over every length desired
		for g in generations: # iterates over every number of SNPs desired
			for gc in GC: # iterates over every GC% desired
				for k in kappa: # iterates over every kappa desired
					for p in phi: # iterates over every phi desired
						for z in range(iterations): # runs it for the desired iterations
							save_path = current_path + '/ID Sim/L = ' + str(l) + '/GC% = ' + str(gc) + '/kappa = ' + str(k) + '/phi = ' + str(p)
							if not os.path.exists(save_path):
								os.makedirs(save_path)
							# save_path = 'C:/Users/Owner/Documents/UNCG REU/Project/Recombination-Rates/Data/ID Sim/L = ' + str(l) + '/GC% = ' + str(gc) + '/kappa = ' + str(int(k))
							file_name = 'ID_sim_data_' + str(l) + '_' + str(g) + '_' + str(gc) + '_' + str(k) + '_' + str(p) + '_' + '{0:04}'.format(z+1)
							full_name = os.path.join(save_path, file_name + '.csv')   
							with open((full_name), 'w', newline = '') as f: # writes the data to a .csv file
							    writer = csv.writer(f)
							    writer.writerow(['Mutations on Each Strain', 'ID%', 'CMs', 'L', 'GC%', 'kappa', 'phi']) # column headers
							    sim = id_percent_sim(l, g, gc, k, p)
							    for key in sim.keys():
							    	writer.writerow([key+1, sim[key][0], sim[key][1], l, gc, k, p]) # writes the data
	if(one_file):
		new_path = current_path + '/ID Sim'
		get_ID_sim_averages(L, generations, GC, kappa, phi, save_path, new_path)

def run_ID_matrix_sim(N, L, generations, mu, kappa, phi, iterations):
	for n in N:
		for l in L:
			for g in generations:
				for m in mu:
					for k in kappa:
						for p in phi:
							with open((id_matrix_sim_data_ + str(n) + '_' + str(l) + '_' + str(m) + '_' + str(k) + '_' + str(p) + '.csv'), 'w', newline = '') as f:
								writer = csv.writer(f)
								c_q = (n-2)*[] # index = q - 3
								for q in range(3,n+1):
									c_q[q-3] = 'c_' + 'str(q)
								writer.writerow(['Iteration Number', c_q, 'n', 'L', 'mu', 'kappa', 'phi'])
								values = i*[None]
								for i in range(iterations):
									values[i] = id_matrix_sim(n, l, g, m, k, p)
									writer.writerow([(i+1),values[i], n,l,m,k,p])
								# data = [list(range(1,iterations+1)), values, iterations*[n], iterations*[l], iterations*[m], iterations*[k], iterations*[p]]
								# data = zip(*data)
								# writer.writerows(data)
								
								
								