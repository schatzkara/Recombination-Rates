#! python 3

# script to run each of the simulations according to use input

from model import expected_cms 
from mutation_sim import sim_equal_mutations
from mutation_sim import sim_unequal_mutations
from mutation_sim import efficient_sim_unequal_mutations
from SNPs_sim import sim_snps
from SNPs_sim import sim_id_percent
from data_averages import get_SNP_sim_averages
from data_averages import get_ID_sim_averages
import csv
import os.path

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
					data = {} # dictionary for data over all iterations (key: iteration number, value: average convergent mutations of all iterations up to that point)
					model_expected = expected_cms(l, m, k, p) # better_expected_cms(m, l, kappa, phi) # number of expected convergent mutations from our model
					for z in range(iterations): # number of iterations to run
						if n == 2:
							actual = efficient_mutations_sim(n, l, m, k, p) # number of actual convergent mutations from simulation
						else: 
							actual = mutations_sim(n, l, m, k, p)
						total_cms += actual
						data[z+1] = (total_cms)/(z+1)
					with open(('better_sim_data_' + str(l) + '_' + str(m) + '.csv'), 'w', newline = '') as f: # writes the data to a .csv file
					    writer = csv.writer(f)
					    writer.writerow(['sim_num',('average_cms_' + str(m))]) # column headers
					    writer.writerows(data.items()) 
					    writer.writerow([str(model_expected)])
					    # writer.writerow([((m)**2 * ((k**2 + 1/2)/(k + 1)**2) * l)]) # writes the expected value
					print('expected: ' + str(model_expected))
					# print('expected: ' + str(((m)**2 * ((k**2 + 1/2)/(k + 1)**2) * l)) + ' model_expected: ' + str(model_expected)) 

def run_SNP_sim(L, generations, kappa, phi, iterations, one_file):
	for l in L: # iterates over every length desired
		for g in generations: # iterates over every number of SNPs desired
			for gc in GC: # iterates over every GC% desired
				for k in kappa: # iterates over every kappa desired
					for p in phi: # iterates over every phi desired
						for z in range(iterations): # runs it for the desired iterations
							current_path = os.path.dirname(os.path.abspath(__file__))
							save_path = current_path + 'SNP Sim/L = ' + str(l) + '/GC% = ' + str(gc) + '/kappa = ' + str(k) + '/phi = ' + str(p)
							# save_path = 'C:/Users/Owner/Documents/UNCG REU/Project/Recombination-Rates/Data/SNP Sim/L = ' + str(l) + '/GC% = ' + str(gc) + '/kappa = ' + str(int(k))
							file_name = 'SNPs_sim_data_' + str(l) + '_' + str(g) + '_' + str(gc) + '_' + str(k) + '_' + str(p) + '_' + '{0:04}'.format(z+1)
							full_name = os.path.join(save_path, file_name + '.csv')   
							with open((full_name), 'w', newline = '') as f: # writes the data to a .csv file
							    writer = csv.writer(f)
							    writer.writerow(['Mutations on Each Strain', 'CMs']) # column headers
							    data = [list(range(1,g+1)), (sim_snps(l, g, gc, k, p))] # one variable is the number of SNPs and the other is the number of convergent mutations
							    data = zip(*data) # formats the data so it can be written row by row and produce two long columns
							    writer.writerows(data) # writes the data
	if(one_file):
		get_SNP_sim_averages(L, generations, GC, kappa, phi, save_path)

def run_ID_sim(L, generations GC, kappa, phi, iterations, one_file):
	for l in L: # iterates over every length desired
		for g in generations: # iterates over every number of SNPs desired
			for gc in GC: # iterates over every GC% desired
				for k in kappa: # iterates over every kappa desired
					for p in phi: # iterates over every phi desired
						for z in range(iterations): # runs it for the desired iterations
							current_path = os.path.dirname(os.path.abspath(__file__))
							save_path = current_path + 'ID Sim/L = ' + str(l) + '/GC% = ' + str(gc) + '/kappa = ' + str(k) + '/phi = ' + str(p)
							# save_path = 'C:/Users/Owner/Documents/UNCG REU/Project/Recombination-Rates/Data/ID Sim/L = ' + str(l) + '/GC% = ' + str(gc) + '/kappa = ' + str(int(k))
							file_name = 'ID_sim_data_' + str(l) + '_' + str(g) + '_' + str(gc) + '_' + str(k) + '_' + str(p) + '_' + '{0:04}'.format(z+1)
							full_name = os.path.join(save_path, file_name + '.csv')   
							with open((full_name), 'w', newline = '') as f: # writes the data to a .csv file
							    writer = csv.writer(f)
							    writer.writerow(['Mutations on Each Strain', 'ID%', 'CMs']) # column headers
							    sim = sim_id_percent(l, g, gc, k, p)
							    for key in sim.keys():
							    	writer.writerow([key+1, sim[key][0], sim[key][1]]) # writes the data
	if(one_file):
		get_ID_sim_averages(L, generations, GC, kappa, phi, save_path)