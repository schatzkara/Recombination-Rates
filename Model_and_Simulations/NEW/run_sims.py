#! python 3

# script to run each of the simulations according to use input

from model import expected_cms 
from simulations import mutation_sim
from simulations import efficient_mutation_sim
from simulations import generation_sim
from simulations import identity_sim
from simulations import id_matrix_sim
from simulations import c_q_sim
from simulations import mutation_sites_sim
from data_averages import get_generation_sim_averages
from data_averages import get_identity_sim_averages
from data_averages import get_mutation_sites_sim_averages
import csv
import os

# function to run the mutation simulation over multiple iterations and compare the actual number of convergent mutations to that from our model
# time complexity for each mu,L,kappa,phi combination: O(n^4), where n is L
# params:
# 	n (int) = number of DNA strands
# 	L (list of ints) = list of DNA strand lengths (units: base pairs)
# 	mu (list of floats) = list of mutation rates (units: mutations per base pair per generation)
#       generations (list of ints) = list of generations
# 	kappa (list of floats) = list of ratios of transitions to transversions
# 	phi (list of floats) = list of probabilities of transversion to its complementary base pair
# 	iterations (int) = number of iterations to run
# return: a .csv file will be produced in the directory where this file is located with the data from all iterations and the expected value
def run_mutation_sim(n, L, mu, generations, kappa, phi, iterations):
	current_path = os.path.dirname(os.path.abspath(__file__))
	for l in L: # iterates over every length desired
		for m in mu: # iterates over every mu desired
			for g in range(generations+1):
				for k in kappa: # iterates over every kappa desired
					for p in phi: # iterates over every phi desired
						total_cms = 0 # counter for total convergent mutations over all iterations
						data = iterations*[None] # dictionary for data over all iterations (key: iteration number, value: average convergent mutations of all iterations up to that point)
						model_expected = expected_cms(l, m, k, p) # better_expected_cms(m, l, kappa, phi) # number of expected convergent mutations from our model
						for i in range(iterations): # number of iterations to run
							if n == 2:
								actual = efficient_mutation_sim(l, m, g, k, p) # number of actual convergent mutations from simulation
							else: 
								actual = mutation_sim(n, l, m, g, k, p)
							total_cms += actual
							data[i] = (total_cms)/(i+1)
							print(i)
						save_path = current_path + '/Mutation Sim/L = ' + str(l)
						if not os.path.exists(save_path):
							os.makedirs(save_path)
						file_name = 'mutation_sim_' + str(l) + '_' + '{0:04}'.format(m) +  '_' + str(g) + '_' + '{0:04}'.format(k) + '_' + '{0:04}'.format(p)
						full_name = os.path.join(save_path, file_name + '.csv')   
						with open((full_name), 'w', newline = '') as f: # writes the data to a .csv file
							writer = csv.writer(f)
							writer.writerow(['Iteration', 'Average c', 'L', 'mu', 'Generation', 'kappa', 'phi', 'E[c]']) # column headers
							data_columns = [(list(range(1,iterations+1))), data, iterations*[l], iterations*[m], list(range(1,generations+1)), iterations*[k], iterations*[p], iterations*[model_expected]]
							data_row = zip(*data_columns)
							writer.writerows(data_row)
						print('expected: ' + str(model_expected))

# function to run the 
def run_generation_sim(L, mu, generations, GC, kappa, phi, iterations, one_file):
	current_path = os.path.dirname(os.path.abspath(__file__))
	for l in L: # iterates over every length desired
		for m in mu: # iterates over every mutation rate desired
			for g in generations: # iterates over every number of generations desired
				for gc in GC: # iterates over every GC% desired
					for k in kappa: # iterates over every kappa desired
						for p in phi: # iterates over every phi desired
							# save_path = 'D:/BIGG DATA/Generation Sim/all GCs with kappa 0.5,1-6/L = 1000/GC% = ' + str(gc) + '/kappa = ' + str(k) + '/phi = ' + str(p)
							save_path = current_path + '/Generation Sim 2/L = ' + str(l) + '/GC% = ' + str(gc) + '/kappa = ' + str(k) + '/phi = ' + str(p)
							for i in range(iterations): # runs it for the desired iterations
								if not os.path.exists(save_path):
									os.makedirs(save_path)
								file_name = 'generation_sim_' + str(l) + '_' + '{0:04}'.format(m) + '_' + str(g) + '_' + '{0:04}'.format(gc) + '_' + '{0:04}'.format(k) + '_' + '{0:04}'.format(p) + '_' + '{0:04}'.format(i+1)
								full_name = os.path.join(save_path, file_name + '.csv')   
								with open((full_name), 'w', newline = '') as f: # writes the data to a .csv file
									writer = csv.writer(f)
									writer.writerow(['Generation', 'Mutations on Each Strain', 'c', 'L', 'mu', 'GC%', 'kappa', 'phi']) # column headers
									mutations_per_generation = []

									data = [list(range(1,g+1)), list(range(int(m*l),int((m*l*g)+1),int(m*l))), (generation_sim(l, m, g, gc, k, p)), g*[l], g*[m], g*[gc], g*[k], g*[p]] # one variable is the number of SNPs and the other is the number of convergent mutations
									data = zip(*data) # formats the data so it can be written row by row and produce two long columns
									writer.writerows(data) # writes the data
								print(i)
							if(one_file):
								new_path = current_path + '/Generation Sim'
								get_generation_sim_averages(l, m, g, gc, k, p, save_path, new_path)

def run_identity_sim(L, mu, generations, GC, kappa, phi, iterations, one_file):
	current_path = os.path.dirname(os.path.abspath(__file__))
	for l in L: # iterates over every length desired
		for m in mu: # iterates over every mutation rate desired
			for g in generations: # iterates over every number of SNPs desired
				for gc in GC: # iterates over every GC% desired
					for k in kappa: # iterates over every kappa desired
						for p in phi: # iterates over every phi desired
							# save_path = 'D:/BIGG DATA/Identity Sim/all GCs with kappa .5,1-6/L = 1000/GC% = ' + str(gc) + '/kappa = ' + str(k) + '/phi = ' + str(p)
							save_path = current_path + '/Identity Sim/L = ' + str(l) + '/GC% = ' + str(gc) + '/kappa = ' + str(k) + '/phi = ' + str(p)
							for i in range(iterations): # runs it for the desired iterations
								if not os.path.exists(save_path):
									os.makedirs(save_path)
								# save_path = 'C:/Users/Owner/Documents/UNCG REU/Project/Recombination-Rates/Data/ID Sim/L = ' + str(l) + '/GC% = ' + str(gc) + '/kappa = ' + str(int(k))
								file_name = 'identity_sim_' + str(l) + '_' + '{0:04}'.format(m) + '_' +  str(g) + '_' + '{0:04}'.format(gc) + '_' + '{0:04}'.format(k) + '_' + '{0:04}'.format(p) + '_' + '{0:04}'.format(i+1)
								full_name = os.path.join(save_path, file_name + '.csv')   
								with open((full_name), 'w', newline = '') as f: # writes the data to a .csv file
									writer = csv.writer(f)
									writer.writerow(['Generation', 'Mutations on Each Strain', 'ID%', 'c', 'L', 'mu', 'GC%', 'kappa', 'phi']) # column headers
									sim = identity_sim(l, m, g, gc, k, p)
									for key in sim.keys():
										writer.writerow([key, int(m*l*key), sim[key][0], sim[key][1], l, m, gc, k, p]) # writes the data
								print(i)
							if(one_file):
								new_path = current_path + '/Identity Sim'
								get_identity_sim_averages(l, m, g, gc, k, p, save_path, new_path)

def run_id_matrix_sim(N, L, mu, generations, kappa, phi, iterations):
    current_path = os.path.dirname(os.path.abspath(__file__))
    for n in N:
        for l in L:
            for m in mu:
                for g in generations:
                    for k in kappa:
                        for p in phi:
                            save_path = current_path + '/ID Matrix Sim/L = ' + str(l)
                            if not os.path.exists(save_path):
                                os.makedirs(save_path)
                            file_name = 'id_matrix_sim_' + '{0:04}'.format(n) + '_' + str(l) + '_' + '{0:04}'.format(m) + '_' + '{0:04}'.format(k) + '_' + '{0:04}'.format(p)
                            full_name = os.path.join(save_path, file_name + '.csv')   
                            with open((full_name), 'w', newline = '') as f: # writes the data to a .csv file
                                writer = csv.writer(f)
                                matrices = id_matrix_sim(n, l, m, g, k, p)
                                id_matrix = matrices['id_matrix']
                                c_matrix = matrices['c_matrix']
                                shape = id_matrix.shape # (rows, columns)
                                writer.writerow(['Iteration'])
                                for i in range(iterations):
                                    writer.writerow([i+1, 'ID%'])
                                    header = [i+1, '']
                                    for strain in range(n):
                                        header.append('strain' + str(strain+1))
                                    writer.writerow(header)
                                    for row in range(shape[0]):
                                        write_row = [i+1, 'strain' + str(row+1)]
                                        write_row.extend(id_matrix[row])
                                        writer.writerow(write_row)
                                    writer.writerow([i+1, 'c'])
                                    writer.writerow(header)
                                    for row in range(shape[0]):
                                        write_row = [i+1, 'strain' + str(row+1)]
                                        write_row.extend(c_matrix[row])
                                        writer.writerow(write_row)
                                    writer.writerow([])
                                    print(i)
                                                        
                                                                
                                                                

def run_c_q_sim(N, L, mu, generations, kappa, phi, iterations):
	current_path = os.path.dirname(os.path.abspath(__file__))
	for n in N:
		for l in L:
			for g in generations:
				for m in mu:
					for k in kappa:
						for p in phi:
							save_path = current_path + '/c_q Sim/n = ' + str(n) + '/L = ' + str(l)
							if not os.path.exists(save_path):
								os.makedirs(save_path)
							file_name = 'c_q_sim_' + '{0:04}'.format(n) + '_' + str(l) + '_' + '{0:04}'.format(m) + '_' + '{0:04}'.format(k) + '_' + '{0:04}'.format(p)
							full_name = os.path.join(save_path, file_name + '.csv')
							with open((full_name), 'w', newline = '') as f: # writes the data to a .csv file
								writer = csv.writer(f)
								header = ['Iteration']
								for q in range(2,n+1):
									header.append('c_' + str(q))
								header.extend(('n', 'L', 'mu', 'kappa', 'phi'))
								writer.writerow(header)
								for i in range(iterations):
									row = [(i+1)]
									row.extend(c_q_sim(n, l, m, g, k, p))
									row.extend((n,l,m,k,p))
									writer.writerow(row)
									print(i)
	
                
								
								
								
def run_mutation_sites_sim(L, mu, generations, kappa, phi, iterations, one_file):
	current_path = os.path.dirname(os.path.abspath(__file__))
	for l in L: # iterates over every length desired
		for m in mu: # iterates over every mutation rate desired
			for g in generations: # iterates over every number of generations desired
				for k in kappa: # iterates over every kappa desired
					for p in phi: # iterates over every phi desired
						# save_path = 'D:/BIGG DATA/Generation Sim/all GCs with kappa 0.5,1-6/L = 1000/GC% = ' + str(gc) + '/kappa = ' + str(k) + '/phi = ' + str(p)
						save_path = current_path + '/Mutation Sites Sim/L = ' + str(l) + '/kappa = ' + str(k) + '/phi = ' + str(p)
						for i in range(iterations): # runs it for the desired iterations
							if not os.path.exists(save_path):
								os.makedirs(save_path)
							file_name = 'mutation_sites_sim_' + str(l) + '_' + '{0:04}'.format(m) + '_' + str(g) + '_' + '{0:04}'.format(k) + '_' + '{0:04}'.format(p) + '_' + '{0:04}'.format(i+1)
							full_name = os.path.join(save_path, file_name + '.csv')   
							with open((full_name), 'w', newline = '') as f: # writes the data to a .csv file
								writer = csv.writer(f)
								writer.writerow(['Mutations', 'Mutation Sites', 'L', 'mu', 'kappa', 'phi']) # column headers
								mutations_per_generation = []

								data = [list(range(int(m*l),int((m*l*g)+1),int(m*l))), (mutation_sites_sim(l, m, g, k, p)), g*[l], g*[m], g*[k], g*[p]] # one variable is the number of SNPs and the other is the number of convergent mutations
								data = zip(*data) # formats the data so it can be written row by row and produce two long columns
								writer.writerows(data) # writes the data
							print(i)
						if(one_file):
							new_path = current_path + '/Mutation Sites Sim'
							get_mutation_sites_sim_averages(l, m, g, k, p, save_path, new_path)
