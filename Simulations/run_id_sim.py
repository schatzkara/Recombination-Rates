#! python 3

# script to run the SNP simulation according to the users inputs
# time complexity for each L,GC,kappa combination: O(n^3)

from SNPs_sim import sim_id_percent
from decimal import Decimal
import csv

# params:
# 	L (list of ints) = list of DNA strand lengths (units: base pairs)
# 	generations (list of ints) = list of the number of SNPs to introduce to each strain
# 	GC (list of floats) = list of the GC% composition of the ancestral strain
# 	kappa (list of floats) = list of the ratios of transitions to transversions
# 	phi (float) = probabiltiy of transversion to its complementary base pair
# 	iterations (int) = number of iterations to run
# return: a .csv file will be produced in the directory where this file is located with the data for a single interation (one .csv file for each iteration)

# enter parameters here
L = [1000]
generations = [300]
GC = [1.0]
kappa = [1.0, 2.0, 3.0]
phi = 1/2
iterations = 1000

for l in L: # iterates over every length desired
	for g in generations: # iterates over every number of SNPs desired
		for gc in GC: # iterates over every GC% desired
			for k in kappa: # iterates over every kappa desired
				for z in range(iterations): # runs it for the desired iterations
					with open(('ID_sim_data_' + str(l) + '_' + str(gc) + '_' + str(k) + '_' + '{0:04}'.format(z+1) + '.csv'), 'w', newline = '') as f: # writes the data to a .csv file
					    writer = csv.writer(f)
					    writer.writerow(['Mutations on Each Strain', 'ID%', 'CMs']) # column headers
					    sim = sim_id_percent(l,g,gc,k,phi)
					    for key in sim.keys():

					    # values = zip(*sim.values())
					    # data = [(sim.keys()), (values)] # one variable is the number of SNPs and the other is the number of convergent mutations
					    # data = zip(*data) # formats the data so it can be written row by row and produce two long columns
					    	writer.writerow([key+1, sim[key][0], sim[key][1]]) # writes the data
					    # writer.writerows(sim.items()) 