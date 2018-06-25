#! python 3

# script to run the SNP simulation according to the users inputs
# time complexity for each L, GC, kappa combination: O(n^3)

# from SNPs_sim import sim_snps
from run_sims import run_SNP_sim
from decimal import Decimal
# import csv

# params:
# 	L (list of ints) = list of DNA strand lengths (units: base pairs)
# 	generations (list of ints) = list of the number of SNPs to introduce to each strain
# 	GC (list of floats) = list of the GC% composition of the ancestral strain
# 	kappa (list of floats) = list of the ratios of transitions to transversions
# 	phi (float) = probabiltiy of transversion to its complementary base pair
# 	iterations (int) = number of iterations to run
# 	one_file (boolean) = True if you just want one output file that contains the averages over all iterations; False if you want one file for each iteration
# return: a .csv file will be produced in the directory where this file is located with the data for a single interation (one .csv file for each iteration)

# enter parameters here
L = [1000]
generations = [300]
GC = [0.5]
kappa = [1.0]
phi = [1/2]
iterations = 1000
one_file = True

run_SNP_sim(L, generations, kappa, phi, iterations, one_file)

# for l in L: # iterates over every length desired
# 	for g in generations: # iterates over every number of SNPs desired
# 		for gc in GC: # iterates over every GC% desired
# 			for k in kappa: # iterates over every kappa desired
# 				for z in range(iterations): # runs it for the desired iterations
# 					with open(('SNPs_sim_data_' + str(l) + '_' + str(gc) + '_' + str(k) + '_' + '{0:04}'.format(z+1) + '.csv'), 'w', newline = '') as f: # writes the data to a .csv file
# 					    writer = csv.writer(f)
# 					    writer.writerow(['Mutations on Each Strain', 'CMs']) # column headers
# 					    data = [list(range(1,g+1)), (sim_snps(l,g,gc,k,phi))] # one variable is the number of SNPs and the other is the number of convergent mutations
# 					    data = zip(*data) # formats the data so it can be written row by row and produce two long columns
# 					    writer.writerows(data) # writes the data
# if(one_file):
# 		get_SNP_sim_averages(L, generations, GC, kappa, path)