#! python 3

# script that contains several functions that take in 1000 .csv files (that correspond to the 1000 iterations of a single sim) and output the average and standard deviations of each row (which corresponds to a generation)

import csv
import os 
import glob
import numpy as np

# this function gets the average and standard deviation of each generation for all iterations of the SNP sim
# params:
# 	L (list of ints) = the genome lengths that the sim ran for 
# 	generations (list of ints) = the number of generations that the sim ran for
# 	GC (list of floats) = the GC%s that the sim ran for
# 	kappa (list of floats) = the kappa values that the sim ran for 
# 	phi (list of floats) = the phi values that the sim ran for
# 	path (string) = the path where all the .csv files are located (make sure that only the files you want the averages for are in the path)
# output: a .csv file containing the CM average and standard deviation of each generation with the following columns: Mutations on Each Strain, Average CMs, STDev
def get_SNP_sim_averages(L, generations, GC, kappa, phi, data_path, save_path):
	for l in L: # iterates over every length the sim ran for
		for g in generations: # iterates over every number of generations the sim ran for
			for gc in GC: # iterates over every GC% the sim ran for
				for k in kappa: # iterates over every kappa the sim ran for
					for p in phi: # iterates over every phi the sim ran for
						# path = ('C:/Users/Owner/Documents/UNCG REU/Project/Recombination-Rates/Data/SNP Sim/L = ' + str(l) + '/GC% = ' + str(gc) + '/kappa = ' + str(int(k)))   # path where the .csv files are located 
						data = {} # a dictionary to hold all the data from all the files; has key: generation number - 1 and value: list of CMs for that generation
						for x in range(g):
							data[x] = []
						file_name = 'averages_for_SNP_sim_' + str(l) + '_' + str(g) + '_' + str(gc) + '_' + str(k) + '_' + str(p)
						full_name = os.path.join(save_path, file_name + '.csv')
						with open((full_name), 'w', newline = '') as f: # opens the file that the data will be written to
							writer = csv.writer(f)
							writer.writerow(['Mutations on Each Strain', 'Average CMs', 'STDev', 'L', 'GC%', 'kappa', 'phi']) # column headers
							j = 1 # counter for the number of files
							for filename in glob.glob(os.path.join(data_path, '*.csv')): # iterates over every .csv file in the path
								with open(filename) as d: # reads the file
									reader = csv.reader(d)
									next(reader) # skips the column headers
									file_data = [r for r in reader] # this is a list of lists; the external list contains all the rows, the internal list contains each row element
									for row in range(len(file_data)): # iterates over every row of the data; the row corresponds to a generation
										ms = int(file_data[row][0]) # extracts the number of mutations on each strand from the data
										cms = int(file_data[row][1]) # extracts the number of convergent mutations from the data
										data[ms-1].append(cms)
								print('FILE NUMBER ' + str(j))
								j += 1
							i = 1 # counter for the generation number 
							for item in data.values(): # writes each row to the new file; each row corresponds to a generation
								writer.writerow([i, np.mean(item), np.std(item), l, gc, k, p])
								i += 1

# this function gets the average and standard deviation of each generation for all iterations of the ID sim
# params:
# 	L (list of ints) = the genome lengths that the sim ran for 
# 	generations (list of ints) = the number of generations that the sim ran for
# 	GC (list of floats) = the GC%s that the sim ran for
# 	kappa (list of floats) = the kappa vales that the sim ran for 
# 	phi (list of floats) = the phi values that the sim ran for
# 	path (string) = the path where all the .csv files are located (make sure that only the files you want the averages for are in the path)
# output: a .csv file containing the ID% and CM averages and standard deviations of each generation with the following columns: Mutations on Each Strain, Average ID%, ID% STDev, Average CMs, CMs STDev
def get_ID_sim_averages(L, generations, GC, kappa, phi, data_path, save_path):
	for l in L: # iterates over every length the sim ran for
		for g in generations: # iterates over every number of generations the sim ran for
			for gc in GC: # iterates over every GC% the sim ran for
				for k in kappa: # iterates over every kappa the sim ran for
					for p in phi: # iterates over every phi the sim ran for
						# path = ('C:/Users/Owner/Documents/UNCG REU/Project/Recombination-Rates/Data/ID Sim/L = ' + str(l) + '/GC% = ' + str(gc) + '/kappa = ' + str(int(k)))   # path where the .csv files are located 
						values = {} # a dictionary to hold all the data for the ID% from all the files; has key: generation number - 1 and value: list of ID%s for that generation
						values2 = {} # a dictionary to hold all the data for the CMs from all the files; has key: generation number - 1 and value: list of CMs for that generation
						for x in range(g):
							values[x] = []
							values2[x] = []
						file_name = 'averages_for_id_sim_' + str(l) + '_' + str(g) + '_' + str(gc) + '_' + str(k) + '_' + str(p)
						full_name = os.path.join(save_path, file_name + '.csv')
						with open((full_name), 'w', newline = '') as f: # opens the file that the data will be written to
							writer = csv.writer(f)
							writer.writerow(['Mutations on Each Strain', 'Average ID%', 'ID% STDev', 'Average CMs', 'CMs STDev', 'L', 'GC%', 'kappa', 'phi']) # column headers
							j = 1 # counter for the number of files
							for filename in glob.glob(os.path.join(data_path, '*.csv')): # iterates over every .csv file in the path
								with open(filename) as d: # reads the file
									reader = csv.reader(d)
									next(reader) # skips the column headers
									data = [r for r in reader] # this is a list of lists; the external list contains all the rows, the internal list contains each row element
									for row in range(len(data)): # iterates over every row of the data; the row corresponds to a generation
										ms = int(data[row][0]) # extracts the number of mutations on each strand from the data
										idp = float(data[row][1]) # extracts the ID% from the data
										cms = int(data[row][2]) # extracts the number of convergent mutations from the data
										values[ms-1].append(idp) 
										values2[ms-1].append(cms)
								print('FILE NUMBER ' + str(j))
								j += 1
							for x in range(g): # writes each row to the new file; each row corresponds to a generation
								writer.writerow([(x+1), np.mean(values[x]), np.std(values[x]), np.mean(values2[x]), np.std(values2[x])], l, gc, k, p)