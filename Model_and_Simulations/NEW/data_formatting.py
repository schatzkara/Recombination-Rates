#! python 3

# this script contains several functions that take existing data and format it to making processing in R easier 

import csv
import os 
import glob
import numpy as np

def combine_mutation_sim_data(path):
	with open(('all_mutation_sim_data.csv'), 'w', newline = '') as f:
		writer = csv.writer(f)
		# writer.writerow(['Iteration Number', 'Average CMs', 'L', 'mu', 'kappa', 'phi', 'Expected CMs'])
		j = 1
		for filename in glob.glob(os.path.join(path, '*.csv')):
			with open(filename) as d:
				reader = csv.reader(d)
				if(j > 1):
					next(reader)
				for row in reader:
					writer.writerow(row)
			print('file number ' + str(j))
			j += 1

# this function takes all the .csv files that contain the averages of all the iterations of the SNP sim (for a certain L, GC%, and kappa) and puts them into one file with columns noting the GC% and kappa values
# this makes it easier to process the data in R beacuse it puts it in tidy format
# params:
# 	L (list of ints) = genome lengths that the sim ran for 
# 	generations (list of ints) = generations that the sim ran for 
# 	path (string) = the path where all the .csv files currently are (make sure that it only contains the files you want to combine)
# output: a .csv file containing all the data with the following columns: Mutations on Each Strain, Average CMs%, STDev, GC%, Kappa
def combine_SNP_sim_data(L, generations, path):
	# for l in L: # iterates over every length the sim ran for 
	# 	for g in generations: # iterates over every number of generations the sim ran for 
			# for gc in GC: # iterates over every GC% desired
				# for k in kappa: # iterates over every kappa desired
			# path = ('C:/Users/Owner/Documents/UNCG REU/Project/Recombination-Rates/Data/SNP Sim/Averages')   # path where the .csv files are located 
			# values2 = {}
			# for x in range(g):
				# values2[x] = []
	with open(('all_averages_for_SNP_sim_data.csv'), 'w', newline = '') as f: # opens the file that the data will be written to
		writer = csv.writer(f)
		# writer.writerow(['Mutations on Each Strain', 'Average CMs', 'STDev', 'L', 'GC%', 'kappa', 'phi']) # column headers
		j = 1 # counter for the number of files
		for filename in glob.glob(os.path.join(path, '*.csv')): # iterates over every .csv file in the path
			# gc = filename[-11:-8] # extracts the GC% from the file name
			# k = filename[-7:-4] # extracts the kappa value from the file name
			with open(filename) as d: # reads the file
				reader = csv.reader(d)
				if(j > 1):
					next(reader) # skips the column headers
				for row in reader: # writes each row to the new file
					writer.writerow(row)
					# writer.writerow([row[0], row[1], row[2], gc, k])
			print('FILE NUMBER ' + str(j))
			j += 1

# this function takes all the .csv files that contain the averages of all the iterations of the ID sim (for a certain L, GC%, and kappa) and puts them into one file with columns noting the GC% and kappa values
# this makes it easier to process the data in R beacuse it puts it in tidy format
# params:
# 	L (list of ints) = genome lengths that the sim ran for 
# 	generations (list of ints) = generations that the sim ran for 
# 	path (string) = the path where all the .csv files currently are (make sure that it only contains the files you want to combine)
# output: a .csv file containing all the data with the following columns: Mutations on Each Strain, Average ID%, ID% STDev, Average CMs, CMs STDev, GC%, Kappa
def combine_ID_sim_data(L, generations, path):
	# for l in L: # iterates over every length the sim ran for 
	# 	for g in generations: # iterates over every number of generations the sim ran for 
			# for gc in GC: # iterates over every GC% desired
				# for k in kappa: # iterates over every kappa desired
			# path = ('C:/Users/Owner/Documents/UNCG REU/Project/Recombination-Rates/Data/ID Sim/Averages')   # path where the .csv files are located 
			# values2 = {}
			# for x in range(g):
			# 	values2[x] = []
	with open(('all_averages_for_ID_sim_data.csv'), 'w', newline = '') as f: # opens the file that the data will be written to
		writer = csv.writer(f)
		# writer.writerow(['Mutations on Each Strain', 'Average ID%', 'ID% STDev', 'Average CMs', 'CMs STDev', 'L', 'GC%', 'kappa', 'phi']) # column headers
		j = 1 # counter for the number of files 
		for filename in glob.glob(os.path.join(path, '*.csv')): # iterates over every .csv file in the path
			# gc = filename[-11:-8] # extracts the GC% from the file name 
			# k = filename[-7:-4] # extracts the kappa value from the file name
			with open(filename) as d: # reads the file
				reader = csv.reader(d)
				if(j > 1):
					next(reader) # skips the column headers
				for row in reader: # writes each row to the new file
					writer.writerow(row)
					# writer.writerow([row[0], row[1], row[2], row[3], row[4], gc, k])
			print('FILE NUMBER ' + str(j))
			j += 1

# this function takes all the .csv files that contain the model cm vs ID data (for a certain L and kappa) and puts them into one file with columns noting the and kappa values
# this makes it easier to process the data in R beacuse it puts it in tidy format
# params:
# 	L (list of ints) = genome lengths that the model ran for 
# 	generations (list of ints) = generations that the model ran for 
# 	path (string) = the path where all the .csv files currently are (make sure that it only contains the files you want to combine)
# output: a .csv file containing all the data with the following columns: Mutations on Each Strain, Expected ID%, Expected CMs, Kappa
def combine_model_ID_data(L, generations, path):
	for l in L: # iterates over every length the model ran for
		for g in generations: # iterates over every number of generations the model ran for 
			# for gc in GC: # iterates over every GC% desired
				# for k in kappa: # iterates over every kappa desired
			# path = ('C:/Users/Owner/Documents/UNCG REU/Project/Recombination-Rates/Data/Model CM vs ID')   # path where the .csv files are located 
			# values2 = {}
			# for x in range(g):
			# 	values2[x] = []
			with open(('all_model_id_vs_cm_data_' + str(l) + '_' + str(g) + '.csv'), 'w', newline = '') as f: # open the file that the data will be written to
				writer = csv.writer(f)
				writer.writerow(['Mutations on Each Strain', 'Expected ID%', 'Expected CMs', 'Kappa']) # column headers
				j = 1 # counter for the number of files
				for filename in glob.glob(os.path.join(path, '*.csv')): # iterates over every .csv file in the path
					# gc = filename[-11:-8]
					k = filename[-7:-4] # extracts the kappa values from the file name
					with open(filename) as d: # reads the file
						reader = csv.reader(d)
						next(reader) # skips the column headers
						for row in reader: # writes each row to the new file
							writer.writerow([row[0], row[1], row[2], k])
					print('FILE NUMBER ' + str(j))
					j += 1