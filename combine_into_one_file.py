#! python 3

# script to get one file of averages for each L

# from SNPs_sim import sim_snps
# from decimal import Decimal
import csv
import os 
import glob
import numpy as np

# enter parameters here
L = [1000]
generations = [300]
GC = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
kappa = [1.0, 2.0, 3.0]
# phi = 1/2
# iterations = 1000


for l in L: # iterates over every length desired
	for g in generations:
# 		for gc in GC: # iterates over every GC% desired
# 			for k in kappa: # iterates over every kappa desired
		path = ('C:/Users/Owner/Documents/UNCG REU/Project/Recombination-Rates/Data/SNP Sim/Averages')   # path where the .csv files are located 
		values2 = {}
		for x in range(g):
			values2[x] = []
		with open(('all_averages_for_SNP_sim_' + str(l) + '.csv'), 'w', newline = '') as f:
			writer = csv.writer(f)
			writer.writerow(['Mutations on Each Strain', 'Average CMs', 'STDev', 'GC%', 'Kappa']) # column headers
			j = 0
			for filename in glob.glob(os.path.join(path, '*.csv')):
				gc = filename[-11:-8]
				k = filename[-7:-4]
				with open(filename) as d:
					reader = csv.reader(d)
					next(reader) # skips the column headers
					for row in reader:
						writer.writerow([row[0], row[1], row[2], gc, k])
				print('FILE NUMBER ' + str(j))
				j += 1

for l in L: # iterates over every length desired
	for g in generations:
# 		for gc in GC: # iterates over every GC% desired
# 			for k in kappa: # iterates over every kappa desired
		path = ('C:/Users/Owner/Documents/UNCG REU/Project/Recombination-Rates/Data/ID Sim/Averages')   # path where the .csv files are located 
		values2 = {}
		for x in range(g):
			values2[x] = []
		with open(('all_averages_for_ID_sim_' + str(l) + '.csv'), 'w', newline = '') as f:
			writer = csv.writer(f)
			writer.writerow(['Mutations on Each Strain', 'Average ID%', 'ID% STDev', 'Average CMs', 'CMs STDev', 'GC%', 'Kappa']) # column headers
			j = 0
			for filename in glob.glob(os.path.join(path, '*.csv')):
				gc = filename[-11:-8]
				k = filename[-7:-4]
				with open(filename) as d:
					reader = csv.reader(d)
					next(reader) # skips the column headers
					for row in reader:
						writer.writerow([row[0], row[1], row[2], row[3], row[4], gc, k])
				print('FILE NUMBER ' + str(j))
				j += 

for l in L: # iterates over every length desired
	for g in generations:
# 		for gc in GC: # iterates over every GC% desired
# 			for k in kappa: # iterates over every kappa desired
		path = ('C:/Users/Owner/Documents/UNCG REU/Project/Recombination-Rates/Data/SNP Sim/Averages')   # path where the .csv files are located 
		values2 = {}
		for x in range(g):
			values2[x] = []
		with open(('all_model_id_vs_cm_data_' + str(l) + '.csv'), 'w', newline = '') as f:
			writer = csv.writer(f)
			writer.writerow(['Mutations on Each Strain', 'Expected ID%', 'Expected CMs', 'Kappa']) # column headers
			j = 0
			for filename in glob.glob(os.path.join(path, '*.csv')):
				# gc = filename[-11:-8]
				k = filename[-7:-4]
				with open(filename) as d:
					reader = csv.reader(d)
					next(reader) # skips the column headers
					for row in reader:
						writer.writerow([row[0], row[1], row[2], k])
				print('FILE NUMBER ' + str(j))
				j += 1