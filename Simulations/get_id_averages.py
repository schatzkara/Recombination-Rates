#! python 3

# script to take in 1000 .csv files and output the average and standard deviations of each row

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
		for gc in GC: # iterates over every GC% desired
			for k in kappa: # iterates over every kappa desired
				path = ('C:/Users/Owner/Documents/UNCG REU/Project/Recombination-Rates/Data/ID Sim/L = ' + str(l) + '/GC% = ' + str(gc) + '/kappa = ' + str(int(k)))   # path where the .csv files are located 
				values = {}
				values2 = {}
				for x in range(g):
					values[x] = []
					values2[x] = []
				# print(values2)

				with open(('averages_for_id_sim_' + str(l) + '_' + str(gc) + '_' + str(k) + '.csv'), 'w', newline = '') as f:
					writer = csv.writer(f)
					writer.writerow(['Mutations on Each Strain', 'Average ID%', 'ID% STDev', 'Average CMs', 'CMs STDev']) # column headers
					j = 0
					for filename in glob.glob(os.path.join(path, '*.csv')):
						k = 0
						with open(filename) as d:
							reader = csv.reader(d)
							next(reader) # skips the column headers
							data = [r for r in reader] # this is a list of lists; the external list contains all the rows, the internal list contains each row element
							# print(len(data))
							for row in range(len(data)): 
								ms = int(data[row][0])
								idp = float(data[row][1])
								cms = int(data[row][2])
								values[ms-1].append(idp)
								values2[ms-1].append(cms)
								k += 1
						print('FILE NUMBER ' + str(j))
						j += 1
					i = 1
					for x in range(g):
						writer.writerow([(x+1), np.mean(values[x]), np.std(values[x]), np.mean(values2[x]), np.std(values2[x])])
					# for item in values2.values():
					# 	writer.writerow([i, np.mean(item), np.std(item)])
					# 	i += 1