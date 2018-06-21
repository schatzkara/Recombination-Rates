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
				# 			# path = ('C:/Users/Owner/Documents/UNCG REU/Project/Recombination-Rates/Data/SNP Sim/L = 1000/GC% = 0.5/kappa = 1')   # path where the .csv files are located 
				# path = 'C:/Users/Owner/Documents/UNCG REU/Project/Recombination-Rates/Data'
				path = ('C:/Users/Owner/Documents/UNCG REU/Project/Recombination-Rates/Data/SNP Sim/L = ' + str(l) + '/GC% = ' + str(gc) + '/kappa = ' + str(int(k)))   # path where the .csv files are located 
				# values = 100*[10*[None]]
				# (values[0])[2] = 5
				# print(values[0])
				values2 = {}
				for x in range(g):
					values2[x] = []
				print(values2)

				# print(values)
							# datas = {}
							# values2 = {}
							# print(values)
							# print(values[0])
							# averages = []
							# stdevs = []
				# with open(('averages_test.csv'), 'w', newline = '') as f:
				with open(('averages_for_SNP_sim_' + str(l) + '_' + str(gc) + '_' + str(k) + '.csv'), 'w', newline = '') as f:
					writer = csv.writer(f)
					writer.writerow(['Mutations on Each Strain', 'Average CMs', 'STDev']) # column headers
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
								cms = int(data[row][1])
								# print(ms)
								# print(cms)
								# values[ms-1][j] = (cms)
								values2[ms-1].append(cms)
								# print(values[ms-1])
								# print(values2)
								# print(k)
								k += 1
						print('FILE NUMBER ' + str(j))
						j += 1
					# for x in range(10):
					# 	print(values[10])
					# 	print(len(values[10]))
					i = 1
					for item in values2.values():
						writer.writerow([i, np.mean(item), np.std(item)])
						i += 1
					# 		datas[k] = []
					# 		for row in data:
					# 			datas[k].append(row[1])
					# 		# print(data)
					# 		# for row in range(len(data)):
					# 		# 	# print(int(data[row][1]))
					# 		# 	values[row].append(int(data[row][1]))
					# 		# 	# print(values[row])
					# 	l = 0
					# 	for item in datas.values():
					# 		values2[l] = datas[l]
					# 		l += 1
					# 	k += 1

					# #print(values[0])
					# j = 1
					# for i in values:
					# 	print(i[0])
					# 	writer.writerow([j, np.average(i), np.std(i)])
					# 	j += 1
	# averages.append(np.mean(i))
	# stdevs.append(np.std(i))
