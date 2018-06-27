
from process_genomes import read_in_strains
from IDMatrix import id_matrix
import numpy as np
import os
import glob
import csv

# n = 4

path = 'C:/Users/Owner/Documents/UNCG REU/Project/Recombination-Rates/concatenates' # path where the .fa files are located 
with open(('species_IDs.csv'), 'w', newline = '') as f: 
	writer = csv.writer(f)
	# writer.writerow(['species', 'pi', 'theta', 'GC%']) # column headers
	for filename in glob.glob(os.path.join(path, '*.fa')): # finds the values for each species
		species = read_in_strains(filename)
		number = len(species.keys())
		name = filename.strip('C:/Users/Owner/Documents/UNCG REU/Project/Recombination-Rates/concatenates/concat_') # strips off everything but the actual species name
		print(name)
		writer.writerow([name])
		# id_matrix = np.matrix([[1,2,3,4], [1,2,3,4], [1,2,3,4], [1,2,3,4]])
		id_matrix = id_matrix(species)
		writer.writerow(species.keys())
		for x in range(number):
			ids = (number+1)*[None]
			ids[0] = species.keys()[x]
			for y in range(number):
				ids[y+1] = id_matrix[x,y]
			writer.writerow(ids)


		# writer.writerow([name, pi, theta, GC_average])