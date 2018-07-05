
from process_genomes import read_in_strains
from process_genomes import id_matrix
# from id_matrix import id_matrix
import numpy as np
import os
import glob
import csv

# n = 4

path = 'C:/Users/Owner/Documents/UNCG REU/Project/concatenates/clonal' # path where the .fa files are located 
for filename in glob.glob(os.path.join(path, '*.fa')): # finds the values for each species
	species = read_in_strains(filename)
	strains = list(species.keys())
	name = filename[len(path)+1:len(filename)-3]
	print(name)
	# number = len(species.keys())
	# name = filename.strip(path).strip('concat_') # strips off everything but the actual species name
	# print(name)
	# writer.writerow([name])
	# id_matrix = np.matrix([[1,2,3,4], [1,2,3,4], [1,2,3,4], [1,2,3,4]])
	identity_matrix = id_matrix(species)
	shape = identity_matrix.shape
	# header = ['']
	# header.extend(strains)
	# writer.writerow(header)
	# for row in range(shape[0]):
		# write_row = [strains[row]]
		# write_row.extend(identity_matrix[row])
		# writer.writerow(write_row)
	# writer.writerow([])
	# print(name)
	with open(('ID_Matrix_' + name + '.csv'), 'w', newline = '') as f: 
		writer = csv.writer(f)
		writer.writerow([name])
		# for filename in glob.glob(os.path.join(path, '*.fa')): # finds the values for each species
		# species = read_in_strains(filename)
		# strains = list(species.keys())
		# number = len(species.keys())
		# name = filename.strip(path).strip('concat_') # strips off everything but the actual species name
		# print(name)
		# writer.writerow([name])
		# id_matrix = np.matrix([[1,2,3,4], [1,2,3,4], [1,2,3,4], [1,2,3,4]])
		# identity_matrix = id_matrix(species)
		# shape = identity_matrix.shape
		header = ['']
		header.extend(strains)
		writer.writerow(header)
		for row in range(shape[0]):
			write_row = [strains[row]]
			write_row.extend(identity_matrix[row])
			writer.writerow(write_row)
		writer.writerow([])
	print(name)
		# writer.writerow(species.keys())
		# for x in range(number):
		# 	ids = (number+1)*[None]
		# 	ids[0] = species.keys()[x]
		# 	for y in range(number):
		# 		ids[y+1] = id_matrix[x,y]
		# 	writer.writerow(ids)


# writer = csv.writer(f)
# matrices = id_matrix_sim(n, l, m, g, k, p)
# id_matrix = matrices['id_matrix']
# c_matrix = matrices['c_matrix']
# shape = id_matrix.shape # (rows, columns)
# writer.writerow(['Iteration'])
# for i in range(iterations):
#         writer.writerow([i+1, 'ID%'])
#         header = [i+1, '']
#         for strain in range(n):
#                 header.append('strain' + str(strain+1))
#         writer.writerow(header)
#         for row in range(shape[0]):
#                 write_row = [i+1, 'strain' + str(row+1)]
#                 write_row.extend(id_matrix[row])
#                 writer.writerow(write_row)
#         writer.writerow([i+1, 'c'])
#         writer.writerow(header)
#         for row in range(shape[0]):
#                 write_row = [i+1, 'strain' + str(row+1)]
#                 write_row.extend(c_matrix[row])
#                 writer.writerow(write_row)
#         writer.writerow([])
#         print(i)