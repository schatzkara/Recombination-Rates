
import csv
import numpy as np
from process_genomes import read_in_strains
from process_genomes import id_matrix

path = '' # where the ID_matrix.csv is
filename = '' # of the .csv

def read_in_matrix(path,filename):
	with open((filename)) as d:
		reader = csv.reader(d)
		next(reader)
		next(reader)
		data = [r for r in reader] # this is a list of lists; the external list contains all the rows, the internal list contains each row element
		ID = np.empty([n,n], dtype = np.float, order = 'C')
		# c = np.empty([n,n], dtype = np.float, order = 'C')
		n = len(data)
		for strain1 in range(n):
			ID[strain1,strain1] = 1
			for strain2 in range(strain1+1,n):
				ID[strain1,strain2] = data[strain1][strain2+1]
				ID[strain2,strain1] = data[strain1][strain2+1]
	return ID

def trendline(identity):
	return c

def apply_trendline(ID):
	c_matrix = np.empty([n,n], dtype = np.float, order = 'C')
	for row in ID:
		for col in row:
			c = trendline(ID[row,col])
			c_matrix[row,col] = c
	return c_matrix

# path and filename or species
def calculate_c(**keyword_parameters):
	if 'path' and 'filename' in keyword_parameters:
		path = keyword_parameters['path']
		filename = keyword_parameters['filename']
		id_matrix = read_in_matrix(path,filename)
		c_matrix = apply_trendline(id_matrix)
		rows = np.sum(c_matrix, axis = 1)
		total = np.sum(rows, axis = 0)
		c = total[0]
		return c
	if 'species' in keyword_parameters:
		species = keyword_parameters[species]
		species = read_in_strains(species)
		id_matrix = id_matrix(species)
		c_matrix = apply_trendline(id_matrix)
		rows = np.sum(c_matrix, axis = 1)
		total = np.sum(rows, axis = 0)
		c = total[0]
	return c





