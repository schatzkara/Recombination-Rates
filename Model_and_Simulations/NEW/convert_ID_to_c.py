
import csv
import numpy as np
from process_genomes import read_in_strains
from process_genomes import id_matrix
from simulations import identity_sim 

path = '' # where the ID_matrix.csv is
filename = '' # name of the .csv

def read_in_matrix(path,filename):
	with open((filename)) as d:
		reader = csv.reader(d)
		next(reader)
		next(reader)
		data = [r for r in reader] # this is a list of lists; the external list contains all the rows, the internal list contains each row element
		n = len(data)
		ID = np.empty([n,n], dtype = np.float, order = 'C')
		# c = np.empty([n,n], dtype = np.float, order = 'C')
		for strain1 in range(n):
			ID[strain1,strain1] = 1
			for strain2 in range(strain1+1,n):
				ID[strain1,strain2] = data[strain1][strain2+1]
				ID[strain2,strain1] = data[strain1][strain2+1]
	return ID

def trendline(identity,regression_coefficients):
	a1 = regression_coefficients[0]
	a2 = regression_coefficients[1]
	a3 = regression_coefficients[2]
	c = a1*(identity**2) + a2*identity + a3
	return c

def calculate_trendline(L,mu,generations,GC_prop,kappa,phi,iterations):
	all_identities = generations*[None]
	all_cs = generations*[None]
	identites = generations*[None]
	cs = generations*[None]
	for g in range(generations):
		all_identities[g] = iterations*[None]
		all_cs[g] = iterations*[None]
	for i in range(iterations):
		values = identity_sim(L, mu, generations, GC_prop, kappa, phi)
		for g in range(generations):
			all_identities[g][i] = values[g][0]
			all_cs[g][i] = values[g][1]
	for g in range(generations):
		identites[g] = np.mean(all_identities[g])
		cs[g] = np.mean(all_cs[g])

	regression_coefficients = np.polyfit(identites, cs, 2)
	return regression_coefficients

def apply_trendline(ID):
	n = shape(ID)[0]
	regression_coefficients = calculate_trendline(L,mu,generations,GC_prop,kappa,phi,iterations)
	c_matrix = np.empty([n,n], dtype = np.float, order = 'C')
	for row in ID:
		for col in row:
			c = trendline(ID[row,col],regression_coefficients)
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





