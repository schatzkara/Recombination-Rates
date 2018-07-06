
import csv
import numpy as np
from process_genomes import read_in_strains
from process_genomes import id_matrix
from simulations import identity_sim 
from model import expected_h_c_mg

# path = '' # where the ID_matrix.csv is
# filename = '' # name of the .csv

def read_in_matrix(path):
	with open((path)) as d:
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
			all_identities[g][i] = values[g+1][0]
			all_cs[g][i] = values[g+1][1]
	for g in range(generations):
		identites[g] = np.mean(all_identities[g])
		cs[g] = np.mean(all_cs[g])

	regression_coefficients = np.polyfit(identites, cs, 2)
	regression_coefficients2 = np.polyfit(identites, list(range(1,generations+1)), 2)
	print(regression_coefficients)
	print(regression_coefficients2)
	return {'regression_coefficients': regression_coefficients, 'regression_coefficients2': regression_coefficients2}

def apply_trendline(n, L, mu, generations, GC_prop, kappa, phi, iterations, ID, regression_coefficients):
	
	c_matrix = np.empty([n,n], dtype = np.float, order = 'C')
	for row in range(n):
		c_matrix[row,row] = 0
		for col in range(1,n):
			c = trendline(ID[row,col],regression_coefficients)
			c_matrix[row,col] = c
	return c_matrix

def get_generations(identity, regression_coefficients2):
	a1 = regression_coefficients2[0]
	a2 = regression_coefficients[1]
	a3 = regression_coefficients[2]
	generations = a1*(identity**2) + a2*identity + a3
	return generations

# path and filename or species
def calculate_c(L, mu, generations, GC_prop, kappa, phi, iterations, **keyword_parameters):
	if 'path' in keyword_parameters:
		path = keyword_parameters['path']
		# filename = keyword_parameters['filename']
		id_matrix = read_in_matrix(path)
		n = (id_matrix.shape)[0]
		print(id_matrix)

		regressions = calculate_trendline(L,mu,generations,GC_prop,kappa,phi,iterations)
		regression_coefficients = regressions['regression_coefficients']
		regression_coefficients2 = regressions['regression_coefficients2']

		c_matrix = apply_trendline(n, L, mu, generations, GC_prop, kappa, phi, iterations, id_matrix, regression_coefficients)
		print(c_matrix)
		rows = np.sum(c_matrix, axis = 1)
		total = np.sum(rows, axis = 0)
		c_sim = total

		total = 0
		for row in range(1,n):
			for col in range(row+1,n):
				total += id_matrix[row,col]
		average_id = total/((n*(n-1))/2)
		gen = get_generations(average_id, regression_coefficients2)
		c_model = expected_h_c_mg(n, L, mu, gen, kappa, phi)

	if 'species' in keyword_parameters:
		species = keyword_parameters[species]
		species = read_in_strains(species)
		id_matrix = id_matrix(species)
		c_matrix = apply_trendline(L, mu, generations, GC_prop, kappa, phi, id_matrix)
		rows = np.sum(c_matrix, axis = 1)
		total = np.sum(rows, axis = 0)
		c = total[0]
	return {'c_model': c_model, 'c_sim': c_sim}

def calculate_h_c(L, mu, generations, GC_prop, kappa, phi, iterations, **keyword_parameters):
	if 'path' in keyword_parameters:
		path = keyword_parameters['path']
		# filename = keyword_parameters['filename']
		# id_matrix = read_in_matrix(path)

		with open((path)) as d:
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
		# return ID
		n = (ID.shape)[0]
		# print(id_matrix)

		regressions = calculate_trendline(L,mu,generations,GC_prop,kappa,phi,iterations)
		regression_coefficients = regressions['regression_coefficients']
		regression_coefficients2 = regressions['regression_coefficients2']

		c_matrix = apply_trendline(n, L, mu, generations, GC_prop, kappa, phi, iterations, id_matrix, regression_coefficients)
		print(c_matrix)
		rows = np.sum(c_matrix, axis = 1)
		total = np.sum(rows, axis = 0)
		c_sim = total[0]

		total = 0
		for row in range(1,n):
			for col in range(row+1,n):
				total += id_matrix[row,col]
		average_id = total/((n*(n-1))/2)
		gen = get_generations(average_id, regression_coefficients2)
		c_model = expected_h_c_mg(n, L, mu, gen, kappa, phi)

	if 'species' in keyword_parameters:
		species = keyword_parameters[species]
		species = read_in_strains(species)
		id_matrix = id_matrix(species)
		c_matrix = apply_trendline(L, mu, generations, GC_prop, kappa, phi, id_matrix)
		rows = np.sum(c_matrix, axis = 1)
		total = np.sum(rows, axis = 0)
		c = total[0]
	return {'c_model': c_model, 'c_sim': c_sim}





