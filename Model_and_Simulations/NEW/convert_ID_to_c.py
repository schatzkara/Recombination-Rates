
import csv
import numpy as np
import math
from process_genomes import read_in_strains
from process_genomes import id_matrix
from simulations import identity_sim 
from simulations import c_q_sim
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
	data = [identites, cs]
	data = zip(*data)
	with open('data.csv', 'w', newline = '') as f:
		writer = csv.writer(f)
		writer.writerows(data)

	# print(identites)
	# print(cs)

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
	a2 = regression_coefficients2[1]
	a3 = regression_coefficients2[2]
	generations = int(a1*(identity**2) + a2*identity + a3)
	return generations

# path and filename or species
def calculate_c(real_L, L, mu, generations, GC_prop, kappa, phi, iterations, **keyword_parameters):
	if 'path' in keyword_parameters:
		path = keyword_parameters['path']
		# filename = keyword_parameters['filename']
		id_matrix = read_in_matrix(path)
		n = 2 # (id_matrix.shape)[0]
		# print(id_matrix)

		regressions = calculate_trendline(L,mu,generations,GC_prop,kappa,phi,iterations)
		regression_coefficients = regressions['regression_coefficients']
		regression_coefficients2 = regressions['regression_coefficients2']

		c_matrix = apply_trendline(n, L, mu, generations, GC_prop, kappa, phi, iterations, id_matrix, regression_coefficients)
		print(c_matrix)
		rows = np.sum(c_matrix, axis = 1)
		total = np.sum(rows, axis = 0)
		c_sim = total * (real_L/L)

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

def calculate_h_c(real_L, L, mu, generations, GC_prop, kappa, phi, iterations, **keyword_parameters):
	if 'path' in keyword_parameters:
		path = keyword_parameters['path']
		# filename = keyword_parameters['filename']
		# id_matrix = read_in_matrix(path)
		n = 0
		with open((path)) as d:
			reader = csv.reader(d)
			next(reader)
			next(reader)
			data = [r for r in reader] # this is a list of lists; the external list contains all the rows, the internal list contains each row element
			n = len(data)
			id_matrix = np.empty([n,n], dtype = np.float, order = 'C')
			# c = np.empty([n,n], dtype = np.float, order = 'C')
			for strain1 in range(n):
				id_matrix[strain1,strain1] = 1
				for strain2 in range(strain1+1,n):
					id_matrix[strain1,strain2] = data[strain1][strain2+1]
					id_matrix[strain2,strain1] = data[strain1][strain2+1]
		# return ID
		# print(id_matrix)

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
		data = [list(range(1,generations+1)), identites, cs]
		data = zip(*data)
		with open('data.csv', 'w', newline = '') as f:
			writer = csv.writer(f)
			writer.writerow(['Generation', 'ID%', 'c'])
			writer.writerows(data)
		regression_coefficients = np.polyfit(identites, cs, 2)
		regression_coefficients2 = np.polyfit(identites, list(range(1,generations+1)), 2)
		print(regression_coefficients)
		print(regression_coefficients2)
		# regressions = calculate_trendline(L,mu,generations,GC_prop,kappa,phi,iterations)
		# regression_coefficients = regressions['regression_coefficients']
		a0 = regression_coefficients[0]
		a1 = regression_coefficients[1]
		a2 = regression_coefficients[2]
		# regression_coefficients2 = regressions['regression_coefficients2']
		b0 = regression_coefficients2[0]
		b1 = regression_coefficients2[1]
		b2 = regression_coefficients2[2]

		c_matrix = np.empty([n,n], dtype = np.float, order = 'C')
		for row in range(n):
			c_matrix[row,row] = 0
			for col in range(row+1,n):
				identity = id_matrix[row,col]
				c = a0*(identity**2) + a1*identity + a2
				c_matrix[row,col] = c
				c_matrix[col,row] = c
		# print(c_matrix)
		# c_matrix = apply_trendline(n, L, mu, generations, GC_prop, kappa, phi, iterations, id_matrix, regression_coefficients)
		# print(c_matrix)
		# rows = np.sum(c_matrix, axis = 1)
		# total = np.sum(rows, axis = 0)
		c_sim = ((c_matrix.sum())/2)*(real_L/L)


		total = 0
		num = 0
		for row in range(1,n):
			for col in range(row+1,n):
				num += 1
				total += id_matrix[row,col]
		average_id = total/num
		gen = int(b0*(average_id**2) + b1*average_id + b2)
		all_cqs = (n-2)*[None]
		mean_cqs = (n-2)*[None]
		for x in range(n-2):
			all_cqs[x] = iterations*[None]
		for i in range(iterations):
			c_qs = c_q_sim(n, real_L, mu, gen, kappa, phi)
			for x in range(n-2):
				all_cqs[x][i] = c_qs[x]
		for x in range(n-2):
			mean_cqs[x] = np.mean(all_cqs[x])
		# print(all_cqs)
		# print(mean_cqs)

		correction_factor = 0
		for q in range(3,n):
			c_q = mean_cqs[q-2]
			overcounts = ((q*(q-1))/2) - 1
			correction_factor += c_q * overcounts * (real_L/L)

		c_sim_corrected = c_sim - correction_factor


		# total = 0
		# num = 0
		# for row in range(1,n):
		# 	for col in range(row+1,n):
		# 		num += 1
		# 		total += id_matrix[row,col]
		# average_id = total/num
		# gen = int(b0*(average_id**2) + b1*average_id + b2)
		# gen = get_generations(average_id, regression_coefficients2)
		c_model = (expected_h_c_mg(n, L, mu, gen, kappa, phi)) * (real_L/L)

		c_model_corrected = c_model - correction_factor


		print(c_model)

	# if 'species' in keyword_parameters:
	# 	species = keyword_parameters[species]
	# 	species = read_in_strains(species)
	# 	id_matrix = id_matrix(species)
	# 	c_matrix = apply_trendline(L, mu, generations, GC_prop, kappa, phi, id_matrix)
	# 	rows = np.sum(c_matrix, axis = 1)
	# 	total = np.sum(rows, axis = 0)
	# 	c = total[0]
	# return {'c_model': c_model, 'c_sim': c_sim_corrected}

	print({'c_model': c_model, 'c_model_corrected': c_model_corrected, 'c_sim': c_sim, 'c_sim_corrected': c_sim_corrected})





