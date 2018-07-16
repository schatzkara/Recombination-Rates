
# script to apply the model given m1 and m2 along the phylogeny and generate a matrix of c values

import numpy as np
from process_genomes import read_in_strains
from process_genomes import species_size
from process_genomes import genome_length
from process_genomes import theta_value
from completely_new_thing import get_min_m
from completely_new_thing import scale_newick_format_tree
from model import expected_c_given_ms

def apply_model_along_phylogeny(species_path, kappa, tree_string):
	# print('here')
	# ancestor = ''
	strains = read_in_strains(species_path)
	print(strains)
	strain_names = list(strains.keys())
	n = species_size(strains)
	L = genome_length(strains)
	# print('L = ' + str(L))
	# print(L)
	theta = theta_value(strains)
	mu = theta/(2*n)
	# print('mu = ' + str(mu))
	min_m = get_min_m(strains, L)
	# scaled_tree_string = tree_string
	scaled_tree_string = scale_newick_format_tree(strains, L, min_m, tree_string, 0)
	# print(strains)
	# SHARED = np.empty([n,n], dtype = np.float, order='C')
	# ANCESTRAL = np.empty([n,n], dtype = np.float, order='C')
	# RECOMBINANT = np.empty([n,n], dtype = np.float, order='C')
	CONVERGENT = np.empty([n,n], dtype = np.float, order='C')

	for s1 in range(n):
		strain1 = strains[strain_names[s1]]
		# SHARED[s1,s1] = L
		# RECOMBINANT[s1,s1] = 0
		# ANCESTRAL[s1,s1] = L
		CONVERGENT[s1,s1] = 0
		for s2 in range(s1,n):
			strain2 = strains[strain_names[s2]]
			# s,r,a,c = 0,0,0,0
			# for site in range(L):
			# 	if strain1[site] == strain2[site]:
			# 		s += 1
			# 		if strain1[site] == ancestor[site]:
			# 			a += 1
			s1_tree_location = scaled_tree_string.find(strain_names[s1][1:])
			# for x in range(10000):
			# 	print(s1_tree_location)
			s2_tree_location = scaled_tree_string.find(strain_names[s2][1:])
			# for x in range(10000):
			# 	print('DDDDDDDDDDDDDDDDDDDDDOOOOOOOOOOOOOOOOOOOOOOOOOOONNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE')

			start_length_1 = scaled_tree_string.find(':', s1_tree_location) + 1
			x1 = scaled_tree_string.find(',', start_length_1)
			y1 = scaled_tree_string.find(')', start_length_1)
			if x1 == -1:
				end_length_1 = y1
			elif y1 == -1:
				end_length_1 = x1
			else:
				end_length_1 = min(x1,y1)

			start_length_2 = scaled_tree_string.find(':', s2_tree_location) + 1
			x2 = scaled_tree_string.find(',', start_length_2)
			y2 = scaled_tree_string.find(')', start_length_2)
			if x2 == -1:
				end_length_2 = y2
			elif y2 == -1:
				end_length_2 = x2
			else:
				end_length_2 = min(x2,y2)

			length_1 = float(scaled_tree_string[start_length_1:end_length_1])
			length_2 = float(scaled_tree_string[start_length_2:end_length_2])
			# for x in range(10000):
			# 	print(length_1)
			# print(L)
			m_1 = int(length_1 * L + 1)
			m_2 = int(length_2 * L + 1)
			# for x in range(10000):
			# print('m_1 = ' + str(m_1))
			# print(mu)
			# print('m_2 = ' + str(m_2))
			generations_1 = int(m_1/mu)
			generations_2 = int(m_2/mu)
			# print('generations_1 = ' + str(generations_1))
			# print('generations_2 = ' + str(generations_2))

			c = expected_c_given_ms(L, m_1, m_2, mu, generations_1, generations_2, kappa, 0.5)

			# SHARED[s1,s2] = s
			# SHARED[s2,s1] = s
			# ANCESTRAL[s1,s2] = a
			# ANCESTRAL[s2,s1] = a
			CONVERGENT[s1,s2] = c
			CONVERGENT[s2,s1] = c
			# print('----------------------------HERE----------------------------')
			# RECOMBINANT[s1,s2] = s - a - c
			# RECOMBINANT[s2,s1] = s - a - c
	# print(strains)
	return {'strain_names': strain_names, 'Convergent': CONVERGENT}
	# return {'strain_names': strain_names, 'Shared': SHARED, 'Recombinant': RECOMBINANT, 'Ancestral': ANCESTRAL, 'Convergent': CONVERGENT}



