
# -*- coding: utf-8 -*-

import subprocess
import pyvolve
import numpy as np
import csv
import os
from process_genomes import read_in_strains
from process_genomes import read_in_reduced_strains
from get_internal_nodes import get_internal_nodes
from process_genomes import species_size
from process_genomes import genome_length
from process_genomes import pi_value
from process_genomes import theta_value
from completely_new_thing import get_min_m
from completely_new_thing import get_max_m
from model_application import find_parents
from model_application import get_branch_lengths
from model_application import find_MRCA
from model_application import get_distances_to_MRCA
# from completely_new_thing import scale_newick_format_tree
from completely_new_thing import scale_branch_lengths
from simulations import generate_ancestor
from completely_new_thing import scale 
from model import expected_c_given_ms

# species_alignment = '/mnt/c/Users/Owner/Documents/UNCG/Project/BIGG_DATA/Useful_Data/Concatenates,Trees,Homoplasies/Aayyy_Clonal/Bacillus_anthracis/concat_universal.fa'
reduced_species_alignment = '/mnt/c/Users/Owner/Documents/UNCG/Project/standard-RAxML/Bacillus_subtilis/concat_universal.fa.reduced'
raxml_path = '/mnt/c/Users/Owner/Documents/UNCG/Project/standard-RAxML/Bacillus_subtilis'
tree_file = 'RAxML_bestTree.tree'
rooted_tree_file = 'RAxML_rootedTree.root'
# ancestral_alignment = 'RAxML_marginalAncestralStates.anc'
ancestral_tree_file = 'RAxML_nodeLabelledRootedTree.anc'
# kappa_file = '/mnt/c/Users/Owner/Documents/UNCG/Project/BIGG_DATA/Useful_Data/Concatenates,Trees,Homoplasies/Aayyy_Clonal/Bacillus_anthracis/kappa.txt'

def get_SCAR_matrices(species_alignment, ancestral_alignment, kappa_file, mu):
	# get_tree_string(species_alignment, raxml_path)
	reduced = os.path.exists(reduced_species_alignment)
	# get_tree_root(tree_file, raxml_path)
	# get_ancestors(rooted_tree_file, species_alignment, raxml_path, reduced)

	if not reduced:
		strains = read_in_strains(species_alignment)
	else:
		strains = read_in_reduced_strains(reduced_species_alignment) # dictionary with the genomes of all the strains; key = strain name, value = genome
	# for strain in strains.keys():
	# 	print(strain)
	# 	print(strains[strain][:10])

	L = genome_length(strains) # number of base pairs in the genome
	n = species_size(strains) # number of extant strains
	strain_names = list(strains.keys()) # list of all the extant strain names

	SHARED = np.empty([n,n], dtype = np.float, order='C') # a matrix of the number of nucleotides shared between two strains; the (i,j) entry is the number of nucleotides that are the same between strain i and strain j
	CONVERGENT = np.empty([n,n], dtype = np.float, order='C') # a matrix of the number of nucleotides that match due to convergent mutation between two strains; the (i,j) entry is the number of convergent mutations between strain i and strain j
	ANCESTRAL = np.empty([n,n], dtype = np.float, order='C') # a matrix of the number of nucleotides that match due to direct inheritence from the ancestor; the (i,j) entry is the number of nucleotides that were inherited by both strain i and strain j
	RECOMBINANT = np.empty([n,n], dtype = np.float, order='C') # a matrix of the number of nucleotides that match due to a recombination event; the (i,j) entry is the number of nucleotides that were recombined between strain i and strain j
	RATES = np.empty([n,n], dtype = np.float, order='C')

	tree_file = open((os.path.join(raxml_path, rooted_tree_file)), 'r')
	rooted_tree_string = list(tree_file)[0]
	tree_file = open((os.path.join(raxml_path, ancestral_tree_file)), 'r')
	ancestral_tree_string = list(tree_file)[0]
	# strains = read_in_strains(species_alignment) # dictionary with the genomes of all the strains; key = strain name, value = genome
	internal_nodes = get_internal_nodes(os.path.join(raxml_path, ancestral_alignment))
	# all_nodes = internal_nodes

	internal_nodes,ancestral_tree_string = rename_ancestors(internal_nodes, strain_names, ancestral_tree_string)

	all_nodes = {}
	for key in strains.keys():
		all_nodes[key] = strains[key]
	for key in internal_nodes.keys():
		all_nodes[key] = internal_nodes[key]

	# strain_names = list(strains.keys()) # list of all the extant strain names
	all_node_names = list(all_nodes.keys())
	# print(strain_names)
	# print(all_node_names)



	# n = species_size(strains) # number of extant strains
	total_pairs = int((n*(n-1)) / 2) # the total number of strain pairs that will be compared
	# L = genome_length(strains) # number of base pairs in the genome
	pi = pi_value(strains)
	theta = theta_value(strains) # proportion of the genome that is polymorphic
	# print(theta)
	# print(n)
	# mu = (theta)/(2*n) # mutation rate in mutations per base pair per generation
	# print(mu)

	# tree_file = open((os.path.join(raxml_path, rooted_tree_file)), 'r')
	# rooted_tree_string = list(tree_file)[0]
	# tree_file = open((os.path.join(raxml_path, ancestral_tree_file)), 'r')
	# ancestral_tree_string = list(tree_file)[0]

	# print(rooted_tree_string)
	# print(ancestral_tree_string)
	complete_tree_string = merge_trees(rooted_tree_string, ancestral_tree_string)
	# print(complete_tree_string)

	kappa_file = open(kappa_file, 'r')
	kappa = float(list(kappa_file)[0])
	min_m = get_min_m(strains, L) # minimum number of mutations that could account for all the polymorphisms in the species
	max_m = get_max_m(strains, L, complete_tree_string)
	
	###############################################################################
	##### CHANGE THIS!!!!!!!!!!!!!!!!!!!!! ########################################
	###############################################################################
	# scaled_tree_string = complete_tree_string

	scaled_tree_string = scale_branch_lengths(L, complete_tree_string, min_m, max_m, pi, theta, kappa, 1) # scale_newick_format_tree(strains, L, min_m, tree_string, 0) # the tree_string scaled by min_m
		# L, tree_string, min_m, max_m, real_pi, real_theta, kappa
	phylogeny = pyvolve.read_tree(tree = scaled_tree_string)
	# pyvolve.print_tree(phylogeny)

	g = open('scaled_tree.txt', 'w')
	g.write(scaled_tree_string)
	g.close()

	# updated_tree_info = name_nodes(tree_string, strain_names) 
	# scaled_tree_string = updated_tree_info['tree_string'] # version of the tree_string where every node is labeled
	# new_nodes = updated_tree_info['new_nodes'] # the new node names that were added
	# all_nodes = all_node_names + new_nodes # list of all the node names in the pyhlogenetic tree
	# print(all_nodes) 

	# parents = find_parents(strain_names, tree_string) # a dictionary of the sequence of parents of each strain; key = strain name, value = list of the parents in order of increasing distance from the strain
	parents = find_parents(all_node_names, scaled_tree_string)
	print('found parents')
	# print(parents)
	distances = get_branch_lengths(all_node_names, scaled_tree_string) # a dictionary of the distances of each strain to its closest ancestor; key = strain name, value = distance to its closest ancestor
	print('found distances')
	# print(distances)

	count = 1 # a counter for the current strain pair number that is being processed
	total = 0
	for s1 in range(n): # allows each strain to be strain 1
		strain1 = strain_names[s1]
		genome1 = strains[strain1]
		SHARED[s1,s1] = L
		CONVERGENT[s1,s1] = 0 # there can be no convergent mutations between a strain and itself
		ANCESTRAL[s1,s1] = L
		RECOMBINANT[s1,s1] = 0
		for s2 in range(s1+1,n): # allows each strain after strain 1 to be strain 2
			strain2 = strain_names[s2]
			genome2 = strains[strain2]

			MRCA = find_MRCA(strain1, strain2, parents) # the Most Recent Common Ancestor between the two strains
			MRCA_genome = all_nodes[MRCA]

			s,a = get_s_a(genome1, genome2, MRCA_genome, L)
			print('got s and a')

			c,pair_distances = get_c(strain1, strain2, MRCA, parents, scaled_tree_string, distances, L, mu, kappa)
			print('got c')

			r = s - c - a

			# fills in the appropriate values to the S,C,A,R matrices for the current strain pair
			SHARED[s1,s2] = s
			SHARED[s2,s1] = s
			CONVERGENT[s1,s2] = c
			CONVERGENT[s2,s1] = c
			ANCESTRAL[s1,s2] = a
			ANCESTRAL[s2,s1] = a
			RECOMBINANT[s1,s2] = r
			RECOMBINANT[s2,s1] = r
			RATES[s1,s2] = r/float(pair_distances['distance_1'])
			RATES[s2,s1] = r/float(pair_distances['distance_2'])
			total += r/float(pair_distances['distance_1'])
			total += r/float(pair_distances['distance_2'])


			print('\n\nCompleted strain pairing ' + str(count) + ' out of ' + str(total_pairs) + '\n\n')
			count += 1
			
	average_rate = total / total_pairs

	return {'strain_names': strain_names, 'Shared': SHARED, 'Convergent': CONVERGENT, 'Ancestral': ANCESTRAL, 'Recombinant': RECOMBINANT, 'Rates': RATES, 'average': average_rate}


def rename_ancestors(internal_nodes, strain_names, tree_with_ancestors):
	# print(tree_with_ancestors)
	old_internal_nodes = dict(internal_nodes)
	for node in old_internal_nodes:
		if ('Q' + node) not in strain_names:
			internal_nodes['Q' + node] = old_internal_nodes[node]
			node_location = tree_with_ancestors.find(node)
			# print(node)
			# print(tree_with_ancestors[node_location:node_location+len(node)])
			left_char = tree_with_ancestors[node_location-1]
			right_char = tree_with_ancestors[node_location+len(node)]
			# print('left_char:' + left_char)
			# print('right_char:' + right_char)
			while left_char not in ['(',')',',',':'] or right_char not in ['(',')',',',':',';']:  # != '(' and right_char != ')' and right_char != ':' and right_char != ',' and right_char != ';': # makes sure this is the not a strain with the same name just longer (i.e. B5 and B57)
				node_location = tree_with_ancestors.find(node,node_location + len(node) + 1)
				left_char = tree_with_ancestors[node_location-1]
				right_char = tree_with_ancestors[node_location+len(node)]
				# print('left_char:' + left_char)
				# print('right_char:' + right_char)
			left = tree_with_ancestors[:node_location]
			right = tree_with_ancestors[node_location + len(node):]
			middle = 'Q' + node
			# print('left:' + left )
			# print('right:' + right)
			# print('middle: ' + middle)
			tree_with_ancestors = left + middle + right
			# print(tree_with_ancestors)

		# elif ('Z' + node) not in strain_names:
		# 	internal_nodes['Z' + node] = internal_nodes[node]
		# elif('X' + node) not in strain_names:
		# 	internal_nodes['X' + node] = internal_nodes[node]

		# del internal_nodes[node]

	return internal_nodes,tree_with_ancestors


	

# STEPS
# 1. take in the extant strain alignments
# 2. get the tree string from raxml
def get_tree_string(species_alignment, raxml_path):
	# subprocess.Popen(args, bufsize=0, executable=None, stdin=None, stdout=None, stderr=None, preexec_fn=None, close_fds=False, shell=False, cwd=None, env=None, universal_newlines=False, startupinfo=None, creationflags=0)
	raxml_command = 'raxmlHPC -s ' + species_alignment + ' -m GTRGAMMA -n tree -p 12345'
	subprocess.Popen(raxml_command.split(), cwd = raxml_path)
	print('Obtained the tree.')
	# subprocess.check_call(command.split())

# get_tree_string(species_alignment, raxml_path)

# 3. root it 
def get_tree_root(tree_file, raxml_path):
	print('here')
	raxml_command = 'raxmlHPC ­-f I -­t ' + tree_file + ' -­m GTRGAMMA -­n root'
	subprocess.Popen(raxml_command.split(), cwd = raxml_path)
	print('Rooted the tree.')

# get_tree_root(tree_file, raxml_path)


# 4. get the ancestral alignments
	# must use the file without the '-'s
def get_ancestors(rooted_tree_file, species_alignment, raxml_path, reduced):
	if not reduced:
		raxml_command = 'raxmlHPC ­-f A -­t ' + rooted_tree_file + ' -­s ' + species_alignment + ' -­m GTRGAMMA ­-n anc'
	else:
		raxml_command = 'raxmlHPC ­-f A -­t ' + rooted_tree_file + ' -­s ' + species_alignment + '.reduced -­m GTRGAMMA ­-n anc'
	subprocess.Popen(raxml_command.split(), cwd = raxml_path)
	print('Obtained the ancestral alignments.')

# get_ancestors(rooted_tree_file, species_alignment, raxml_path)

# 5. get the S & A matrices
def get_s_a(genome1, genome2, MRCA, L):
	# print(len(genome1))
	# print(len(genome2))
	# print(len(MRCA))
	l = min(len(MRCA), L)
	s,a = 0,0 # initializes the shared and ancestral values for the pair of strains to 0
	for site in range(l): # goes through every site along the genome
		if genome1[site] == genome2[site]: # counts up the number of shared sites
			s += 1
			if genome1[site] == MRCA[site]: # counts up the number of shared sites that were inherited from the ancestor
				a += 1

	return s,a


# # 6. scale the branch lengths
# def scale_branch_lengths(tree_string, min_m, max_m, real_pi, real_theta, kappa):
# 	scaled_trees = []
# 	average_pi = []
# 	average_theta = []
# 	ms = list(range(min_m, max_m+1))
# 	for m in ms:
# 		scaled_tree_string = scale_newick_format_tree(strains, L, m, max_m, tree_string)
# 		scaled_trees.append(scaled_tree_string)
# 		pis = []
# 		thetas = []

# 		phylogeny = pyvolve.read_tree(tree = scaled_tree_string)
# 		# pyvolve.print_tree(phylogeny)

# 		freqs = [0.25,0.25,0.25,0.25]
# 		nuc_model = pyvolve.Model('nucleotide', {'kappa':kappa, 'state_freqs':freqs})

# 		ancestor = generate_ancestor(L)
# 		for i in ranage(5):
# 			my_partition = pyvolve.Partition(models = nuc_model, root_sequence = ancestor)
# 			my_evolver = pyvolve.Evolver(partitions = my_partition, tree = phylogeny)
# 			my_evolver(ratefile = None, infofile = None, seqfile = "simulated_alignment_" + m + "_universal_" + str(iteration + 1) + ".fasta" )
# 			simulated_strains = my_evolver.get_sequences()
# 			# strains = my_evolver.get_sequences(anc = True)
# 			pi = pi_value(simulated_strains)
# 			theta = theta_value(simulated_strains)
# 			pis.append(pi)
# 			thetas.append(theta)

# 			average_pi.append(np.mean(pis))
# 			average_theta.append(np.mean(thetas))

# 	pi_margins = []
# 	theta_margins = []
# 	for item in range(len(ms)):
# 		pi_margins[item] = average_pi[item] - real_pi
# 		theta_margins[item] = average_theta[item] - real_theta

# 	min_margin, index = min((val, idx) for (idx, val) in enumerate(theta_margins))

# 	accurate_tree = scaled_trees[index]

# 	return accurate_tree

def scale_newick_format_tree(strains, L, desired_m, tree_string):
	l = len(tree_string)
	branch_lengths = []
	current = 0
	while current < l:
		start = tree_string.find(':', current) + 1
		if start == 0:
			current = l
		else:
			x = tree_string.find(',', start)
			y = tree_string.find(')', start)
			if x == -1:
				end = y
			elif y == -1:
				end = x
			else:
				end = min(x,y)
			# if end == -1:
			# 	end = -2
			branch_lengths.append([start, end])
			current = end
	scaled_tree_string = tree_string
	for branch in branch_lengths:
		branch_length = tree_string[branch[0]:branch[1]]
		scaled_tree_string = scaled_tree_string.replace(branch_length, scale(branch_length, desired_m, max_m))

	return scaled_tree_string

def scale(branch_length, desired_m, max_m):
	scaled_branch_length = str(float(branch_length) * (int(min_m)) / max_m) # (float(total_branch_length) * int(L)))
	return scaled_branch_length


def merge_trees(tree_with_distances, tree_with_nodes):

	left = tree_with_distances[:-1]
	tree_with_distances = left + 'ROOT:0.0;'

	current_1 = 0
	current_2 = 0
	l_1 = len(tree_with_distances)
	l_2 = len(tree_with_nodes)

	while current_1 != -1:
		# while current_2 < l_2:

		current_1 = tree_with_distances.find(':', current_1+1)
		if(tree_with_distances[current_1-1] == ',' or tree_with_distances[current_1-1] == '(' or tree_with_distances[current_1-1] == ')'):
			left = tree_with_distances[:current_1]
			right = tree_with_distances[current_1:]
			left_count = left.count('(')
			right_count = left.count(')')
			comma_count = left.count(',')
			while (left_count != 0 or right_count != 0 or comma_count != 0):
				if tree_with_nodes[current_2] == '(':
					left_count -= 1
				elif tree_with_nodes[current_2] == ')':
					right_count -= 1
				elif tree_with_nodes[current_2] == ',':
					comma_count -= 1
				current_2 += 1
			name_start = current_2
			while not (tree_with_nodes[current_2] == '(' or tree_with_nodes[current_2] == ')' or tree_with_nodes[current_2] == ','):
				current_2 += 1
			name_end = current_2

			middle = tree_with_nodes[name_start:name_end]
			tree_with_distances = left + middle + right
			current_2 = 0
			# print(tree_with_distances)
	return tree_with_distances


	# current = 0 # the current position along the tree_string we are accessing
	# count = 1 # counter for the number of new nodes that have been added, for naming purposes
	# new_nodes = [] # will be populated as a list of the node names that get added 
	# while current != -1: # continues while the end of the string has not been reached
	# 	# print(tree_string)
	# 	current = tree_string.find(':', current+1) # finds the next ':' in the tree_string, the characters proceeding the ':' should always be the strain name
	# 	# print(tree_string[current-1])
	# 	if(tree_string[current-1] == ',' or tree_string[current-1] == '(' or tree_string[current-1] == ')'): # will return True if the node is not yet named and False if it is already named 
	# 		# print(tree_string)
	# 		left = tree_string[:current] # the portion of the string to the left of the current ':'
	# 		right = tree_string[current:] # the portion of the string to the right of the current ':'
	# 		while 'Q' + str(count) in strain_names: # will continue as long as the name we plan to use is already used for another node, will terminate as soon as it generates a unique name
	# 			count += 1
	# 		middle = 'Q' + str(count) # the new node name
	# 		tree_string = left + middle + right # the updated tree_string
	# 		new_nodes.append(middle) # keeps track of the node that was added
	# 		count += 1
	# return {'tree_string': tree_string, 'new_nodes': new_nodes}

# # 7. get the sequences of parents 
# def get_parents():

# # 8. get the MRCA
# def get_MRCA():

# 9. get the C matrix
def get_c(strain1, strain2, MRCA, parents, tree_string, distances, L, mu, kappa):
	pair_distances = get_distances_to_MRCA(strain1, strain2, MRCA, tree_string, parents, distances) # gets the total lengths of the branches back to the MRCA of strain 1 and strain 2
	distance_1 = pair_distances['distance_1'] # the distance from strain 1 to the MRCA
	distance_2 = pair_distances['distance_2'] # the distance from strain 2 to the MRCA
	m_1 = int(distance_1 * L) # the number of mutations that occurred on strain 1
	m_2 = int(distance_1 * L) # the number of mutations that occurred on strain 2
	###################### CHANGE #####################
	generations_1 = int(m_1/mu) # the number of generations over which these mutations occurred on strain 1
	generations_2 = int(m_2/mu) # the number of generations over which these mutations occurred on strain 2
	# print('m_1 = ' + str(m_1))
	# print('m_2 = ' + str(m_2))
	# print('generations_1 = ' + str(generations_1))
	# print('generations_2 = ' + str(generations_2))

	c = expected_c_given_ms(L, m_1, m_2, mu, generations_1, generations_2, kappa, 0.5) # the expected number of convergent mutations between strain 1 and strain 2

	return c, pair_distances


# # 10. get the R matrix
# def get_R():

# 11. print the matrices to a csv
def print_matrices(output_file, S,C,A,R,RATES, strain_names, average_rate):
	shape = S.shape

	with open((output_file + '.csv'), 'w', newline = '') as f: 
		writer = csv.writer(f)
		# writer.writerow([species[x][:-1]])
		writer.writerow(['SHARED'])
		header = ['']
		header.extend(strain_names)
		writer.writerow(header)
		for row in range(shape[0]):
			write_row = [strain_names[row]]
			write_row.extend(S[row])
			writer.writerow(write_row)

		writer.writerow(['CONVERGENT'])
		writer.writerow(header)
		for row in range(shape[0]):
			write_row = [strain_names[row]]
			write_row.extend(C[row])
			writer.writerow(write_row)

		writer.writerow(['ANCESTRAL'])
		writer.writerow(header)
		for row in range(shape[0]):
			write_row = [strain_names[row]]
			write_row.extend(A[row])
			writer.writerow(write_row)

		writer.writerow(['RECOMBINANT'])
		writer.writerow(header)
		for row in range(shape[0]):
			write_row = [strain_names[row]]
			write_row.extend(R[row])
			writer.writerow(write_row)

		writer.writerow(['RATES'])
		writer.writerow(header)
		for row in range(shape[0]):
			write_row = [strain_names[row]]
			write_row.extend(RATES[row])
			writer.writerow(write_row)
		
		writer.writerow([average_rate])


# # 12. print average recombs? with stddev?