
# filename = 

import operator
import pyvolve
from simulations import generate_ancestor
from process_genomes import pi_value
from process_genomes import theta_value
import numpy as np

def scale_branch_lengths(L, tree_string, min_m, max_m, real_pi, real_theta, kappa):
	scaled_trees = []
	average_pi = []
	average_theta = []
	ms = list(range(min_m, max_m+1))
	for m in ms:
		print('trying all m values')
		scaled_tree_string = scale_newick_format_tree(m, max_m, tree_string)
		scaled_trees.append(scaled_tree_string)
		pis = []
		thetas = []

		phylogeny = pyvolve.read_tree(tree = scaled_tree_string)
		pyvolve.print_tree(phylogeny)

		freqs = [0.25,0.25,0.25,0.25]
		nuc_model = pyvolve.Model('nucleotide', {'kappa':kappa, 'state_freqs':freqs})

		ancestor = generate_ancestor(L)
		###############CHANGE THIS##############
		for i in range(1):
			my_partition = pyvolve.Partition(models = nuc_model, root_sequence = ancestor)
			my_evolver = pyvolve.Evolver(partitions = my_partition, tree = phylogeny)
			my_evolver(ratefile = None, infofile = None, seqfile = "simulated_alignment_" + str(m) + "_universal_" + str(i + 1) + ".fasta" )
			print('simulated')
			simulated_strains = my_evolver.get_sequences()
			# strains = my_evolver.get_sequences(anc = True)
			pi = pi_value(simulated_strains)
			theta = theta_value(simulated_strains)
			pis.append(pi)
			thetas.append(theta)

		average_pi.append(np.mean(pis))
		average_theta.append(np.mean(thetas))
	print('got all m values')

	pi_margins = []
	theta_margins = []
	for item in range(len(ms)):
		pi_margins[item] = average_pi[item] - real_pi
		theta_margins[item] = average_theta[item] - real_theta

	print(theta_margins)

	min_theta_margin_index, min_theta_margin = min(enumerate(theta_margins), key=operator.itemgetter(1))
	min_pi_margin_index, min_pi_margin = min(enumerate(pi_margins), key=operator.itemgetter(1))

	# min_margin, index = min((val, idx) for (idx, val) in enumerate(theta_margins))
	index_difference = min_theta_margin_index - min_pi_margin_index
	if index_difference == 0:
		print('theta and pi have both been optimized to within ' + min_theta_margin + ' and ' + min_pi_margin + ', respectively.')
	elif index_difference < 0:
		print('Theta has been optimized to within ' + min_theta_margin + ', but to optimize pi would require ' + index_difference + ' more mutations.')
	elif index_difference > 0:
		print('Theta has been optimized to within ' + min_theta_margin + ', but to optimize pi would require ' + index_difference + ' fewer mutations.')

	# print('found most accurate m')

	accurate_tree = scaled_trees[min_theta_margin_index]

	return accurate_tree

def scale_newick_format_tree(desired_m, max_m, tree_string):
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
			if end == -1:
				end = -2
			branch_lengths.append([start, end])
			current = end
	scaled_tree_string = tree_string
	for branch in branch_lengths:
		branch_length = tree_string[branch[0]:branch[1]]
		scaled_tree_string = scaled_tree_string.replace(branch_length, scale(branch_length, desired_m, max_m))

	return scaled_tree_string

def scale(branch_length, desired_m, max_m):
	scaled_branch_length = str(float(branch_length) * (int(desired_m)) / max_m) # (float(total_branch_length) * int(L)))
	return scaled_branch_length

# def read_in_strains(filename):
# 	f = open(filename, 'r')
# 	f = list(f)
# 	strains = {} # dictionary to hold all the strains (key: strain name, value: strain genome)

# 	for line in f: 
# 		if(line[0] == '>'): # separates out the strain names
# 			key = line.strip('\n')
# 			strains[key] = ''
# 		else:
# 			strains[key] += line.strip('\n') # concatenates all the lines of genome
# 	values = list(strains.values())
# 	length = len(values[0])
# 	for strain in strains.values(): # makes sure that the genome length of each strain is uniform
# 		# print(len(strain))
# 		if len(strain) != length:
# 			print('ABORT: THE LENGTHS ARE NOT THE SAME')
# 	return strains

# def genome_length(species):
# 	strains = list(species.values())
# 	length = len(strains[0])
# 	for i in range(1, len(strains)):
# 		if len(strains[i]) != length:
# 			raise ValueError('The strains do not all have the same genome length')
# 	return length

def get_min_m(strains, L):
	print('Determining the minimum number of mutations.\n')
	sequences = list(strains.values())

	alleles = L*[None]
	for site in range(L):
		alleles[site] = []

	for site in range(L):
		for sequence in sequences:
			if sequence[site] in ['A', 'T', 'G', 'C'] and sequence[site] not in alleles[site]:
				alleles[site].append(sequence[site])
		# print(site)

	min_m = 0
	for site in range(L):
		min_m += (len(alleles[site]) - 1)

	print('complete')

	return min_m

def get_max_m(strains, L, tree_string):
	l = len(tree_string)
	total_branch_length = 0
	# print(l)
	# branch_lengths = []
	current = 0
	while current < l:
		start = tree_string.find(':', current) + 1
		if start == 0:
			current = l
		else:
			# print('start = ' + str(start))
			x = tree_string.find(',', start)
			y = tree_string.find(')', start)
			if x == -1:
				end = y
			elif y == -1:
				end = x
			else:
				end = min(x,y)
			if end == -1:
				end = -2
			# print(end)
			# print('end = ' + str(end))
			# branch_lengths.append([start, end])
			# print(tree_string[start:end])
			total_branch_length += float(tree_string[start:end])
			# print(tree_string[start:end])
			current = end
			# print('current = ' + str(current))
	print('found max_m')
	return int((float(total_branch_length) * int(L)))

# increment = True/False
# def scale(branch_length, total_branch_length, L, desired_m, max_m):
# 	# max_m = (float(total_branch_length) * int(L))
# 	scaled_branch_length = str(float(branch_length) * (int(min_m)) / max_m) # (float(total_branch_length) * int(L)))
# 	# print(scaled_branch_length)
# 	return scaled_branch_length
# 	# return min_m / L

# def scale_newick_format_tree(strains, L, desired_m, tree_string):
# 	print('Scaling the branch lengths.\n')
# 	l = len(tree_string)
# 	# total_branch_length = 0
# 	# print(l)
# 	branch_lengths = []
# 	max_m = get_max_m(strains, L, tree_string)
# 	current = 0
# 	while current < l:
# 		start = tree_string.find(':', current) + 1
# 		if start == 0:
# 			current = l
# 		else:
# 			# print('start = ' + str(start))
# 			x = tree_string.find(',', start)
# 			y = tree_string.find(')', start)
# 			if x == -1:
# 				end = y
# 			elif y == -1:
# 				end = x
# 			else:
# 				end = min(x,y)
# 			if end == -1:
# 				end = -2
# 			# print(end)
# 			# print('end = ' + str(end))
# 			branch_lengths.append([start, end])
# 			# print(tree_string[start:end])
# 			# total_branch_length += float(tree_string[start:end])
# 			# print(tree_string[start:end])
# 			current = end
# 			# print('current = ' + str(current))
# 	scaled_tree_string = tree_string
# 	for branch in branch_lengths:
# 		branch_length = tree_string[branch[0]:branch[1]]
# 		# print(length)
# 		scaled_tree_string = scaled_tree_string.replace(branch_length, scale(branch_length, total_branch_length, L, desired_m, max_m))
# 	# for x in range(10000):
# 	# 	print(scaled_tree_string)
# 	return scaled_tree_string


	# str.index(str, beg = 0 end = len(string))

	# str.replace(old, new[, count])
	# Return a copy of the string with all occurrences of substring old replaced by new. If the optional argument count is given, only the first count occurrences are replaced.

# strains = {1:'AAAA', 2:'TTTT', 3:'GGGG', 4:'CCCC'}
# L = 4
# min_m = get_min_m(strains, L)
# # print(min_m)
# tree_string = '((t1:0.5,t2:0.5)i1:0.5,(t3:0.5,t4:0.5)i2:0.5)root;'

# print(scale_newick_format_tree(strains, L, min_m, tree_string))