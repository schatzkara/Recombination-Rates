
# script to apply the model given m1 and m2 along the phylogeny and generate a matrix of c values

import numpy as np
from process_genomes import read_in_strains
from get_internal_nodes import get_internal_nodes
from process_genomes import species_size
from process_genomes import genome_length
from process_genomes import theta_value
from completely_new_thing import get_min_m
from completely_new_thing import scale_newick_format_tree
from model import expected_c_given_ms


# function to apply our probabilistic model along all the branches of the phylogeny and get S = C + A + R 
# params:
# 	species_path (string) = complete path where the .fa file containing the species' genomes can be found
# 	kappa (float) = the ratio of transitions to transversions specifically for this species
# 	tree_string (string) = a Newick formatted phylogenetic tree of the species
# return: a dictionary with the strain names and each of the S,C,A,R matrices; key = 'strain_names', the matrix names, value = the list of names, the matrices
def apply_model_along_phylogeny(species_path, kappa, tree_string):
	# ancestor = ''
	# print('Reading in the strains.\n\n')
	strains = read_in_strains(species_path) # dictionary with the genomes of all the strains; key = strain name, value = genome
	internal_nodes = get_internal_nodes(anc_path)
	all_nodes = strains + internal_nodes
	strain_names = list(strains.keys()) # list of all the extant strain names
	all_node_names = list(all_nodes.keys())
	print(strain_names)
	n = species_size(strains) # number of extant strains
	total_pairs = (n*(n-1))/2 # the total number of strain pairs that will be compared
	L = genome_length(strains) # number of base pairs in the genomd
	theta = theta_value(strains) # proportion of the genome that is polymorphic
	mu = theta/(2*n) # mutation rate in mutations per base pair per generation
	min_m = get_min_m(strains, L) # minimum number of mutations that could account for all the polymorphisms in the species
	# print('Scaling the branch lengths of the tree.\n\n')
	# scaled_tree_string = tree_string
	scaled_tree_string = scale_newick_format_tree(strains, L, min_m, tree_string, 0) # the tree_string scaled by min_m

	SHARED = np.empty([n,n], dtype = np.float, order='C') # a matrix of the number of nucleotides shared between two strains; the (i,j) entry is the number of nucleotides that are the same between strain i and strain j
	CONVERGENT = np.empty([n,n], dtype = np.float, order='C') # a matrix of the number of nucleotides that match due to convergent mutation between two strains; the (i,j) entry is the number of convergent mutations between strain i and strain j
	ANCESTRAL = np.empty([n,n], dtype = np.float, order='C') # a matrix of the number of nucleotides that match due to direct inheritence from the ancestor; the (i,j) entry is the number of nucleotides that were inherited by both strain i and strain j
	RECOMBINANT = np.empty([n,n], dtype = np.float, order='C') # a matrix of the number of nucleotides that match due to a recombination event; the (i,j) entry is the number of nucleotides that were recombined between strain i and strain j

	updated_tree_info = name_nodes(tree_string, strain_names) 
	tree_string = updated_tree_info['tree_string'] # version of the tree_string where every node is labeled
	new_nodes = updated_tree_info['new_nodes'] # the new node names that were added
	all_nodes = all_node_names + new_nodes # list of all the node names in the pyhlogenetic tree
	# print(all_nodes) 

	# parents = find_parents(strain_names, tree_string) # a dictionary of the sequence of parents of each strain; key = strain name, value = list of the parents in order of increasing distance from the strain
	parents = find_parents(all_node_names, tree_string)
	distances = get_branch_lengths(all_nodes, tree_string) # a dictionary of the distances of each strain to its closest ancestor; key = strain name, value = distance to its closest ancestor
	# print(parents)
	print(distances)

	# parents = parents_and_distances['parents']
	# distances = parents_and_distances['distances']

	count = 1 # a counter for the current strain pair number that is being processed
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

			s,a = 0,0 # initializes the shared and ancestral values for the pair of strains to 0
			for site in range(L): # goes through every site along the genome
				if genome1[site] == genome2[site]: # counts up the number of shared sites
					s += 1
					if genome1[site] == ancestor[site]: # counts up the number of shared sites that were inherited from the ancestor
						a += 1

			# s1_tree_location = scaled_tree_string.find(strain_names[s1])
			
			# s2_tree_location = scaled_tree_string.find(strain_names[s2])

			# start_length_1 = scaled_tree_string.find(':', s1_tree_location) + 1
			# x1 = scaled_tree_string.find(',', start_length_1)
			# y1 = scaled_tree_string.find(')', start_length_1)
			# if x1 == -1:
			# 	end_length_1 = y1
			# elif y1 == -1:
			# 	end_length_1 = x1
			# else:
			# 	end_length_1 = min(x1,y1)

			# start_length_2 = scaled_tree_string.find(':', s2_tree_location) + 1
			# x2 = scaled_tree_string.find(',', start_length_2)
			# y2 = scaled_tree_string.find(')', start_length_2)
			# if x2 == -1:
			# 	end_length_2 = y2
			# elif y2 == -1:
			# 	end_length_2 = x2
			# else:
			# 	end_length_2 = min(x2,y2)

			# length_1 = float(scaled_tree_string[start_length_1:end_length_1])
			# length_2 = float(scaled_tree_string[start_length_2:end_length_2])
			# MRCA = find_MRCA(strain1, strain2, parents) # the Most Recent Common Ancestor between the two strains
			pair_distances = get_distances_to_MRCA(strain1, strain2, MRCA, tree_string, strain_names, parents, distances) # gets the total lengths of the branches back to the MRCA of strain 1 and strain 2
			distance_1 = pair_distances['distance_1'] # the distance from strain 1 to the MRCA
			distance_2 = pair_distances['distance_2'] # the distance from strain 2 to the MRCA
			m_1 = int(distance_1 * L + 1) # the number of mutations that occurred on strain 1
			m_2 = int(distance_1 * L + 1) # the number of mutations that occurred on strain 2
			generations_1 = int(m_1/mu) # the number of generations over which these mutations occurred on strain 1
			generations_2 = int(m_2/mu) # the number of generations over which these mutations occurred on strain 2

			c = expected_c_given_ms(L, m_1, m_2, mu, generations_1, generations_2, kappa, 0.5) # the expected number of convergent mutations between strain 1 and strain 2

			# fills in the appropriate values to the S,C,A,R matrices for the current strain pair
			SHARED[s1,s2] = s
			SHARED[s2,s1] = s
			CONVERGENT[s1,s2] = c
			CONVERGENT[s2,s1] = c
			ANCESTRAL[s1,s2] = a
			ANCESTRAL[s2,s1] = a
			RECOMBINANT[s1,s2] = s - a - c
			RECOMBINANT[s2,s1] = s - a - c

			count += 1
			print('\n\nCompleted strain pairing ' + str(count) + ' out of ' + str(n**2) + '\n\n')

	# return {'strain_names': strain_names, 'Convergent': CONVERGENT}
	return {'strain_names': strain_names, 'Shared': SHARED, 'Convergent': CONVERGENT, 'Ancestral': ANCESTRAL, 'Recombinant': RECOMBINANT}

# a function to provide a name to each node along the Newick formatted tree_string
# params:
# 	tree_string (string) = a Newick formatted phylogenetic tree of the species
# 	strain_names (list) = a list of the names of all the strains in tree
# return: a dictionary as follows: {'tree_string': tree_string, 'new_nodes': new_nodes}
# 	tree_string (string) = the new tree string with fully labeled nodes
# 	new_nodes (list) = a list of the node names that were added
def name_nodes(tree_string, strain_names):
	# print('Making sure all the nodes are named.\n')
	current = 0 # the current position along the tree_string we are accessing
	count = 1 # counter for the number of new nodes that have been added, for naming purposes
	new_nodes = [] # will be populated as a list of the node names that get added 
	while current != -1: # continues while the end of the string has not been reached
		# print(tree_string)
		current = tree_string.find(':', current+1) # finds the next ':' in the tree_string, the characters proceeding the ':' should always be the strain name
		# print(tree_string[current-1])
		if(tree_string[current-1] == ',' or tree_string[current-1] == '(' or tree_string[current-1] == ')'): # will return True if the node is not yet named and False if it is already named 
			# print(tree_string)
			left = tree_string[:current] # the portion of the string to the left of the current ':'
			right = tree_string[current:] # the portion of the string to the right of the current ':'
			while 'Q' + str(count) in strain_names: # will continue as long as the name we plan to use is already used for another node, will terminate as soon as it generates a unique name
				count += 1
			middle = 'Q' + str(count) # the new node name
			tree_string = left + middle + right # the updated tree_string
			new_nodes.append(middle) # keeps track of the node that was added
			count += 1
	return {'tree_string': tree_string, 'new_nodes': new_nodes}

# a function to get the distance of each strain to its closest ancestor
# params:
# 	strain_names (list) = a list of the names of all the strains in tree
# 	tree_string (string) = a Newick formatted phylogenetic tree of the species
# return: distances (dict) = a dictionary of the distances of each strain to its closest ancestor; key = strain name, value = distance to its closest ancestor
def get_branch_lengths(strain_names, tree_string):
	# parents = {} # dictionary of all the parents for each strain; key = strain name, value = list of the parents in order of increasing distance from the strain
	distances = {} # dictonary of all the distances to the closest ancestors; key = strain name, value = distance to its closest ancestor
	for child in strain_names: # iterates over all strains given
		# print(child)
		child_location = tree_string.find(child) # the index of the starting position of the strain name in the tree_string
		while tree_string[child_location+len(child)] != ':': # makes sure this is the not a strain with the same name just longer (i.e. B5 and B57)
			child_location = tree_string.find(child,child_location + len(child) + 1)
		child_location_end = child_location + len(child) # the index of the ending position of the strain name in the tree_string
		length_start = child_location_end + 1 # the index of the starting position of the branch length
		x = tree_string.find(',', length_start) # the index of the next ','
		y = tree_string.find(')', length_start) # the index of the next ')'
		if x == -1: # if there are no more ','s, then the length ends at the ')'
			length_end = y
		elif y == -1: # if there are no more ')'s, then the length ends at the ')'
			length_end = x
		else: # otherwise, the length ends at whichever of ',' and ')' is reached first (or it ends at the end of the string (-1); both x and y are -1 in this case)
			length_end = min(x,y)
		length = tree_string[length_start:length_end] # the branch length associated with the current node
		distances[child] = length
		# print(child)
	return distances

# a function to determine the sequence of parents for each strain in the tree and the distance of each strain to its ancestor
# params:
# 	strain_names (list) = a list of the names of all the strains in tree
# 	tree_string (string) = a Newick formatted phylogenetic tree of the species
# return: a dictionary of the sequence of parents of each strain; key = strain name, value = list of the parents in order of increasing distance from the strain
def find_parents(strain_names, tree_string):
	parents = {} # dictionary of all the parents for each strain; key = strain name, value = list of the parents in order of increasing distance from the strain
	# distances = {} # dictonary of all the distances to the closest ancestors; key = strain name, value = distance to its closest ancestor
	for child in strain_names: # iterates over all strains given
		child_location = tree_string.find(child) # the index of the starting position of the strain name in the tree_string
		if tree_string[child_location+len(child)] != ':': # makes sure this is the not a strain with the same name just longer (i.e. B5 and B57)
			child_location = tree_string.find(child,child_location + len(child) + 1)
		child_location_end = child_location + len(child) # the index of the ending position of the strain name in the tree_string
	# 	length_start = child_location_end + 1 # the index of the starting position of the branch length
	# 	x = tree_string.find(',', length_start) # the index of the next ','
	# 	y = tree_string.find(')', length_start) # the index of the next ')'
	# 	if x == -1: # if there are no more ','s, then the length ends at the ')'
	# 		length_end = y
	# 	elif y == -1: # if there are no more ')'s, then the length ends at the ')'
	# 		length_end = x
	# 	else: # otherwise, the length ends at whichever of ',' and ')' is reached first (or it ends at the end of the string (-1); both x and y are -1 in this case)
	# 		length_end = min(x,y)
	# 	length = tree_string[length_start:length_end] # the branch length associated with the current node
	# 	distances[child] = length

		parent_tally = 0 # tally for the number of parents that the strain has
		for character in range(0,child_location): # counts up the number of '(' present before the strain that do not have a closing ')'; this is the number of parents that the strain has
			parent_tally += (tree_string[character] == '(')
			parent_tally -= (tree_string[character] == ')')
		parents[child] = [child] # puts the child as the first 'parent' so that every possible pairing of strains can occur later

		found = 0 # counter for the number of parents that have been found so far
		left_tally = 0 # tally for the number of extra '(' that are seen when searching for the parents; these are other internal subtrees that must be ignored
		current = child_location + 1 # the current position along the tree we are accessing
		while found != parent_tally: # continues while all the parents have not been found; terminates once they have all been found
			while current < len(tree_string): # makes sure that the end of the string is not reached
				# print(current)
				if tree_string[current] == '(': # marks that an extra '(' has been seen
					left_tally += 1
					current += 1
				elif tree_string[current] == ')': # finds where ')'s occur
					if left_tally != 0: # if there are extra '(' that have not been closed, this closes it
						left_tally -= 1
						current += 1
					else: # otherwise, the parent has been found
						parent_start = current + 1 # the index of the starting position of the parent node
						parent_end = tree_string.find(':', parent_start) # the index of the ending position of the parent node
						parents[child].append(tree_string[parent_start:parent_end]) # adds the parent to the list
						found += 1
						current += 1
				else: # passes over any other characters in the string
					current += 1
			# print(children)
	return parents
	# return {'parents': parents, 'distances': distances}
			
# tree_string = '' # '(A1:0.1,A11:0.2,(A111:0.3,D:0.4)E:0.5)F;'
# # '((((D:0.5,(E:0.5,F:0.5)T:0.5)U:0.5,(G:0.5,H:0.5)R:0.5)S:0.5,(A:0.5,(B:0.5,C:0.5)V:0.5)W:0.5)Y:0.5,((I:0.5(J:0.5,K:0.5)PP:0.5)QQ:0.5,(L:0.5,(M:0.5,(N:0.5,(O:0.5,(P:0.5,Q:0.5)KK:0.5)LL:0.5)MM:0.5)NN:0.5)OO:0.5)X:0.5)Z:0.5;' 
# # '(F:0.5,((C:0.5,D:0.5):0.5,E:0.5)H:0.5,(A:0.5,B:0.5)G:0.5)I:0.5;' 
# # '(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;'
# # '(A1:0.1,A11:0.2,(A111:0.3,D:0.4)E:0.5)F;'
# tree_string = name_nodes(tree_string)
# print(find_parents(['A1','A11','A111','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z','QQ','PP','OO','NN','MM','LL','KK'], tree_string))


			# if tree_string[current] == '(':
			# 	current = tree_string.find(')', current)
			# elif tree_string[current] == ')':
			# 	parent_start = current + 1
			# 	parent_end = tree_string.find(':', parent_start)
			# 	children[child][found] = tree_string[parent_start:parent_end]
			# 	found += 1
			# else:
			# 	current += 1

# a function to find the Most Recent Common Ancestor between any two strains
# params:
# 	strain1, strain2 (strings) = the two strains between which to find the MRCA
# 	parents (dict) = a dictionary of the sequence of parents of each strain; key = strain name, value = list of the parents in order of increasing distance from the strain
# return: a string corresponding to the MRCA between the two strains 
def find_MRCA(strain1, strain2, parents):
	parents_1 = parents[strain1] # the parents of strain 1
	parents_2 = parents[strain2] # the parents of strain 2
	found = False # False if the MRCA has not yet been found, True if it has
	# i = 0
	
	if len(parents_1) < len(parents_2): # goes through the shorter of the two parents lists
		i = 0 # the index of the current parent that is being checked
		while not found: # continues while the MRCA has not been found, terminates once it has
			if parents_1[i] in parents_2: # if the parent is also a parent of the other strain, then the MRCA has been found
				MRCA = parents_1[i]
				found = True
			i += 1 # otherwise, the next parent needs to be checked
	else:
		i = 0
		while not found:
			if parents_2[i] in parents_1:
				MRCA = parents_2[i]
				found = True
			i += 1

	return MRCA

# a function to get the distances of two strains to their MRCA
# params: 
# 	strain1, strain2 (strings) = the two strains between which to find distances to the MRCA
# 	tree_string (string) = a Newick formatted phylogenetic tree of the species
# 	strain_names (list) = a list of the names of all the strains in tree
# return: a dictionary with the distances from strain 1 and strain 2 to their MRCA; keys: 'distance_1' and 'distance_2'
def get_distances_to_MRCA(strain1, strain2, MRCA, tree_string, parents, distances):
	print('Finding the distances to the MRCA.\n')

	distance_1 = 0 # float(distances[strain1]) # the sum of the distances from strain 1 to the MRCA
	distance_2 = 0 # float(distances[strain2]) # the sum of the distances from strain 2 to the MRCA

	parents_1 = parents[strain1] # a list of the parents of strain 1
	parents_2 = parents[strain2] # a list of the parents of strain 2

	separation_1 = parents_1.index(MRCA) # the number of "generations" that separate strain 1 from the MRCA; e.g. the number of nodes between them, the number of branch lengths to sum up
	separation_2 = parents_2.index(MRCA) # the number of "generations" that separate strain 2 from the MRCA; e.g. the number of nodes between them, the number of brance lengths to sum up 

	for i in range(separation_1): # adds up all the branch lengths between strain 1 and the MRCA
		distance_1 += float(distances[parents_1[i]])
		# print(parents_1[i])
		# distance_1 += sum(parents_1[0:separation_1])
	for i in range(separation_2): # adds up all the branch lengths between strain 2 and the MRCA
		# print(parents_2[i])
		distance_2 += float(distances[parents_2[i]])

	return {'distance_1': distance_1, 'distance_2': distance_2}


# tree_string = '(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;'
# # '((((D:0.5,(E:0.5,F:0.5)T:0.5)U:0.5,(G:0.5,H:0.5)R:0.5)S:0.5,(A:0.5,(B:0.5,C:0.5)V:0.5)W:0.5)Y:0.5,((I:0.5(J:0.5,K:0.5)PP:0.5)QQ:0.5,(L:0.5,(M:0.5,(N:0.5,(O:0.5,(P:0.5,Q:0.5)KK:0.5)LL:0.5)MM:0.5)NN:0.5)OO:0.5)X:0.5)Z:0.5;' 
# # '(F:0.5,((C:0.5,D:0.5):0.5,E:0.5)H:0.5,(A:0.5,B:0.5)G:0.5)I:0.5;' 
# # '(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;'
# # '(A1:0.1,A11:0.2,(A111:0.3,D:0.4)E:0.5)F;'
# strain_names = ['A','B','C','D','E','F']
# for s1 in range(len(strain_names)):
# 	strain1 = strain_names[s1]
# 	for s2 in range(s1+1, len(strain_names)):
# 		strain2 = strain_names[s2]
# 		print(strain1)
# 		print(strain2)
# 		print(get_distances_to_MRCA(strain1,strain2,tree_string,strain_names))

	# s1_tree_location = tree_string.find(strain1)
	# s2_tree_location = tree_string.find(strain2)
	# distance_1 = 0
	# distance_2 = 0
	# start_MRCA_1 = 0
	# start_MRCA_2 = 1
	# i = 0
	# while start_MRCA_1 != start_MRCA_2:
	# 	print('iteration: ' + str(i))
	# 	start_length_1 = tree_string.find(':', s1_tree_location) + 1
	# 	end_length_1 = start_length_1
	# 	# print(end_length_1)
	# 	# print(tree_string[start_length_1:end_length_1])
	# 	while not ((tree_string[end_length_1] == ',') or (tree_string[end_length_1] == ')') or (tree_string[end_length_1] == ';')):
	# 	# while tree_string[end_length_1] != (',' or ')' or ';'):
	# 	# while tree_string[end_length_1] != ',' 
	# 		end_length_1 += 1
	# 		# print(str(end_length_1) + ',' + str(tree_string[end_length_1]))
	# 	# print(tree_string[start_length_1:end_length_1])
	# 	# x1 = tree_string.find(',', start_length_1)
	# 	# y1 = tree_string.find(')', start_length_1)
	# 	# if x1 == -1:
	# 	# 	end_length_1 = y1
	# 	# elif y1 == -1:
	# 	# 	end_length_1 = x1
	# 	# else:
	# 	# 	end_length_1 = min(x1,y1)
	# 	print('length_1 = ' + str(tree_string[start_length_1:end_length_1]))
	# 	distance_1 += float(tree_string[start_length_1:end_length_1])

	# 	current = end_length_1
	# 	left_count = 1
	# 	right_count = 0
	# 	temp1 = tree_string[current] == ')'
	# 	temp2 = tree_string[current+1] != ','
	# 	temp3 = left_count == right_count
	# 	print(str(temp1) + ',' + str(temp2) + ',' + str(temp3))
	# 	while not (temp1 and temp2 and temp3): # ((tree_string[current] == ')') and (tree_string[current+1] != ',') and (left_count == right_count)):
	# 		old_current = current 
	# 		current = tree_string.find(')', current)
	# 		for character in range(old_current,current+1):
	# 			left_count += (tree_string[character] == '(')
	# 			right_count += (tree_string[character] == ')')
	# 		# print(current)
	# 		temp1 = tree_string[current] == ')'
	# 		temp2 = tree_string[current+1] != ','
	# 		temp3 = left_count == right_count
	# 		print(str(temp1) + ',' + str(temp2) + ',' + str(temp3))
	# 	# check = 0
	# 	# while check: # not ( (tree_string[current] == ')') and (tree_string[current+1] != ',') and (left_count == right_count)):
	# 	# 	old_current = current
	# 	# 	current = tree_string.find(')', current+1)
	# 	# 	# print(current)
	# 	# 	for character in range(old_current,current):
	# 	# 		left_count += (tree_string[character] == '(')
	# 	# 		right_count += (tree_string[character] == ')')
	# 	# 	check = tree_string[current] == ')' and tree_string[current+1] != ',' and left_count == right_count
	# 	# 	print(current)
	# 	left_count = 1
	# 	right_count = 0
	# 	start_MRCA_1 = current + 1
	# 	end_MRCA_1 = tree_string.find(':', start_MRCA_1)
	# 	MRCA_1 = tree_string[start_MRCA_1:end_MRCA_1]


	# 	start_length_2 = tree_string.find(':', s2_tree_location) + 1
	# 	end_length_2 = start_length_2 
	# 	while not ((tree_string[end_length_2] == ',') or (tree_string[end_length_2] == ')') or (tree_string[end_length_2] == ';')):
	# 	# while tree_string[end_length_2] != (',' or ')' or ';'):
	# 		end_length_2 += 1
	# 	# x2 = tree_string.find(',', start_length_2)
	# 	# y2 = tree_string.find(')', start_length_2)
	# 	# if x2 == -1:
	# 	# 	end_length_2 = y2
	# 	# elif y2 == -1:
	# 	# 	end_length_2 = x2
	# 	# else:
	# 	# 	end_length_2 = min(x2,y2)

	# 	print('length_2 = ' + str(tree_string[start_length_2:end_length_2]))
	# 	distance_2 += float(tree_string[start_length_2:end_length_2])

	# 	current = end_length_2
	# 	print(current)
	# 	left_count = 1
	# 	right_count = 0
	# 	check = 0
	# 	temp1 = tree_string[current] == ')'
	# 	temp2 = tree_string[current+1] != ','
	# 	temp3 = left_count == right_count
	# 	print(str(temp1) + ',' + str(temp2) + ',' + str(temp3))
	# 	while not (temp1 and temp2 and temp3): # ((tree_string[current] == ')') and (tree_string[current+1] != ',') and (left_count == right_count)):
	# 		old_current = current 
	# 		current = tree_string.find(')', current+1)
	# 		for character in range(old_current,current):
	# 			left_count += (tree_string[character] == '(')
	# 			right_count += (tree_string[character] == ')')
	# 		print(current)
	# 		temp1 = tree_string[current] == ')'
	# 		temp2 = tree_string[current+1] != ','
	# 		temp3 = left_count == right_count
	# 		print(str(temp1) + ',' + str(temp2) + ',' + str(temp3))
	# 		# check = tree_string[current] == ')' and tree_string[current+1] != ',' and left_count == right_count
	# 	left_count = 1
	# 	right_count = 0
	# 	start_MRCA_2 = current + 1
	# 	end_MRCA_2 = tree_string.find(':', start_MRCA_2)
	# 	MRCA_2 = tree_string[start_MRCA_2:end_MRCA_2]
	# 	print(str(start_MRCA_1) + ',' + str(start_MRCA_2))
	# 	print(tree_string[start_MRCA_1:end_MRCA_1])

	# 	s1_tree_location = end_MRCA_1
	# 	s2_tree_location = end_MRCA_2
	# 	i += 1

	# print('distance_1 = ' + str(distance_1) + '\ndistance_2 = ' + str(distance_2))


# tree_string = '(B57:0.5,B66:0.5):0.5,(B4:0.5,((B7:0.5,B9:0.5):0.5):0.5,(B6:0.5,B15:0.5):0.5):0.5):0.5):0.5):0.5,(((((B83:0.5,B45:0.5):0.5,B21:0.5):0.5,B36:0.5):0.5,B84:0.5):0.5,(B89:0.5,B99:0.5):0.5):0.5):0.5):0.5):0.5,B58:0.5):0.0;'
# tree_string = '(B1:0.5,(B2:0.5,(B3:0.5,B4:0.5):0.5):0.5):0.0;'
# phylogeny = pyvolve.read_tree(tree = tree_string)
# pyvolve.print_tree(phylogeny)
# find_distances_to_MRCA(tree_string, 'B1', 'B2')