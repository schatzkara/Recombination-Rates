
# filename = 

def read_in_strains(filename):
	f = open(filename, 'r')
	f = list(f)
	strains = {} # dictionary to hold all the strains (key: strain name, value: strain genome)

	for line in f: 
		if(line[0] == '>'): # separates out the strain names
			key = line.strip('\n')
			strains[key] = ''
		else:
			strains[key] += line.strip('\n') # concatenates all the lines of genome
	values = list(strains.values())
	length = len(values[0])
	for strain in strains.values(): # makes sure that the genome length of each strain is uniform
		# print(len(strain))
		if len(strain) != length:
			print('ABORT: THE LENGTHS ARE NOT THE SAME')
	return strains

def genome_length(species):
	strains = list(species.values())
	length = len(strains[0])
	for i in range(1, len(strains)):
		if len(strains[i]) != length:
			raise ValueError('The strains do not all have the same genome length')
	return length

def get_min_m(strains, L):
	sequences = list(strains.values())

	alleles = L*[None]
	for site in range(L):
		alleles[site] = []

	for site in range(L):
		for sequence in sequences:
			if sequence[site] not in alleles[site]:
				alleles[site].append(sequence[site])

	min_m = 0
	for site in range(L):
		min_m += (len(alleles[site]) - 1)

	return min_m

def scale(branch_length, total_branch_length, L, min_m, increment):
	scaled_branch_length = str(float(branch_length) * (int(min_m) + increment) / (float(total_branch_length) * int(L)))
	# print(scaled_branch_length)
	return scaled_branch_length
	# return min_m / L

def scale_newick_format_tree(strains, L, min_m, tree_string, increment):
	l = len(tree_string)
	total_branch_length = 0
	# print(l)
	branch_lengths = []
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
			branch_lengths.append([start, end])
			# print(tree_string[start:end])
			total_branch_length += float(tree_string[start:end])
			# print(tree_string[start:end])
			current = end
			# print('current = ' + str(current))
	scaled_tree_string = tree_string
	for branch in branch_lengths:
		branch_length = tree_string[branch[0]:branch[1]]
		# print(length)
		scaled_tree_string = scaled_tree_string.replace(branch_length, scale(branch_length, total_branch_length, L, min_m, increment))

	return scaled_tree_string


	# str.index(str, beg = 0 end = len(string))

	# str.replace(old, new[, count])
	# Return a copy of the string with all occurrences of substring old replaced by new. If the optional argument count is given, only the first count occurrences are replaced.

# strains = {1:'AAAA', 2:'TTTT', 3:'GGGG', 4:'CCCC'}
# L = 4
# min_m = get_min_m(strains, L)
# # print(min_m)
# tree_string = '((t1:0.5,t2:0.5)i1:0.5,(t3:0.5,t4:0.5)i2:0.5)root;'

# print(scale_newick_format_tree(strains, L, min_m, tree_string))