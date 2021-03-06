#! python3

# script to estimate Pi, Theta, and the nucleotide composition (GC%) of each clonal species

# import numpy as np
import os 
import glob
import csv
import numpy as np

# function to read in the genomes and separate the strains from the .fa files
# params: 
# 	filename (str) = name of .fa file
# return: dictionary with all the strains (key: strain name, value: strain genome)
# time complexity: O(n), where n is the number of lines in the file
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

def species_size(species):
	strains = list(species.keys())
	return len(strains)
# this function finds the value of pi, which is the average proportion of site differences b/w each 2 genomes
# units: average number of differences b/w each 2 genomes / length
# params: 
# 	species (dict) = where the values are the genomes of different strains
# return: float that equals the pi value of the species
# time complexity: O(n^3), where n is the length of genomes
def pi_value(species):
	diff_prop = [] # list of the proportion of site differences b/w each 2 genomes
	differences = 0 # counter for the number of differences
	strains = list(species.values())
	for i in range(len(strains)-1): # allows each strain except the last to be strain1
		strain1 = strains[i]
		for j in range(i+1,len(strains)): # allows each strain following strain1 to be strain2
			strain2 = strains[j]
			for k in range(len(strain1)):
				if strain1[k] != strain2[k] and strain1[k] != '-' and strain2[k] !='-': # ignores any genome sites that are '-'
					differences+=1
			diff_prop.append(differences/len(strain1))
			differences = 0 # resets the difference counter before moving on to a new pair of strains
	# pi = np.mean(diff_prop)
	pi = sum(diff_prop)/len(diff_prop)
	# print(pi)
	return pi

# this function finds the value of theta, which is the proportion of of sites that are polymorphic (have mutations) 
# units: # of polymorphic sites / length
# params: 
# 	species (dict) = where the values are the genomes of different strains
# return: float that equals the theta value of the species
# time complexity: O(n^2), where n is the length of genomes
def theta_value(species):
	diff_sites = [] # list of the sites that are polymorphic
	strains = list(species.values())
	strain1 = strains[0] 
	for k in range(len(strain1)): # loops through each genome site
		for i in range(1,len(strains)): # compares each other strain to the first in the species dict
			strain2 = strains[i]
			if strain1[k] != strain2[k] and strain1[k] != '-' and strain2[k] != '-' and k not in diff_sites: # ignores any genome sites that are '-' or are already marked as polymorphic
				diff_sites.append(k)
	theta = len(diff_sites)/len(strains[0])
	# print(theta)
	return theta

# this function finds the nucleotide composition (GC%) of the entire species
# units: # of GC sites / length
# params: 
# 	species (dict) = where the values are the genomes of different strains
# return: list that contains the average GC% of all the strains and the standard deviation
# time compleity: O(n^2), where n is the length of the strains
def nucleotide_composition(species):
	GC_comp = [] # list of the GC% of each strain
	strains = list(species.values())
	# values = [] # list containing GC% average and standard deviation over the entire species
	for strain in strains:
		GC = 0 # counter for the number of Gs and Cs in the strain
		for k in range(len(strain)):
			if(strain[k] == 'G' or strain[k] == 'C'):
				GC+=1
		GC_comp.append(GC/len(strain))
	# print('The average GC composition is ' + str(np.mean(GC_comp)))
	# print('The STD of GC composition is ' + str(np.std(GC_comp)))
	# return (str(np.mean(GC_comp)) + ' with STD ' + str(np.std(GC_comp)))
	# return np.mean(GC_comp)
	# values.append(np.mean(GC_comp))
	# values.append(np.std(GC_comp))
	# print(values)
	# return values
	GC = sum(GC_comp)/len(GC_comp)
	print(GC)
	return GC


def id_matrix(species): # Given ordered list of String sequences seqList
	strain_names = list(species.keys())
	n = len(strain_names)
	ID = np.empty([n,n], dtype = np.float, order='C')
	for s1 in range(n):
		strain1 = species[strain_names[s1]]
		ID[s1,s1] = 1
		for s2 in range(s1+1,n):
			strain2 = species[strain_names[s2]]
			identity = calc_id(strain1, strain2)
			ID[s1,s2] = identity
			ID[s2,s1] = identity
	return ID

def calc_id (s1, s2):
	num_diff = 0
	num_same = 0
	if(len(s1) != len(s2)):
		raise(ValueError('Strand Lengths not equal'))
	else:
		for i in range(len(s1)):
			if s1[i] != '-':
				num_diff += (s1[i]!=s2[i]) # Uses boolean as int, +=1 if true, +=0 if false
				num_same += (s1[i]==s2[i])
	# if (num_diff + num_same) != len(s1):
		# raise(ValueError('numDif + numSame != length'))
	return num_same/len(s1)
	
# filename or species	
def average_id(**keyword_parameters):
	if 'filename' in keyword_parameters:
		filename = keyword_parameters['filename']
		with open(filename) as d:
			reader = csv.reader(d)
			next(reader)
			next(reader)
			n = 0
			total = 0
			i = 2
			for row in reader:
				n = len(row)-1
				for j in range(i,len(row)):
					total += row[j]
				i += 1
		average_id = total/((n*(n-1))/2)
	if 'species' in keyword_parameters:
		species = keyword_parameters['species']
		strain_names = list(species.keys())
		n = len(strain_names)
		total = 0
		# ID = np.empty([n,n], dtype = np.float, order='C')
		for s1 in range(n):
			strain1 = species[strain_names[s1]]
			# ID[s1,s1] = 1
			for s2 in range(s1+1,n):
				strain2 = species[strain_names[s2]]
				identity = calc_id(strain1, strain2)
				total += identity
				# ID[s1,s2] = identity
				# ID[s2,s1] = identity
		average_id = total/((n*(n-1))/2)
	return average_id

# # runs the functions to get pi, theta, GC% average, and GC% standard deviation for each species and write them into a .csv file
# # time complexity: O(n^4), where n is the length of the strains
# path = 'C:/Users/Owner/Documents/UNCG REU/Project/Recombination-Rates/concatenates' # path where the .fa files are located 
# with open(('species_params3.csv'), 'w', newline = '') as f: 
# 	writer = csv.writer(f)
# 	writer.writerow(['species', 'pi', 'theta', 'GC%']) # column headers
# 	for filename in glob.glob(os.path.join(path, '*.fa')): # finds the values for each species
# 		species = read_in_strains(filename)
# 		name = filename.strip('C:/Users/Owner/Documents/UNCG REU/Project/Recombination-Rates/concatenates').strip('/concat_') # strips off everything but the actual species name
# 		print(name)
# 		pi = str(pi_value(species))
# 		theta = str(theta_value(species))
# 		# GC_comp = list(nucleotide_composition(species))
# 		GC_average = nucleotide_composition(species)
# 		# GC_average = GC_comp[0]
# 		# GC_stdev = GC_comp[1]


# 		writer.writerow([name, pi, theta, GC_average])

# 		# print(filename)
# 		# print('pi = ' + str(pi))
# 		# print('theta = ' + str(theta))
# 		# print('GC% = ' + str(GC_comp))
# 		print('\n')
