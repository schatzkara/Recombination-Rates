# script to estimate Pi, Theta and the nucleotide composition (GC%) of each clonal species

import numpy as np
import os 
import glob
import csv

# time complexity: O(n), where n is the number of lines in the file
def read_in_strains(filename):
	f = open(filename, 'r')
	f = list(f)
	# print(f)
	strains = {}

	for line in f: 
		if(line[0] == '>'):
			key = line.strip('\n')
			strains[key] = ''
		else:
			strains[key] += line.strip('\n')
	values = list(strains.values())
	length = len(values[0])
	for strain in strains.values():
		# print(len(strain))
		if len(strain) != length:
			print('ABORT: THE LENGTHS ARE NOT THE SAME')
	return strains

# this function finds the value of pi, which is the average proportion of site differences b/w each 2 genomes
# units: average number of differences b/w each 2 genomes / length
# time complexity: O(n^3), where n is the bigger of the number of strains and the length of the strains; technically it's n^2 * L
def pi_value(species):
	diff_prop = []
	differences = 0
	strains = list(species.values())
	for i in range(len(strains)-1):
		strain1 = strains[i]
		for j in range(i+1,len(strains)):
			strain2 = strains[j]
			for k in range(len(strain1)):
				# print(k)
				if strain1[k] != strain2[k] and strain1[k] != '-' and strain2[k] !='-':
					differences+=1
			diff_prop.append(differences/len(strain1))
			differences = 0
	pi = np.mean(diff_prop)
	return pi

# this function finds the value of theta, which is the proportion of of sites that are polymorphic (have mutations) 
# units: # of polymorphic sites / length
# time complexity: O(n^2), where n is the bigger of the length of the strains and the number of strains; technically it's L * n
def theta_value(species):
	diff_sites = []
	strains = list(species.values())
	strain1 = strains[0]
	for k in range(len(strain1)):
		for i in range(1,len(strains)):
			strain2 = strains[i]
			if strain1[k] != strain2[k] and strain1[k] != '-' and strain2[k] != '-' and k not in diff_sites:
				diff_sites.append(k)
	theta = len(diff_sites)/len(strains[0])
	return theta

# this function finds the nucleotide composition (GC%)
# units: # of GC sties / length
# time compleity: O(n^2), where n is the bigger of the number of strains and the length of the strains; technically, it's n * L
def nucleotide_composition(species):
	GC_comp = []
	strains = list(species.values())
	for strain in strains:
		GC = 0
		for k in range(len(strain)):
			if(strain[k] == 'G' or strain[k] == 'C'):
				GC+=1
		GC_comp.append(GC/len(strain))
	# print('The average GC composition is ' + str(np.mean(GC_comp)))
	# print('The STD of GC composition is ' + str(np.std(GC_comp)))
	# return (str(np.mean(GC_comp)) + ' with STD ' + str(np.std(GC_comp)))
	# return np.mean(GC_comp)
	return [str(np.mean(GC_comp)), str(np.std(GC_comp))]

# time complexity: O(n^4), where n is the biggest of the number of files, the number of lines in the files, the number of strains, and the length of the strains; technically, it's n^4 + 2n^3 + n^2
path = 'C:/Users/Owner/Documents/UNCG REU/Project/Recombination-Rates/concatenates'
with open(('species_params.csv'), 'w', newline = '') as f:
	writer = csv.writer(f)
	writer.writerow(['species', 'pi', 'theta', 'GC%', 'STDev'])
	for filename in glob.glob(os.path.join(path, '*.fa')):
		species = read_in_strains(filename)
		name = filename.strip('C:/Users/Owner/Documents/UNCG REU/Project/Recombination-Rates/concatenates')
		pi = str(pi_value(species))
		theta = str(theta_value(species))
		GC_comp = str(nucleotide_composition(species))
		GC_average = GC_comp[0]
		GC_stdev = GC_comp[1]


		writer.writerow([name, pi, theta, GC_average, GC_stdev])

		print(filename)
		print('pi = ' + str(pi))
		print('theta = ' + str(theta))
		print('GC% = ' + str(GC_comp))
		print('\n')

# for filename in glob.glob(os.path.join(path, '*.fa')):
# 		species = read_in_strains(filename)

# 		print(filename)
# 		print('pi = ' + str(pi_value(species)))
# 		print('theta = ' + str(theta_value(species)))
# 		print('GC% = ' + str(nucleotide_composition(species)))
# 		print('\n')


# print('pi = ' + str(pi_value(read_in_strains(0))))
# print('theta = ' + str(theta_value(read_in_strains(0))))
# print('GC% = ' + str(nucleotide_composition(read_in_strains(0))))