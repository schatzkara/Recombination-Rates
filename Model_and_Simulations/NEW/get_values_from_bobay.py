
import csv
import glob
import os

path = 'C:/Users/Owner/Documents/UNCG REU/Project/Data from Bobay'
with open(('homoplaies.csv'), 'w', newline = '') as f:
	writer = csv.writer(f)
	writer.writerow(['Species', 'kappa', 'h_core', 'h_universal'])
	species = os.listdir(path)
	# folders = [x[0] for x in os.walk(path)]
	# for folder in glob.glob(path):
	# print(folders)
	# print(glob.glob('C:/Users/Owner/Documents/UNCG REU/Project/Data from Bobay/export_REU/*/'))
	for x in range(len(species)):
		species_name = species[x]
		# print(folders[x])
		# print(species_name)
		full_path = os.path.join(path, species_name)
		# print(full_path)
		concat_core_file = open(os.path.join(full_path, 'concat_core.fa'), 'r') # open('concat_core.fa', 'r')
		concat_universal_file = open(os.path.join(full_path, 'concat_universal.fa'), 'r')
		kappa_file = open(os.path.join(full_path, 'kappa.txt'), 'r')
		core_file = open(os.path.join(full_path, 'all_core.txt'), 'r')
		universal_file = open(os.path.join(full_path, 'all_universal.txt'), 'r')
		
		concat_core_file = list(concat_core_file)
		concat_universal_file = list(concat_universal_file)
		kappa_file = list(kappa_file)
		core_file = list(core_file)
		universal_file = list(universal_file)
		
		kappa = kappa_file[0]
		print('kappa = ' + str(kappa))

		core_values = core_file[0].split('\t')
		h_core = core_values[3]
		print('h_core = ' + str(h_core))

		universal_values = universal_file[0].split('\t')
		h_universal = universal_values[3]
		print('h_universal = ' + str(h_universal))

		writer.writerow([species_name, kappa, h_core, h_universal])



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





# path = 'C:/Users/Owner/Documents/UNCG REU/Project/concatenates/other' # path where the .fa files are located 
# with open(('species_params2.csv'), 'w', newline = '') as f: 
# 	writer = csv.writer(f)
# 	writer.writerow(['Species', 'Genome Length', 'pi', 'theta', 'GC%']) # column headers
# 	for filename in glob.glob(os.path.join(path, '*.fa')): # finds the values for each species
# 		species = read_in_strains(filename)
# 		name = filename[len(path)+1:len(filename)-3] # filename.strip('C:/Users/Owner/Documents/UNCG REU/Project/Recombination-Rates/concatenates').strip('/concat_') # strips off everything but the actual species name
# 		# name = filename.strip('C:/Users/Owner/Documents/UNCG REU/Project/Recombination-Rates/concatenates').strip('/concat_') # strips off everything but the actual species name
# 		print(name)
# 		size = species_size(species)
# 		length = genome_length(species)
# 		pi = str(pi_value(species))
# 		theta = str(theta_value(species))
# 		# GC_comp = list(nucleotide_composition(species))
# 		GC_average = nucleotide_composition(species)
# 		# GC_average = GC_comp[0]
# 		# GC_stdev = GC_comp[1]


# 		writer.writerow([name, length, size, pi, theta, GC_average])

# 		# print(filename)
# 		# print('pi = ' + str(pi))
# 		# print('theta = ' + str(theta))
# 		# print('GC% = ' + str(GC_comp))
# 		# print('\n')
# 		print(name)


# # for filename in glob.glob(os.path.join(path, '*.fa')): # finds the values for each species
# # 	species = read_in_strains(filename)
# # 	name = filename[len(path)+1:len(filename)-3] # filename.strip('C:/Users/Owner/Documents/UNCG REU/Project/Recombination-Rates/concatenates').strip('/concat_') # strips off everything but the actual species name
# # 	print(name)
# # 	print(genome_length(species))
# # 	print(species_size(species))
