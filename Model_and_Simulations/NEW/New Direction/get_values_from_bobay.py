
import csv
import glob
import os
from process_genomes import read_in_strains
from process_genomes import id_matrix
from process_genomes import genome_length
from process_genomes import species_size
from process_genomes import pi_value
from process_genomes import theta_value

path = 'C:/Users/Owner/Documents/UNCG REU/Project/BIGG DATA/Useful Data/Concatenates, Trees, Homoplasies/Aaayyy Complete 3/'
species_names = []
ns = []
Ls = []
kappas = []
pis = []
thetas = []
h_cores = []
h_universals = []
species = os.listdir(path)
for x in range(len(species)):
	species_name = species[x]
	print(species_name)
	folder_path = os.path.join(path, species_name)

	full_path = os.path.join(folder_path, 'concat_universal.fa')
	if os.path.exists(full_path):
		concat_core_file = open((full_path), 'r') # open('concat_core.fa', 'r')
		concat_core_file = list(concat_core_file)
		strains = read_in_strains(full_path)
		strain_names = list(strains.keys())
		name = full_path[len(path):len(full_path)-3]
		print(name)
		# identity_matrix = id_matrix(strains)
		# shape = identity_matrix.shape
		# with open(('ID_Matrix_' + species_name + '.csv'), 'w', newline = '') as f: 
		# 	writer = csv.writer(f)
		# 	writer.writerow([species_name])
		# 	header = ['']
		# 	header.extend(strain_names)
		# 	writer.writerow(header)
		# 	for row in range(shape[0]):
		# 		write_row = [strain_names[row]]
		# 		write_row.extend(identity_matrix[row])
		# 		writer.writerow(write_row)
		# 	writer.writerow([])
		print(name)
		L = genome_length(strains)
		n = species_size(strains)
		pi = pi_value(strains)
		theta = theta_value(strains)

	else: 
		L = 'N/A'
		n = 'N/A'

	# if (os.path.join(full_path, 'concat_universal.fa'):
	# 	concat_universal_file = open(os.path.join(full_path, 'concat_universal.fa'), 'r')
	# 	concat_universal_file = list(concat_universal_file)

	# print(full_path)
	# concat_core_file = open(os.path.join(full_path, 'concat_core.fa'), 'r') # open('concat_core.fa', 'r')
	# concat_universal_file = open(os.path.join(full_path, 'concat_universal.fa'), 'r')
	full_path = os.path.join(folder_path, 'kappa.txt')
	if os.path.exists(full_path):
		kappa_file = open(full_path, 'r')
		kappa_file = list(kappa_file)
		kappa = kappa_file[0]
	else:
		kappa = 'N/A'

	full_path = os.path.join(folder_path, 'all_core.txt')
	if os.path.exists(full_path):
		core_file = open(full_path, 'r')
		core_file = list(core_file)
		core_values = core_file[0].split('\t')
		h_core = core_values[1]
	else:
		h_core = 'N/A'

	full_path = os.path.join(folder_path, 'all_universal.txt')
	if os.path.exists(full_path):
		universal_file = open(full_path, 'r')
		universal_file = list(universal_file)
		universal_values = universal_file[0].split('\t')
		h_universal = universal_values[1]
	else: 
		h_universal = 'N/A'

	species_names.append(species_name)
	ns.append(n)
	Ls.append(L)
	pis.append(pi)
	thetas.append(theta)
	kappas.append(kappa)
	h_cores.append(h_core)
	h_universals.append(h_universal)


with open(('species_values_universal2.csv'), 'w', newline = '') as f:
	writer = csv.writer(f)
	writer.writerow(['Species', 'n', 'L', 'pi', 'theta', 'kappa', 'h_core', 'h_universal'])
	data = [species_names, ns, Ls, pis, thetas, kappas, h_cores, h_universals]
	data = zip(*data)
	writer.writerows(data)
		# concat_core_file = list(concat_core_file)
		# concat_universal_file = list(concat_universal_file)
		# kappa_file = list(kappa_file)
		# core_file = list(core_file)
		# universal_file = list(universal_file)
		
		# kappa = kappa_file[0]
		# print('kappa = ' + str(kappa))

		# core_values = core_file[0].split('\t')
		# h_core = core_values[1]
		# print('h_core = ' + str(h_core))

		# universal_values = universal_file[0].split('\t')
		# h_universal = universal_values[1]
		# print('h_universal = ' + str(h_universal))

	# writer.writerow([species_name, n, L, kappa, h_core, h_universal])