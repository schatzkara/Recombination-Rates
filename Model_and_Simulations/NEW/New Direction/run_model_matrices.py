
from get_matrices import get_SCAR_matrices
from get_matrices import print_matrices
import csv

species = ['Acinetobacter_pittii', 'Aeromonas_hydrophila', 'Bacillus_anthracis', 'Bacillus_coagulans','Bacillus_licheniformis', 'Bacillus_pumilus', 'Bacillus_subtilis', 'Bordetella_pertussis'] # 'Mycobacterium_tuberculosis'

# species_alignment = '/mnt/c/Users/Owner/Documents/UNCG/Project/BIGG_DATA/Useful_Data/Concatenates,Trees,Homoplasies/' + species + '/concat_universal.fa'
# # raxml_path = '/mnt/c/Users/Owner/Documents/UNCG/Project/standard-RAxML'
# # tree_file = 'RAxML_bestTree.tree'
# # rooted_tree_file = 'RAxML_rootedTree.root'
# ancestral_alignment = 'RAxML_marginalAncestralStates.anc'
# # ancestral_tree_file = 'RAxML_nodeLabelledRootedTree.anc'
# kappa_file = '/mnt/c/Users/Owner/Documents/UNCG/Project/BIGG_DATA/Useful_Data/Concatenates,Trees,Homoplasies/' + species + '/kappa.txt'
# output_file = 'SCAR_' + species + '_universal'
mu = 0.0000000009665 # this is the average of all the known mutation rates
# 0.0000000001150 # Vibrio_cholerae
# 0.0000000001330 # Burkholderia_cenocepacia
# 0.0000000003350 # Bacillus_subtilis
# 0.0000000009665 # this is the average of all the known mutation rates
average_rates = []

for s in range(len(species)):
	species_alignment = '/mnt/c/Users/Owner/Documents/UNCG/Project/BIGG_DATA/Useful_Data/Concatenates,Trees,Homoplasies/' + species[s] + '/concat_universal.fa'
	# raxml_path = '/mnt/c/Users/Owner/Documents/UNCG/Project/standard-RAxML'
	# tree_file = 'RAxML_bestTree.tree'
	# rooted_tree_file = 'RAxML_rootedTree.root'
	ancestral_alignment = 'RAxML_marginalAncestralStates.anc'
	# ancestral_tree_file = 'RAxML_nodeLabelledRootedTree.anc'
	kappa_file = '/mnt/c/Users/Owner/Documents/UNCG/Project/BIGG_DATA/Useful_Data/Concatenates,Trees,Homoplasies/' + species[s] + '/kappa.txt'
	output_file = 'SCAR_' + species[s] + '_universal'

	print(species[s])

	SCAR = get_SCAR_matrices(species_alignment, ancestral_alignment, kappa_file, mu, species[s])
	strain_names,S,C,A,R,RATES,average_rate = SCAR['strain_names'],SCAR['Shared'], SCAR['Convergent'], SCAR['Ancestral'], SCAR['Recombinant'], SCAR['Rates'], SCAR['average']

	print_matrices(output_file, S,C,A,R,RATES, strain_names, average_rate)
	average_rates.append(average_rate)


with open('average_rates2.csv', 'w', newline = '') as f:
	writer = csv.writer(f)
	data = [species, average_rates]
	data = zip(*data)
	writer.writerows(data)
	
