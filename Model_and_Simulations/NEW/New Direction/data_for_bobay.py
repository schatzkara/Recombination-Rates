
import os
import csv
from process_genomes import read_in_strains
from process_genomes import read_in_reduced_strains
from process_genomes import genome_length
from process_genomes import species_size
from process_genomes import pi_value
from process_genomes import theta_value
from completely_new_thing import get_min_m
from completely_new_thing import get_max_m
from completely_new_thing import scale_branch_lengths

species = ['Bacillus_anthracis']

for s in range(len(species)): 
	sp = species[s]
	print(sp)
	species_alignment = '/mnt/c/Users/Owner/Documents/UNCG/Project/BIGG_DATA/Useful_Data/Concatenates,Trees,Homoplasies/' + sp + '/concat_universal.fa'
	ancestral_alignment = 'RAxML_marginalAncestralStates.anc'
	kappa_file = '/mnt/c/Users/Owner/Documents/UNCG/Project/BIGG_DATA/Useful_Data/Concatenates,Trees,Homoplasies/' + sp + '/kappa.txt'
	output_file = 'scaling_trials_' + sp + '_universal'
	reduced_species_alignment = '/mnt/c/Users/Owner/Documents/UNCG/Project/standard-RAxML/done_species/'  + sp + 'concat_universal.fa.reduced'
	raxml_path = '/mnt/c/Users/Owner/Documents/UNCG/Project/standard-RAxML/done_species/' + sp
	tree_file = 'RAxML_bestTree.tree'
	rooted_tree_file = 'RAxML_rootedTree.root'
	ancestral_tree_file = 'RAxML_nodeLabelledRootedTree.anc'
	reduced = os.path.exists(reduced_species_alignment)
	if not reduced:
		strains = read_in_strains(species_alignment)
	else:
		strains = read_in_reduced_strains(reduced_species_alignment) # dictionary with the genomes of all the strains; key = strain name, value = genome

	L = genome_length(strains) # number of base pairs in the genome
	n = species_size(strains) # number of extant strains
	strain_names = list(strains.keys()) # list of all the extant strain names

	tree_file = open((os.path.join(raxml_path, rooted_tree_file)), 'r')
	rooted_tree_string = list(tree_file)[0]
	# tree_file = open((os.path.join(raxml_path, ancestral_tree_file)), 'r')
	# ancestral_tree_string = list(tree_file)[0]
	# internal_nodes = get_internal_nodes(os.path.join(raxml_path, ancestral_alignment))
	# internal_nodes,ancestral_tree_string = rename_ancestors(internal_nodes, strain_names, ancestral_tree_string)
	# all_nodes = {}
	# for key in strains.keys():
		# all_nodes[key] = strains[key]
	# for key in internal_nodes.keys():
		# all_nodes[key] = internal_nodes[key]

	# all_node_names = list(all_nodes.keys())

	total_pairs = int((n*(n-1)) / 2) # the total number of strain pairs that will be compared
	pi = pi_value(strains)
	theta = theta_value(strains) # proportion of the genome that is polymorphic
	
	# complete_tree_string = merge_trees(rooted_tree_string, ancestral_tree_string)

	kappa_file = open(kappa_file, 'r')
	kappa = float(list(kappa_file)[0])
	min_m = get_min_m(strains, L) # minimum number of mutations that could account for all the polymorphisms in the species
	max_m = get_max_m(strains, L, rooted_tree_string)

	accurate_tree, trials = scale_branch_lengths(L, rooted_tree_string, min_m, max_m, pi, theta, kappa, 1)

	with open(output_file + '.csv', 'w', newline='') as f:
		writer = csv.writer(f)
		writer.writerow(sp)
		writer.writerow(['pi','theta'])
		writer.writerow([pi, theta])
		writer.writerow('')
		writer.writerow(['desired_m','max_m','scaling_factor','pi_sim','theta_sim'])
		writer.writerows(trials)