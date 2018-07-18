
from process_genomes import read_in_strains
from process_genomes import genome_length
from completely_new_thing import get_min_m
from completely_new_thing import scale_newick_format_tree
from simulations import generate_ancestor
from process_genomes import pi_value
from process_genomes import theta_value
import pyvolve

def get_random_tree(species, scaled_tree_string, kappa, iteration):
	# strains = read_in_strains(filename)
	# L = genome_length(strains)
	# min_m = get_min_m(strains, L)
	# max_m = get_max_m(strains, L, tree_string)
	# pis = []
	# thetas = []

	# scaled_trees = []

	# for x in range(min_m,max_m+1):
	# 	scaled_tree_string = scale_newick_format_tree(strains, L, x, tree_string, increment)
	# 	scaled_trees.append(scaled_tree_string)

	# for tree in scaled_trees:
	phylogeny = pyvolve.read_tree(tree = scaled_tree_string)
	print('read in the tree')
	pyvolve.print_tree(phylogeny)

	freqs = [0.25,0.25,0.25,0.25]

	nuc_model = pyvolve.Model('nucleotide', {'kappa':kappa, 'state_freqs':freqs})

	ancestor = generate_ancestor(L)
	print('generated an ancestor')
# 	# print(ancestor)

	my_partition = pyvolve.Partition(models = nuc_model, root_sequence = ancestor)

	my_evolver = pyvolve.Evolver(partitions = my_partition, tree = phylogeny)
	my_evolver(ratefile = None, infofile = None, seqfile = "simulated_alignment_" + str(species[:-1]) + "_universal_" + str(iteration + 1) + ".fasta" )
# 	# my_evolver() 
	print('evolved the sequences')
# 	# my_evolver(write_anc = True)
	simulated_strains = my_evolver.get_sequences()
# 	# strains = my_evolver.get_sequences(anc = True)
# 	# strain_names = list(strains.keys())
	pi = pi_value(simulated_strains)
	theta = theta_value(simulated_strains)
# 	pis.append(pi)
# 	thetas.append(theta)

	# # print('pi: ' + str(pi))
	# # print('theta: ' + str(theta))

	# return {'pi': pis, 'theta': thetas}

	return pi,theta

	