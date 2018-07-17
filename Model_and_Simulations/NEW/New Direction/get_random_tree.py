
from completely_new_thing import read_in_strains
from completely_new_thing import genome_length
from completely_new_thing import get_min_m
from completely_new_thing import scale_newick_format_tree
from simulations import generate_ancestor
from process_genomes import pi_value
from process_genomes import theta_value
import pyvolve

def get_random_tree(species, filename, tree_string, L, kappa, iteration, increment):

	strains = read_in_strains(filename)
	# L = genome_length(strains)
	min_m = get_min_m(strains, L)
	scaled_tree_string = tree_string # scale_newick_format_tree(strains, L, min_m, tree_string, increment)

	phylogeny = pyvolve.read_tree(tree = tree_string)
	pyvolve.print_tree(phylogeny)

	freqs = [0.25,0.25,0.25,0.25]

	nuc_model = pyvolve.Model('nucleotide', {'kappa':kappa, 'state_freqs':freqs})

	ancestor = generate_ancestor(L)
	# print(ancestor)

	my_partition = pyvolve.Partition(models = nuc_model, root_sequence = ancestor)

	my_evolver = pyvolve.Evolver(partitions = my_partition, tree = phylogeny)
	my_evolver(ratefile = None, infofile = None, seqfile = "simulated_alignment_" + str(species[:-1]) + "_universal_" + str(iteration + 1) + ".fasta" )
	# my_evolver() 
	# my_evolver(write_anc = True)
	simulated_strains = my_evolver.get_sequences()
	# strains = my_evolver.get_sequences(anc = True)
	# strain_names = list(strains.keys())
	pi = pi_value(simulated_strains)
	theta = theta_value(simulated_strains)

	# print('pi: ' + str(pi))
	# print('theta: ' + str(theta))

	return {'pi': pi, 'theta': theta}


	