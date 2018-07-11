import pyvolve

phylogeny = pyvolve.read_tree(tree = "(t4:0.785,(t3:0.380,(t2:0.806,(t5:0.612,t1:0.660):0.762):0.921):0.207);")

pyvolve.print_tree(phylogeny)

custom_mu = {"AC":0.5, "AG":0.25, "AT":1.23, "CG":0.55, "CT":1.22, "GT":0.47}

freqs = [0.25,0.25,0.25,0.25]

# nuc_model = pyvolve.Model( "nucleotide", {"mu":custom_mu, "state_freqs":freqs} )
nuc_model = pyvolve.Model( "nucleotide", {"kappa":1.86836732388, "state_freqs":freqs} )

my_model = pyvolve.Model("nucleotide")

my_partition = pyvolve.Partition(models = nuc_model, size = 1000000) # 1303140)

my_evolver = pyvolve.Evolver(partitions = my_partition, tree = phylogeny)
my_evolver() 

pyvolve.print_tree(phylogeny)