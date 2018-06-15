# Recombination-Rates
Program to determine the rate of recombination in bacterial species

This program does a few different tasks:
	- calculates the expected number of convergent mutations between strands of DNA
	- simulates actual mutations along strands of DNA and detects the number of convergent mutations that occurred
	- calculates certain important parameters for any given .fa file of a bacterial genome

# Part 1: Calculating the Expected Number of Convergent Mutations
	Based on a probablistic model that we developed, we can determine the number of convergent mutations that we would expect to arise between and pair of DNA sequences. Eventually this value will be used to estimate the number of recombinant events that have occurred in a species so that we can then determine the recombination rates. We have another program that can detect the number of homoplasies along a phylogenetic tree. These homoplasies have arisen exclusively due to convergent mutations and recombinant events. Thus, since we know the number of homoplasies and the expected number of covnergent mutations, we can indirectly determine the number of recombinant events. 

# Part 2: Simulating Mutations Along Strands of DNA
	This simulation takes arbitrary DNA strands and mutates them according to five inputs from the user: 
		n (int) = number of DNA strands to mutate
		L (int) = length of the DNA strands
		mu (float) = mutation rate (units: mutations per base pair per generation)
		kappa (float) = ratio of transitions to transversions
		phi (float) = probability of transversion to its complementary base pair 
	Then, it counts up and outputs the total number of convergent mutations that have occurred between all pairs of DNA strands upon which the simulation acted.

	The simulation begins by randomly generating a single 'ancestor' strand of DNA of length L and duplicating this strand until there are n identical strands. Then, it goes along each nucleotide in each strand of DNA and selects a 'new' nucleotide based on a weigthed random choice (the nucleotide is not necessarily new since the random choice can be what it was originally). The weights are as follows:
		(1- mu) = probability of the original nucleotide
		(mu * kappa)/(kappa + 1) = probability of transition (mutation to the nucleotide of the same 'class' (purines and pyrimidines))
		(mu/(kappa + 1)*phi) = probability of transversion to its complementary base pair
		(mu/(kappa + 1)*(1-phi)) = probability of transversion to the remaining nucleotide

	After all the 'mutations' have occurred for each nucleotide on each strand, the number of convergent mutations between each pair of strands is tallied and outputted.

# Part 3: Calculating Species Parameters
	The program calculates 3 important parameters for any bacterial species when the user inputs the .fa file of the genome. These parameters are:
		- pi: the average proportion of site differences b/w each 2 genomes (units: average number of differences b/w each 2 genomes / length)
		- theta: the proportion of of sites that are polymorphic (have mutations) (units: # of polymorphic sites / length)
		- GC%: percentage of the total genome that are Gs and Cs (units: # of GC sites / length)
	It outputs these to a .csv file in the same directory as the program.