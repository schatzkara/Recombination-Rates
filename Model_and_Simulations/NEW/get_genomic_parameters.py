
import csv
from process_genomes import read_in_strains
from process_genomes import pi_value
from process_genomes import theta_value
from process_genomes import nucleotide_composition

# runs the functions to get pi, theta, GC% average, and GC% standard deviation for each species and write them into a .csv file
# time complexity: O(n^4), where n is the length of the strains
path = 'C:/Users/Owner/Documents/UNCG REU/Project/Recombination-Rates/concatenates' # path where the .fa files are located 
with open(('species_params3.csv'), 'w', newline = '') as f: 
	writer = csv.writer(f)
	writer.writerow(['species', 'pi', 'theta', 'GC%']) # column headers
	for filename in glob.glob(os.path.join(path, '*.fa')): # finds the values for each species
		species = read_in_strains(filename)
		name = filename.strip('C:/Users/Owner/Documents/UNCG REU/Project/Recombination-Rates/concatenates').strip('/concat_') # strips off everything but the actual species name
		print(name)
		pi = str(pi_value(species))
		theta = str(theta_value(species))
		# GC_comp = list(nucleotide_composition(species))
		GC_average = nucleotide_composition(species)
		# GC_average = GC_comp[0]
		# GC_stdev = GC_comp[1]


		writer.writerow([name, pi, theta, GC_average])

		# print(filename)
		# print('pi = ' + str(pi))
		# print('theta = ' + str(theta))
		# print('GC% = ' + str(GC_comp))
		print('\n')
