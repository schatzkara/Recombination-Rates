
import numpy
from process_genomes import read_in_strains
from process_genomes import species_size
from process_genomes import genome_length
from process_genomes import theta_value
from completely_new_thing import get_min_m
from completely_new_thing import scale_newick_format_tree

SHARED = np.empty([n,n], dtype = np.float, order='C')
ANCESTRAL = np.empty([n,n], dtype = np.float, order='C')
RECOMBINANT = np.empty([n,n], dtype = np.float, order='C')
CONVERGENT = np.empty([n,n], dtype = np.float, order='C')

path = 'C:/Users/Owner/Documents/UNCG REU/Project/Data from Bobay/Aayyy Clonal/Bacillus_anthracis/concat_core.fa'
kappa = 1.86836732388
tree_string = '(B57:0.00000536707125664551,((((B83:0.00002318173951070520,B45:0.00003365367407355436):0.00001970126751875890,(((B84:0.00000948118264665464,B99:0.00001914943883007716):0.00001143909786772645,(((B52:0.00005185110847249898,B89:0.00001633969855456498):0.00004537804571909241,B21:0.00011406645637562866):0.00001953691402093377,(B72:0.00008420804257297685,((B51:0.00008978912752602243,B7:0.00002287080714303825):0.00029324771986313119,((B34:0.00002730728276759067,(B9:0.00003316646850023659,B29:0.00002086920547771337):0.00008771969195870552):0.00002492360304559525,(B6:0.00000881808638051164,B15:0.00006370050086096601):0.00001778245864060021):0.00012946764031754985):0.00011766209279446951):0.00001247173836728674):0.00000100000050002909):0.00000844349545619994,B36:0.00002509526808151878):0.00000667150366752081):0.00009163824755429087,B4:0.00002332293862512587):0.00000616528362070958,(B66:0.00003411246880153140,B64:0.00003762063742532716):0.00000653195343569950):0.00001043697034630331,B58:0.00004497951428308054):0.0;'

ancestor = ??????????
strains = read_in_strains(path)
strain_names = list(strains.keys())
n = species_size(strains)
L = genome_length(strains)
theta = theta_value(strains)
mu = theta/(2*n)
min_m = get_min_m(strains, L)
scaled_tree_string = scale_newick_format_tree(strains, L, min_m, tree_string)

for s1 in range(n):
	strain1 = strains[strain_names[s1]]
	for s2 in range(s1,n):
		strain2 = strains[strain_names[s2]]
		s,a,r,c = 0,0,0,0
		for site in range(L):
			if strain1[site] == strain2[site]:
				s += 1
				if strain1[site] == ancestor[site]:
					a += 1
		s1_tree_location = scaled_tree_string.find(strain_names[s1])
		s2_tree_location = scaled_tree_string.find(strain_names[s2])

		start_length_1 = scaled_tree_string.find(':', s1_tree_location) + 1
		x1 = scaled_tree_string.find(',', start_length_1)
		y1 = scaled_tree_string.find(')', start_length_1)
		if x1 == -1:
			end_length_1 = y1
		elif y1 == -1:
			end_length_1 = x1
		else:
			end_length_1 = min(x1,y1)

		start_length_2 = scaled_tree_string.find(':', s2_tree_location) + 1
		x2 = scaled_tree_string.find(',', start_length_2)
		y2 = scaled_tree_string.find(')', start_length_2)
		if x2 == -1:
			end_length_2 = y2
		elif y2 == -1
			end_length_2 = x2
		else:
			end_length_2 = min(x2,y2)

		length_1 = scaled_tree_string[start_length_1:end_length_1]
		length_2 = scaled_tree_string[start_length_2:end_length_2]

		m_1 = length_1 * L
		m_2 = length_2 * L


		c = expected_c_given_ms(L,m_1,m_2,generations,kappa,phi)


		SHARED[s1,s2] = s
		ANCESTRAL[s1,s2] = a
		CONVERGENT[s1,s2] = c
		RECOMBINANT = s - a - c




