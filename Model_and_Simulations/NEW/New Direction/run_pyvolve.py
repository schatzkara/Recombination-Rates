
import csv
from get_random_tree import get_random_tree


L = 46935 # length to simulate
kappa = 1.86836732388 # kappa value for the species
increments = 1
iterations = 100 # number of iterations to simulate
# tree in Newick format (from RAxML)
tree_string = '(B57:0.00000536707125664551,((((B83:0.00002318173951070520,B45:0.00003365367407355436):0.00001970126751875890,(((B84:0.00000948118264665464,B99:0.00001914943883007716):0.00001143909786772645,(((B52:0.00005185110847249898,B89:0.00001633969855456498):0.00004537804571909241,B21:0.00011406645637562866):0.00001953691402093377,(B72:0.00008420804257297685,((B51:0.00008978912752602243,B7:0.00002287080714303825):0.00029324771986313119,((B34:0.00002730728276759067,(B9:0.00003316646850023659,B29:0.00002086920547771337):0.00008771969195870552):0.00002492360304559525,(B6:0.00000881808638051164,B15:0.00006370050086096601):0.00001778245864060021):0.00012946764031754985):0.00011766209279446951):0.00001247173836728674):0.00000100000050002909):0.00000844349545619994,B36:0.00002509526808151878):0.00000667150366752081):0.00009163824755429087,B4:0.00002332293862512587):0.00000616528362070958,(B66:0.00003411246880153140,B64:0.00003762063742532716):0.00000653195343569950):0.00001043697034630331,B58:0.00004497951428308054):0.0;'
# complete path for the .fa file that contains the strains' genomes
filename = 'C:/Users/Owner/Documents/UNCG REU/Project/BIGG DATA/Useful Data/Concatenates, Trees, Homoplasies/Aayyy Clonal/Bacillus_anthracis/concat_universal.fa'

pis = iterations*[None] # list that will be populated with the simulated pi values; index = iteration - 1
thetas = iterations*[None] # list that will be populated with the simulated theta values; index = iteration - 1
for i in range(iterations):
	values = get_random_tree(filename, tree_string, L, kappa, i, 0) # 
	# pi = values['pi']
	# theta = values['theta']
	pis[i] = values['pi']
	thetas[i] = values['theta']
	# print('increment: ' + str(i))
	# print('pi_sim: ' + str(pi))
	# print('theta_sim: ' + str(theta))
	print('pis: ' + str(pis))
	print('thetas: ' + str(thetas))
print('pi_real: ' + str(0.000269877))
print('theta_real: ' + str(0.001704485))

