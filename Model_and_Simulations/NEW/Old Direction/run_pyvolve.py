
import csv
from get_random_tree import get_random_tree


L = 1000
kappa = 1.86836732388
iterations = 1000
tree_string = '(B57:0.00000536707125664551,((((B83:0.00002318173951070520,B45:0.00003365367407355436):0.00001970126751875890,(((B84:0.00000948118264665464,B99:0.00001914943883007716):0.00001143909786772645,(((B52:0.00005185110847249898,B89:0.00001633969855456498):0.00004537804571909241,B21:0.00011406645637562866):0.00001953691402093377,(B72:0.00008420804257297685,((B51:0.00008978912752602243,B7:0.00002287080714303825):0.00029324771986313119,((B34:0.00002730728276759067,(B9:0.00003316646850023659,B29:0.00002086920547771337):0.00008771969195870552):0.00002492360304559525,(B6:0.00000881808638051164,B15:0.00006370050086096601):0.00001778245864060021):0.00012946764031754985):0.00011766209279446951):0.00001247173836728674):0.00000100000050002909):0.00000844349545619994,B36:0.00002509526808151878):0.00000667150366752081):0.00009163824755429087,B4:0.00002332293862512587):0.00000616528362070958,(B66:0.00003411246880153140,B64:0.00003762063742532716):0.00000653195343569950):0.00001043697034630331,B58:0.00004497951428308054):0.0;'
filename = 'C:/Users/Owner/Documents/UNCG REU/Project/Data from Bobay/Aayyy Clonal/Bacillus_anthracis/concat_core.fa'

pis = iterations*[None]
thetas = iterations*[None]
for i in range(iterations):
	values = get_random_tree(filename, tree_string, L, kappa)
	pis[i] = values['pi']
	thetas[i] = values['theta']

print('pis: ' + str(pis))
print('thetas: ' + str(thetas))
print(sum(pis))
print(sum(thetas))




# with open(('c_from_pyvolve_' + str(L) + '_' + str(kappa) + '.csv'), 'w', newline = '') as f:
# 	writer = csv.writer(f)
# 	writer.writerow(['Iteration', 'c', 'L', 'kappa'])
# 	for i in range(iterations):
# 		writer.writerow([i+1, get_c(L, kappa), L, kappa])