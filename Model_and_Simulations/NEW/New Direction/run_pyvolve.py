
import csv
from get_random_tree import get_random_tree
from process_genomes import read_in_strains
from process_genomes import pi_value
from process_genomes import theta_value
from process_genomes import genome_length
from completely_new_thing import get_min_m
from completely_new_thing import get_max_m
from completely_new_thing import scale_branch_lengths

# num_species = 1
# L = [50586] # 50922 # 51003 # 50862 # 46935 # length to simulate
# kappa = [3.68727779895] # 2.20303254463 # 1.57113790459 # 1.53913593397 # 1.86836732388 # kappa value for the species

# pi = [0.014881669, 0.015538238, 0.029502539, 0.017970231, 0.022646838, 0.000269877, 0.030714994, 0.022653747, 0.00063736, 0.007376288, 0.033235612, 0.048060406, 0.020307055]
# theta = [0.090211703, 0.040497307, 0.12800047, 0.047206737, 0.174748575, 0.001704485, 0.175504076, 0.06053807, 0.002753762, 0.027932603, 0.106013944, 0.33979427, 0.088020967]
 # [0.004698989, 0.005666333, 0.004025765]
# increments = 11
iterations = 10 # number of iterations to simulate
# tree in Newick format (from RAxML)
# tree_string = '(B48:0.00000100000050002909,((B45:0.00022386286104843084,(B15:0.00006734345780892552,(B47:0.00010146791215414724,B21:0.00003573837930569803):0.00004371436151746642):0.00015916453666892807):0.00020127175555802429,((B11:0.00197082117122606391,B44:0.00045022166930896947):0.00025515447156947791,(((((B49:0.00001111048843933425,(B34:0.00001644704589964510,B17:0.00001443458119727544):0.00000458342822225343):0.00010435101929535647,B29:0.00008132003817758088):0.00000627978758453736,B12:0.00005050690227474338):0.00010530055326551897,(B37:0.00003327008481985708,(B14:0.00000911948037445038,(B39:0.00008431433646157074,B40:0.00002392328779324190):0.00005634460179426050):0.00002542404904296348):0.00009325017620007546):0.00000100000050002909,((B26:0.00005081173789799183,B35:0.00005504641358681077):0.00001181973065473387,(B22:0.00005955744478709506,B18:0.00005610375667308811):0.00000265852973118371):0.00007513430775781621):0.00080241182640416531):0.00037669332142794664):0.00050840286649229309,B51:0.00000976768913018151):0.0;'
# '((B48:0.00003909902957184906,(B64:0.00006892654277591841,(B65:0.00001169991304830036,(B11:0.00003782894234314989,(B23:0.00002257075520606242,(B7:0.00000843355467027542,B8:0.00002969574066608416):0.00000191458048020756):0.00000114768068593211):0.00001901301109015054):0.00004658929818308435):0.00001260678138773626):0.00001298046228115704,((((B13:0.00185708862613292829,B16:0.00127468956845497609):0.00077654656679172697,(B50:0.00015035037321162059,(B35:0.00013932352859303751,B59:0.00016365921800495971):0.00001707308656041147):0.00046994545309815409):0.00001689569008576868,(((((B58:0.00000100000050002909,B1:0.00024618848935342751):0.00003090697607393632,B3:0.00002583462115646557):0.00003143247214078938,((B34:0.00002644468640633620,B32:0.00005122531236540793):0.00009233804188852781,B30:0.00003756920264291922):0.00001441005666344624):0.00002033210637641083,B31:0.00005083365746347321):0.00029395904877008568,(B40:0.00017978669405508115,((B47:0.00010016344353546103,B38:0.00009438663601346181):0.00000791516086059906,B45:0.00018614652710386760):0.00005315243753183885):0.00016805788868540044):0.00015007465648871559):0.00032534959985941312,B49:0.00008499965335107176):0.00002138269003621933,B51:0.00005906076276635031):0.0;'
# '(B120:0.00001727671197251226,((((((((((B156:0.00006118558071319341,B47:0.00005000747686577315):0.00005200193331054379,((B4:0.00005420521787144568,(B146:0.00022036257690699626,B7:0.00000603363129057989):0.00008539027711667438):0.00000279074971457407,(B86:0.00006149063280793353,B101:0.00005811599342720844):0.00007013415227529132):0.00005642908874909969):0.00000273877141752970,((B130:0.00004212546227056886,B26:0.00004698846328649035):0.00000652394891136137,((B134:0.00003741824187020084,B152:0.00002937581307937704):0.00000314101435964639,B125:0.00004150530041064203):0.00001480239940026830):0.00012615117127992825):0.00008718547844915346,(((B89:0.00012124581927030777,B133:0.00013219131246434946):0.00073909433376289778,(B22:0.00055652421825995110,(B23:0.00026286433869053387,B21:0.00083001206067694056):0.00125530325876930127):0.00055710067312337522):0.00016363307710725106,((B6:0.00016599584125089184,(B84:0.00007086142991773894,(B81:0.00001967800509982177,B37:0.00001555668942755749):0.00007855415978994956):0.00010176636102377980):0.00012354438876951111,((B83:0.00005544379677381936,B82:0.00006058862415253206):0.00005272770322681279,B44:0.00008149586162578571):0.00022162947103005786):0.00043303640909177709):0.00040828638286804342):0.00008044835634892197,B74:0.00011594747680343101):0.00011032638362269623,((B57:0.00003335758287239469,B24:0.00004174535346459301):0.00001378362830570778,(B30:0.00001218571785994021,B161:0.00002162408842181342):0.00003299213684217340):0.00002864472595351346):0.00003057472902893423,(B113:0.00007146796550951428,(B25:0.00003067936280305035,(B36:0.00000798194898092488,B3:0.00001396767282842557):0.00001496313953540248):0.00000100000050002909):0.00000100000050002909):0.00000100000050002909,(B143:0.00002580644782510808,B114:0.00009265548363809927):0.00000754504937241930):0.00001250214875364111,B128:0.00001490001639763295):0.00000144789600368727,B124:0.00005244898997248456):0.00000100000050002909,B118:0.00001098566859886696):0.0;'
# '(B57:0.00000536707125664551,((((B83:0.00002318173951070520,B45:0.00003365367407355436):0.00001970126751875890,(((B84:0.00000948118264665464,B99:0.00001914943883007716):0.00001143909786772645,(((B52:0.00005185110847249898,B89:0.00001633969855456498):0.00004537804571909241,B21:0.00011406645637562866):0.00001953691402093377,(B72:0.00008420804257297685,((B51:0.00008978912752602243,B7:0.00002287080714303825):0.00029324771986313119,((B34:0.00002730728276759067,(B9:0.00003316646850023659,B29:0.00002086920547771337):0.00008771969195870552):0.00002492360304559525,(B6:0.00000881808638051164,B15:0.00006370050086096601):0.00001778245864060021):0.00012946764031754985):0.00011766209279446951):0.00001247173836728674):0.00000100000050002909):0.00000844349545619994,B36:0.00002509526808151878):0.00000667150366752081):0.00009163824755429087,B4:0.00002332293862512587):0.00000616528362070958,(B66:0.00003411246880153140,B64:0.00003762063742532716):0.00000653195343569950):0.00001043697034630331,B58:0.00004497951428308054):0.0;'
# complete path for the .fa file that contains the strains' genomes
path = 'C:/Users/Owner/Documents/UNCG/Project/BIGG_DATA/Useful_Data/Concatenates,Trees,Homoplasies/AAAYYY Run Pyvolve/'
species = ['Acinetobacter_pittii/', 'Actinobacillus_pleuropneumoniae/', 'Aeromonas_hydrophila/', 'Aggregatibacter_actinomycetemcomitans/', 'Bacillus_amyloliquefaciens/', 'Bacillus_anthracis/', 'Bacillus_ceres/', 'Bacillus_coagulans/', 'Bacillus_licheniformis/','Bacillus_methylotrophicus/', 'Bacillus_pumilus/', 'Bacillus_subtilis/', 'Bacillus_thuringiensis/']
concat = 'concat_universal.fa'
# filename = 'C:/Users/Owner/Documents/UNCG REU/Project/BIGG DATA/Useful Data/Concatenates, Trees, Homoplasies/Aayyy Clonal/Brucella_melitensis/concat_universal.fa'

for s in range(len(species)):
	# pis = iterations*[None] # list that will be populated with the simulated pi values; index = iteration - 1
	# thetas = iterations*[None] # list that will be populated with the simulated theta values; index = iteration - 1

	tree_file = open((path + species[s] + 'Universal Tree/RAxML_bestTree.tree'), 'r')
	tree_string = list(tree_file)[0]
	print(tree_string)
	print('got tree string')
	kappa_file = open((path + species[s] + 'kappa.txt'), 'r')
	kappa = float(list(kappa_file)[0])
	print('got kappa')
	strains = read_in_strains(path+species[s]+concat)
	L = genome_length(strains)
	print('read in strains')
	real_pi = pi_value(strains)
	real_theta = theta_value(strains)
	min_m = get_min_m(strains, L)
	print('min_m = ' + str(min_m))
	max_m = get_max_m(strains, L, tree_string)
	print('max_m = ' + str(max_m))
	print('got min and max m')

	tree_string = scale_branch_lengths(L, tree_string, min_m, max_m, real_pi, real_theta, kappa)
	print('got the appropriately scaled tree')

	pis = iterations*[None]
	thetas = iterations*[None]
	for i in range(iterations):
		print(i)
		# tree_file = open((path + species[s] + 'Universal Tree/RAxML_bestTree.tree'), 'r')
		# tree_string = list(tree_file)[0]
		# print(tree_string)
		# kappa_file = open((path + species[s] + 'kappa.txt'), 'r')
		# kappa = float(list(kappa_file)[0])
		# strains = read_in_strains(filename)
		# pi = pi_value(strains)
		# theta = theta_value(strains)
		# values = get_random_tree(filename, tree_string, L, kappa, i, 0) # 
		# for increment in range(increments):
		pi,theta = get_random_tree(species[s], tree_string, kappa, i) # 
		pis[i] = pi
		thetas[i] - theta
		# pi = values['pi']
		# theta = values['theta']
		# pis[i] = values['pis']
		# thetas[i] = values['thetas']
		# pis = values['pis']
		# thetas = values['thetas']
		# print('increment: ' + str(i))
		# print('pi_sim: ' + str(pi))
		# print('theta_sim: ' + str(theta))
		print('iteration complete')
	with open(('pyvolve_data_with_' + species[s][:-1] + '_universal.csv'), 'w', newline = '') as f:
		writer = csv.writer(f)
		writer.writerow(['increment', 'pi', 'theta'])
		data = [list(range(1,iterations+1)) , pis, thetas]
		data = zip(*data)
		writer.writerows(data)
		writer.writerow(['concat_universal', 'pi_real', 'theta_real'])
		# writer.writerow(['concat_universal', pi[s], theta[s]])

	print('pis: ' + str(pis))
	print('thetas: ' + str(thetas))
	print(species)
	# print('pi_real: ' + str(0.000810016)) #  str(0.000738584)) # str(0.000547989)) # str(0.000269877))
	# print('theta_real: ' + str(0.004025765)) # str(0.005666333)) # str(0.004698989)) # str(0.001704485))

