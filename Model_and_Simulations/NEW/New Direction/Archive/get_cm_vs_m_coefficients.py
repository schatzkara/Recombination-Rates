

from model import new_prob_c
import csv

L = [1000]
mu = [1/1000]
kappa = [1.0,2.0,3.0]
phi = [1/2]
generations = 300


with open(('all_cm_vs_m_coefficients.csv'), 'w', newline = '') as f:
	writer = csv.writer(f)
	writer.writerow(['L', 'kappa', 'phi', 'generation', 'a'])
	for l in L:
		for m in mu:
			for k in kappa:
				for p in phi:
					for g in range(generations):
						print(g)
						writer.writerow([l, k, p, g+1, new_prob_c(m, k, p, g+1)])



	# 				with open(('all_averages_for_SNP_sim_' + str(l) + '_' + str(g) + '.csv'), 'w', newline = '') as f: # opens the file that the data will be written to
	# 				writer = csv.writer(f)
	# 				writer.writerow(['Mutations on Each Strain', 'Average CMs', 'STDev', 'GC%', 'Kappa']) # column headers
	# 				j = 1 # counter for the number of files
	# 				for filename in glob.glob(os.path.join(path, '*.csv')): # iterates over every .csv file in the path
	# 					gc = filename[-11:-8] # extracts the GC% from the file name
	# 					k = filename[-7:-4] # extracts the kappa value from the file name
	# 					with open(filename) as d: # reads the file
	# 						reader = csv.reader(d)
	# 						next(reader) # skips the column headers
	# 						for row in reader: # writes each row to the new file
	# 							writer.writerow([row[0], row[1], row[2], gc, k])
	# 					print('FILE NUMBER ' + str(j))
	# 					j += 1


	# prob_c(mu, kappa, phi, generations)