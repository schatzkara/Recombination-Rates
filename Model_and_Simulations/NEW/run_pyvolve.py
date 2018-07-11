
import csv
from get_c_from_pyvolve import get_c


L = 10000 # 1303140
kappa = 1.86836732388
iterations = 1000

with open(('c_from_pyvolve_' + str(L) + '_' + str(kappa) + '.csv'), 'w', newline = '') as f:
	writer = csv.writer(f)
	writer.writerow(['Iteration', 'c', 'L', 'kappa'])
	for i in range(iterations):
		writer.writerow([i+1, get_c(L, kappa), L, kappa])