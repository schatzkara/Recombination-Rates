
import csv
import glob
import os

def concatenate(output_file, data_path):
	with open((output_file + '.csv'), 'w', newline = '') as f:
		writer = csv.writer(f)
		j = 1
		for filename in glob.glob(os.path.join(data_path, '*.csv')):
			with open(filename) as d:
				reader = csv.reader(d)
				if(j > 1):
					next(reader)
				for row in reader:
					writer.writerow(row)
			print('file number ' + str(j))
			j += 1