
file_path = 'C:/Users/Owner/Documents/UNCG REU/Project/BIGG DATA/Useful Data/Concatenates, Trees, Homoplasies/Aayyy Clonal/Bacillus_anthracis/concat_universal.fa'

# print('Reading in the DNA sequences.\n')
f = open(file_path, 'r')
f = list(f)
strains = {} # dictionary to hold all the strains (key: strain name, value: strain genome)

for line in f: 
	if(line[0] == '>'): # separates out the strain names
		key = line.strip('\n').strip('>')
		strains[key] = ''
	else:
		strains[key] += line.strip('\n') # concatenates all the lines of genome
values = list(strains.values())
length = len(values[0])
for strain in strains.values(): # makes sure that the genome length of each strain is uniform
	# print(len(strain))
	if len(strain) != length:
		print('ABORT: THE LENGTHS ARE NOT THE SAME')

g = open('new_universal_genome.fa', 'w')
remove = True
for site in range(length):
	for strain in values:
		if strain[site] != '-':
			remove = False

	if remove:
		for strain in strains:
			left = strain[:site]
			right = strain[site+1:]
			new_strain = left+right
			strains[strain] = new_strain

for key in strains.keys():
	g.write('>' + key + '\n' + strains[key])

