
import pyvolve
from simulations import generate_ancestor


def get_c(L, kappa):

	ancestor = generate_ancestor(L)
	print(ancestor)

	phylogeny = pyvolve.read_tree(tree = '((t1:0.5,t2:0.5)i1:0.5,(t3:0.5,t4:0.5)i2:0.5)root;')
		# '(t4:0.785,(t3:0.380,(t2:0.806,(t5:0.612,t1:0.660)i1:0.762)i2:0.921)i3:0.207)root;')
	# ((s1,s2)n1,(s3,s4)n2)n3
	pyvolve.print_tree(phylogeny)

	freqs = [0.25,0.25,0.25,0.25]

	nuc_model = pyvolve.Model('nucleotide', {'kappa':1.86836732388, 'state_freqs':freqs})

	my_partition = pyvolve.Partition(models = nuc_model, root_sequence = ancestor)

	my_evolver = pyvolve.Evolver(partitions = my_partition, tree = phylogeny)
	my_evolver() 
	# my_evolver(write_anc = True)
	strains = my_evolver.get_sequences()
	# strains = my_evolver.get_sequences(anc = True)
	strain_names = list(strains.keys())

	n = len(strain_names)
	site_counts = L*[None] # list of dictionaries to keep track of which nucleotides are at each convergent site; index = site; key = nucleotide, value = number of strains with that nucleotide
	strains_with_site = L*[None] # list of the strains that have a convergent mutation at each site; index = site
	for x in range(L):
		site_counts[x] = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
		strains_with_site[x] = []
	# c_list_matrix = [[{} for x in range(n)] for y in range(n)] # matrix of the convergent mutation sites; the (i,j) entry is a dictionary of the convergent mutation sites between strain i and strain j; key = site, value = nucleotide

	for s1 in range(n):
		strain1 = strains[strain_names[s1]]
		for s2 in range(s1,n):
			strain2 = strains[strain_names[s2]]
			for site in range(L):
				if strain1[site] == strain2[site] and strain1[site] != ancestor[site]:
					if strain1 not in strains_with_site[site]: # avoids double counting strain1 as convergent at that site
						strains_with_site[site].append(strain1)
						site_counts[site][strain1[site]] += 1
					if strain2 not in strains_with_site[site]: # avoids double counting strain2 as convergent at that site
						strains_with_site[site].append(strain2)
						site_counts[site][strain2[site]] += 1

	c_q = (n-1)*[None] # list of the number of convergent mutations between q strains; index = q - 2
	nucleotides = ['A','T','G','C']
	for x in range(n-1):
		c_q[x] = 0
	for site in site_counts:
		for base in nucleotides:
			for q in range(2,n+1):
				if site[base] == (q):
					c_q[q-2] += 1

	c = sum(c_q)
	print(c)
	return c




					# c_list_matrix[s1][s2][site] = strain1[site]
					# c_list_matrix[s2][s1][site] = strain1[site]




	# site_counts = L*[None] # list of dictionaries to keep track of which nucleotides are at each convergent site; index = site; key = nucleotide, value = number of strains with that nucleotide
	# strains_with_site = L*[None] # list of the strains that have a convergent mutation at each site; index = site
	# for x in range(L):
	# 	site_counts[x] = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
	# 	strains_with_site[x] = []
	# for strain1 in range(n): # allows each strain to be strain1
	# 	for strain2 in range(strain1+1, n): # allows each strain after strain1 to be strain2; avoids doing any combination of strains twice
	# 		for site in c_list_matrix[strain1][strain2].keys(): # goes through every sites that is convergent between strains 1 and 2
	# 			nucleotide = c_list_matrix[strain1][strain2][site] # extracts the nucleotide that strain1 and strain2 have at the site
	# 			if strain1 not in strains_with_site[site]: # avoids double counting strain1 as convergent at that site
	# 				strains_with_site[site].append(strain1)
	# 				site_counts[site][nucleotide] += 1
	# 			if strain2 not in strains_with_site[site]: # avoids double counting strain2 as convergent at that site
	# 				strains_with_site[site].append(strain2)
	# 				site_counts[site][nucleotide] += 1

	# c_q = (n-1)*[None] # list of the number of convergent mutations between q strains; index = q - 2
	# for x in range(n-1):
	# 	c_q[x] = 0
	# for site in site_counts:
	# 	for base in nucleotides:
	# 		for q in range(2,n+1):
	# 			if site[base] == (q):
	# 				c_q[q-2] += 1
	# print(c_q)




def get_accurate_c(L, kappa):

	ancestor = generate_ancestor(L)
	print(ancestor)

	# phylogeny = pyvolve.read_tree(tree = '(  (t1:0.5,t2:0.5)i1:0.5, (t3:0.5,t4:0.5)i2:0.5 ,  (t5:0.5,t6:0.5)i3:0.5, (t7:0.5,t8:0.5)i4:0.5  ) root;')
	phylogeny = pyvolve.read_tree(tree = '(  ((t7:0.5,t8:0.5)i4:0.5,(t5:0.5,t6:0.5)i3:0.5)i1:0.5, (t3:0.5,t4:0.5)i2:0.5  ) root;')

	pyvolve.print_tree(phylogeny)

	freqs = [0.25,0.25,0.25,0.25]

	nuc_model = pyvolve.Model('nucleotide', {'kappa':1.86836732388, 'state_freqs':freqs})

	my_partition = pyvolve.Partition(models = nuc_model, root_sequence = ancestor)

	my_evolver = pyvolve.Evolver(partitions = my_partition, tree = phylogeny)
	# my_evolver() 
	my_evolver(write_anc = True)
	# strains = my_evolver.get_sequences()
	strains = my_evolver.get_sequences(anc = True)
	strain_names = list(strains.keys()) # pre-order traversal of the tree
	n = len(strain_names)

	print(strain_names)

	c_sites = {}
	for key in strain_names:
		c_sites[key] = []

	site_counts = L*[None] # list of dictionaries to keep track of which nucleotides are at each convergent site; index = site; key = nucleotide, value = number of strains with that nucleotide
	strains_with_site = L*[None] # list of the strains that have a convergent mutation at each site; index = site
	for x in range(L):
		site_counts[x] = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
		strains_with_site[x] = []
	# c_list_matrix = [[{} for x in range(n)] for y in range(n)] # matrix of the convergent mutation sites; the (i,j) entry is a dictionary of the convergent mutation sites between strain i and strain j; key = site, value = nucleotide

	c = 0
	strain_names
	# a = 0
	# MRCA = [ancestor]
	# x = 1
	# y = 2
	# while 'i' + str(x) and 'i' + str(y) in strain_names:
	# 	ancestor = MRCA[a]
	# 	strain1 = strains['i' + str(x)]
	# 	strain2 = strains['i' + str(y)]
	# 	MRCA.append('i' + str(x))
	# 	MRCA.append('i' + str(y))
	# 	print(strain1)
	# 	print(strain2)
	# 	for site in range(L):
	# 		if strain1[site] == strain2[site] and strain1[site] != MRCA[site]:
	# 			c += 1
	# 			c_sites['i' + str(x)].append(site)
	# 			c_sites['i' + str(y)].append(site)
	# 	a += 1
	# 	i += 2
	# 	j += 2

get_accurate_c(10,1)


	# n = len(strain_names)
	# site_counts = L*[None] # list of dictionaries to keep track of which nucleotides are at each convergent site; index = site; key = nucleotide, value = number of strains with that nucleotide
	# strains_with_site = L*[None] # list of the strains that have a convergent mutation at each site; index = site
	# for x in range(L):
	# 	site_counts[x] = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
	# 	strains_with_site[x] = []
	# # c_list_matrix = [[{} for x in range(n)] for y in range(n)] # matrix of the convergent mutation sites; the (i,j) entry is a dictionary of the convergent mutation sites between strain i and strain j; key = site, value = nucleotide

	# for s1 in range(n):
	# 	strain1 = strains[strain_names[s1]]
	# 	for s2 in range(s1,n):
	# 		strain2 = strains[strain_names[s2]]
	# 		for site in range(L):
	# 			if strain1[site] == strain2[site] and strain1[site] != ancestor[site]:
	# 				if strain1 not in strains_with_site[site]: # avoids double counting strain1 as convergent at that site
	# 					strains_with_site[site].append(strain1)
	# 					site_counts[site][strain1[site]] += 1
	# 				if strain2 not in strains_with_site[site]: # avoids double counting strain2 as convergent at that site
	# 					strains_with_site[site].append(strain2)
	# 					site_counts[site][strain2[site]] += 1

	# c_q = (n-1)*[None] # list of the number of convergent mutations between q strains; index = q - 2
	# nucleotides = ['A','T','G','C']
	# for x in range(n-1):
	# 	c_q[x] = 0
	# for site in site_counts:
	# 	for base in nucleotides:
	# 		for q in range(2,n+1):
	# 			if site[base] == (q):
	# 				c_q[q-2] += 1

	# c = sum(c_q)
	# print(c)
	# return c




	# 				# c_list_matrix[s1][s2][site] = strain1[site]
	# 				# c_list_matrix[s2][s1][site] = strain1[site]




	# # site_counts = L*[None] # list of dictionaries to keep track of which nucleotides are at each convergent site; index = site; key = nucleotide, value = number of strains with that nucleotide
	# # strains_with_site = L*[None] # list of the strains that have a convergent mutation at each site; index = site
	# # for x in range(L):
	# # 	site_counts[x] = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
	# # 	strains_with_site[x] = []
	# # for strain1 in range(n): # allows each strain to be strain1
	# # 	for strain2 in range(strain1+1, n): # allows each strain after strain1 to be strain2; avoids doing any combination of strains twice
	# # 		for site in c_list_matrix[strain1][strain2].keys(): # goes through every sites that is convergent between strains 1 and 2
	# # 			nucleotide = c_list_matrix[strain1][strain2][site] # extracts the nucleotide that strain1 and strain2 have at the site
	# # 			if strain1 not in strains_with_site[site]: # avoids double counting strain1 as convergent at that site
	# # 				strains_with_site[site].append(strain1)
	# # 				site_counts[site][nucleotide] += 1
	# # 			if strain2 not in strains_with_site[site]: # avoids double counting strain2 as convergent at that site
	# # 				strains_with_site[site].append(strain2)
	# # 				site_counts[site][nucleotide] += 1

	# # c_q = (n-1)*[None] # list of the number of convergent mutations between q strains; index = q - 2
	# # for x in range(n-1):
	# # 	c_q[x] = 0
	# # for site in site_counts:
	# # 	for base in nucleotides:
	# # 		for q in range(2,n+1):
	# # 			if site[base] == (q):
	# # 				c_q[q-2] += 1
	# # print(c_q)


