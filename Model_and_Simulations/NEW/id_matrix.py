# Calculates Identity matrix

import numpy as np
# f = open("file.fa")


def id_matrix(seqList): # Given ordered list of String sequences seqList
	ID = np.empty([n,n], dtype = np.float, order='C')
	strain_names = seqList.keys()
	strains = len(strain_names)
	for strain1 in range(strains):
		ID[strain1,strain1] = 1
		for strain2 in range(strain1+1,strains):
			identity = calcID(seqList[strain1], seqList[strain2])
			ID[strain1,strain2] = identity
			ID[strain2,strain1] = identity
	return ID

	# ID = []
	# strainNames = seqList.keys()
	# for n in range(len(strainNames)):
	# 	ID[n] = []
	# 	for m in range(n, len(strainNames)):
	# 		if(m==n):
	# 			ID[n][m]=1
	# 		else:
	# 			ID[n][m] = calcID(seqList[n], seqList[m])
	# 			ID[m][n] = ID[n][m] 
	# return ID

	# ID = {}
	# strainNames = seqList.keys()
	# for n in range(len(strainNames)):
		# ID[strainNames[n]] = {}
		# for m in range(n, len(strainNames)):
			# if(m==n):
				# ID[strainNames[n]][strainNames[m]]=1
			# else:
				# ID[strainNames[n]][strainNames[m]] = calcID(seqList[strainNames[n]], seqList[strainNames[m]])
				# ID[strainNames[m]][strainNames[n]] = ID[strainNames[n]][strainNames[m]] 
				
			# ID = {}
	# for strain1 in seqList.keys():
		# ID[strain1] = {}
		# for strain2 in seqList.keys():
			# ID[strain1][strain2] = calcID(seqList[strain1], seqList[strain2])
			# ID[strain2][strain1] = ID[strain1][strain2]
			# if(strain1 == strain2):
				# continue

	# return ID
def calcID (s1, s2):
	numDif = 0
	numSame = 0
	if(len(s1) != len(s2)):
		raise(ValueError('Strand Lengths not equal'))
	else:
		for i in range(s1):
			numDif += (s1[i]!=s2[i]) # Uses boolean as int, +=1 if true, +=0 if false
			numSame += (s1[i]==s2[i])
	if (numDif + numSame) != len(s1):
		raise(ValueError('numDif != numSame'))
	return numDif/len(s1)
