import math
import numpy as np

def n_choose_k(n, k):
	if n >= k:
		if k < 0: 
			return 0
		else:
			num = math.factorial(n)
			denom = math.factorial(k) * math.factorial(n-k)
			answer = num/denom
			return answer
	else:
		return "undefined"

def n_change_k(n, k, x):
	return ((x**k) * n_choose_k(n, k))

# def pi_wrong(c,o):
# 	answer = 0
# 	for k in range(c,(o+1)):
# 		answer += (((-1)**c) * (n_choose_k(o, k)) * ((-1/3)**k))
# 	return answer

def pi_formula(c, o):
	combos = n_choose_k(o,c)
	match = (1/3)**c
	diff = (2/3)**(o-c)
	if(type(combos) is not str):
		prob = 1/((n_choose_k(o,c))**2)
		return (combos * prob * match * diff)
	else:
		return 'undefined'

def overlapping_prob(o, m, n, L):
	if(m < n):
		temp = m
		m = n
		n = temp
	overlaps = n_choose_k(L, o)
	strain2 = n_choose_k((L-o), (n-o))
	strain1 = n_choose_k((L-n), (m-o))
	combos1 = n_choose_k(L, m)
	combos2 = n_choose_k(L, n)
	# logo = math.log10(overlaps)
	# logs1 = math.log10(strain1)
	# logs2 = math.log10(strain2)
	# logc1 = math.log10(combos1)
	# logc2 = math.log10(combos2)
	# print('overlaps: ' + str(overlaps) + 'strain2: ' + str(strain2) + 'strain1: ' + str(strain1) + 'combos1: ' + str(combos1) + 'combos2: ' + str(combos2))
	if(type(overlaps) is not str and type(strain1) is not str and type(strain2) is not str and type(combos1) is not str and type(combos2) is not str):
		prob = np.longdouble(1/(combos1 * combos2))
		# print('prob: ' + str(prob))
		# print('answer: ' + str(prob * overlaps * strain2 * strain1))
		answer = np.longdouble(prob * overlaps * strain2 * strain1)
		return answer
		# return (logo + logs1 + logs2 + logc1 + logc2)
	else: 
		return 'undefined'

def matching_prob(c, m, n, L, x):
	if(m < n):
		temp = m
		m = n
		n = temp
	matches = n_change_k(L, c, x)
	strain2 = n_change_k((L-c), (n-c), x)
	strain1 = n_change_k((L-n), (m-c), x)
	combos1 = n_change_k(L, m, x)
	combos2 = n_change_k(L, n, x)
	# logm = math.log10(matches)
	# logs1 = math.log10(strain1)
	# logs2 = math.log10(strain2)
	# logc1 = math.log10(combos1)
	# logc2 = math.log10(combos2)
	if(type(matches) is not str and type(strain1) is not str and type(strain2) is not str and type(combos1) is not str and type(combos2) is not str):
		prob = np.longdouble(1/(combos1 * combos2))
		answer1 = np.longdouble(prob * matches * strain2 * strain1)
		# loganswer = (logm + logs1 + logs2 + logc1 + logc2)
	else: 
		answer1 = 'undefined'
		# loganswer = 'undefined'
	# if(type(overlapping_prob(c, m, n, L)) is not str):
	# 	answer2 = (((1/3) ** c) * overlapping_prob(c, m, n, L))
	# else:
	# 	answer2 = 'undefined'
	# # print(str(answer1) + "    " + str(answer2))
	# if(type(answer1) is not str and type(answer2) is not str):
	# 	if(abs(answer1 - answer2) > .000001):
	# 		return '			YOU DONE GOOFED'
	# 	return answer1 - answer2
	# else:
	# 	return 'undefined'
	return answer1
	# return loganswer

def expected_overlaps(m, n, L):
	expected = np.longdouble(0)
	x = min(m,n)
	for o in range(0,(x+1)):
		# print(overlapping_prob(o, m, n, L))
		expected += np.longdouble(o * overlapping_prob(o, m, n, L))
		# expected += (o * (10 ** overlapping_prob(o, m, n, L)))
		# print(expected)
		#log = math.log10(expected)
	return expected

def expected_matches(m, n, L):
	expected = np.longdouble(0)
	y = min(m,n)
	for c in range(0,(y+1)):
		# print(matching_prob(c, m, n, L, 3))
		expected += np.longdouble(c * matching_prob(c, m, n, L, 3))
	return expected
