import math

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
	if(type(overlaps) is not str and type(strain1) is not str and type(strain2) is not str and type(combos1) is not str and type(combos2) is not str):
		prob = 1/(combos1 * combos2)
		return (overlaps * strain2 * strain1 * prob)
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
	if(type(matches) is not str and type(strain1) is not str and type(strain2) is not str and type(combos1) is not str and type(combos2) is not str):
		prob = 1/(combos1 * combos2)
		answer1 =  (matches * strain2 * strain1 * prob)
	else: 
		answer1 = 'undefined'
	if(type(overlapping_prob(c, m, n, L)) is not str):
		answer2 = (((1/3) ** c) * overlapping_prob(c, m, n, L))
	else:
		answer2 = 'undefined'
	# print(str(answer1) + "    " + str(answer2))
	if(type(answer1) is not str and type(answer2) is not str):
		if(abs(answer1 - answer2) > .000001):
			return '			YOU DONE GOOFED'
	# 	return answer1 - answer2
	# else:
	# 	return 'undefined'
	return answer1