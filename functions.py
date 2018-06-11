import math
import numpy as np

# time complexity: O(1)
def n_choose_k(n, k):
	if n >= k:
		if k < 0: 
			return 0
		else:
			num = 1
			for x in range(n-k+1,n+1):
				num = num * x
			denom = math.factorial(k)
			return num/denom
	elif k < 0:
		return 0
	else:
		return 0

# time complexity: O(1)
def multiple(n, w,x,y):
	if n >= (w + x + y):
		if w < 0 or (x-w) < 0 or (y-w) < 0:
			return 0
		else:
			num = 1
			for m in range((n-w-x-y+1),n+1):
				num = num * m
			denom = math.factorial(w) * math.factorial(x-w) * math.factorial(y-w)
			return num/denom
	elif w < 0 or (x-w) < 0 or (y-w) < 0:
		return 0
	else: 
		return 0

# time complexity: O(1)
def n_change_k(n, k, x):
	return ((x**k) * n_choose_k(n, k))

# time complexity: O(1)
def poisson_prob(m, expected):
	return (math.exp(-expected) * ((expected ** m)/(math.factorial(m))))

# time complexity: O(1)
def pi_formula(c, o):
	combos = n_choose_k(o, c)
	match = (1/3)**c
	diff = (2/3)**(o-c)
	return (combos * match * diff)

# time complexity: O(1)
def overlapping_prob(o, m, n, L):
	if L <= 0:
		return 0
	if(m < n):
		temp = m
		m = n
		n = temp
	choose = multiple(L, o, m-o, n-o)
	overlaps = n_choose_k(L, o)
	strain2 = n_choose_k((L-o), (n-o))
	strain1 = n_choose_k((L-n), (m-o))
	combos1 = n_choose_k(L, m)
	combos2 = n_choose_k(L, n)
	if(combos1 == 0 or combos2 == 0):
		return 0
	else: 
		return (overlaps * strain2 * strain1 * (1/(combos1 * combos2)))

# time complexity: O(1)
def matching_prob(c, m, n, L, x):
	if L <= 0:
		return 0
	if(m < n):
		temp = m
		m = n
		n = temp
	matches = n_change_k(L, c, x)
	strain2 = n_change_k((L-c), (n-c), x)
	strain1 = n_change_k((L-n), (m-c), x)
	combos1 = n_change_k(L, m, x)
	combos2 = n_change_k(L, n, x)
	if(combos1 == 0 or combos2 == 0):
		return 0
	else:
		return (matches * strain2 * strain1 * (1/(combos1 * combos2)))
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
	# return answer1

# time complexity: O(n), where n is the smaller of m and n
def expected_overlaps(m, n, L):
	expected = np.longdouble(0)
	x = min(m,n)
	for o in range(0,(x+1)):
		expected += np.longdouble(o * overlapping_prob(o, m, n, L))
	return expected

# time complexity: O(n), where n is the smaller of m and n
def expected_matches(m, n, L):
	expected = np.longdouble(0)
	y = min(m,n)
	for c in range(0,(y+1)):
		expected += np.longdouble(c * matching_prob(c, m, n, L, 3))
	return expected

# time complexity: O(n^3), where n is L
def prob_cm(c, mu, L):
	expected_value = mu*L
	innersum = 0
	outersum = 0
	prob = []
	for o in range(L+1):
		for m1 in range(L+1):
			for m2 in range(L+1):
				if L > 0:
					x = poisson_prob(m1, expected_value)
					y = poisson_prob(m2, expected_value)
					z = overlapping_prob(o, m1, m2, L)
					innersum += (x * y * z)
		outersum = (innersum * pi_formula(c, o))
		prob.append(outersum)
		innersum = 0
	return sum(prob)

# time complexity: O(n^4), where n is L
def expected_cms(mu, L):
	value = 0
	total = 0
	for c in range(L+1):
		value += (c * prob_cm(c, mu, L))
		total += prob_cm(c, mu, L)
	return value
