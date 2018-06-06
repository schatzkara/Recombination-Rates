import math

def n_change_k(n , k , x):
	return ((x**k) * n_choose_k(n, k))

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
