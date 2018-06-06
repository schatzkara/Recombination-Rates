from n_change_k import n_choose_k

# def pi_wrong(c,o):
# 	answer = 0
# 	for k in range(c,(o+1)):
# 		answer += (((-1)**c) * (n_choose_k(o, k)) * ((-1/3)**k))
# 	return answer

def pi_formula(c,o):
	combos = n_choose_k(o,c)
	match = (1/3)**c
	diff = (2/3)**(o-c)
	if(type(combos) is not str):
		prob = 1/((n_choose_k(o,c))**2)
		return (combos * prob * match * diff)
	else:
		return 'undefined'