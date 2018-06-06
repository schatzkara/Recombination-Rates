from n_change_k import n_choose_k
from n_change_k import n_change_k

x=0
n = 0
k = 0
while 0 <= n <= 10:
	while 0 <= k <= 10 and n >= k:
		while 0 <= x <= 10:
			print(str(n) + ' change ' + str(k) + ' in ' + str(x) + ' ways')
			print(n_change_k(n, k, x))
			x += 1
		x = 0
		k += 1
	k = 0
	n += 1
