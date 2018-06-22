from functions import overlapping_prob
from functions import matching_prob

o = 0
m = 0
n = 0
L = 1
while 0 <= L <= 10:
	while 0 <= m <= L:
		while 0 <= n <= L:
			while 0 <= o <= L:
				print('P(' + str(o) + '; ' + str(m) + ', ' + str(n) + ', ' + str(L) + ') <- for overlaps')
				print(overlapping_prob(o, m, n , L))
				o += 1
			o = 0
			n +=1
		n = 0
		m += 1
	m = 0
	L += 1

o = 0
m = 0
n = 0
L = 4
x = 3
# while 0 <= L <= 10:
while 0 <= m <= L:
	while 0 <= n <= L:
		while 0 <= o <= L:
			print('P(' + str(o) + '; ' + str(m) + ', ' + str(n) + ', ' + str(L) + ') <- for matches')
			print(matching_prob(o, m, n , L, x))
			o += 1
		o = 0
		n +=1
	n = 0
	m += 1
# m = 0
# L += 1


# matching_prob(0, 3, 2, 4, 3)