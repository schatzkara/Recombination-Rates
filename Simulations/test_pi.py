from pi import pi_formula
from pi import pi2

c = 0
o = 0
while 0 <= o <= 10:
	while 0 <= c <= 10:
		print('pi(' + str(c) + ',' + str(o) + ')')
		print(pi2(c,o))
		c += 1
	c = 0
	o += 1