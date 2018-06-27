mu = 10/1000

s = 0
for n in range(300):
	s += (mu**n) * n
	print(s)
print(mu/(1-mu)**2)