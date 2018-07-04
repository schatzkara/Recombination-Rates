
from model import expected_h_c_mu_star
from model import expected_h_c_mg

n = 5
L = 1000
mu = 1/1000
generations = 300
kappa = 3
phi = 0.5


print(expected_h_c_mu_star(n, L, mu, generations, kappa, phi))

print(expected_h_c_mg(n, L, mu, generations, kappa, phi))