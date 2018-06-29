from run_sims import run_ID_matrix_sim

n = [5,10,15,20]
L = [1000]
generations = [300]
mu = [1/1000]
kappa = [1.0,2.0,3.0,4.0,5.0,6.0]
phi = [0.25,0.50,0.75,1.00]
iterations = 1000

run_ID_matrix_sim(n, L, generations, mu, kappa, phi, iterations)