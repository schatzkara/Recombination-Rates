from run_sims import run_c_q_sim

n = [5,10,15,20]
L = [1000]
mu = [1/1000]
generations = [300]
kappa = [1.00,2.00,3.00,4.00,5.00,6.00]
phi = [0.25,0.50,0.75,1.00]
iterations = 1000

run_c_q_sim(n, L, mu, generations, kappa, phi, iterations)

