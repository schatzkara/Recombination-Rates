
from convert_ID_to_c import calculate_c
from convert_ID_to_c import calculate_h_c

real_L = 1303140
L = 1000 # 1303140
mu = 1/L
generations = 300
GC_prop = 0.50
kappa = 1.86836732388 # bacillus anthracis
phi = 0.50
iterations = 1000
path = 'C:/Users/Owner/Documents/UNCG REU/Project/Recombination-Rates/Model_and_Simulations/NEW/Genomic ID Matrices/ID_Matrix_concat_bacillus_anthracis.csv'
# filename = 
# species_file = 
print('Bacillus_anthracis')
print(calculate_h_c(real_L, L, mu, generations, GC_prop, kappa, phi, iterations, path = path))

# calculate_c(species = species_file)


# real_L = 1142022
# L = 1000
# mu = 1/L
# generations = 300
# GC_prop = 0.50
# kappa = 1.53913593397
# phi = 0.50
# iterations = 1000
# path = 'C:/Users/Owner/Documents/UNCG REU/Project/Recombination-Rates/Model_and_Simulations/NEW/Genomic ID Matrices/ID_Matrix_concat_brucella_abortus.csv'
# # filename = 
# # species_file = 
# print('Brucella_abortus')
# print(calculate_c(real_L, L, mu, generations, GC_prop, kappa, phi, iterations, path = path))

# # calculate_c(species = species_file)


# # L = 1000
# # mu = 1/L
# # generations = 300
# # GC_prop = 0.50
# # kappa = 1.57113790459
# # phi = 0.50
# # iterations = 1000
# # path = 'C:/Users/Owner/Documents/UNCG REU/Project/Recombination-Rates/Model_and_Simulations/NEW/Genomic ID Matrices/ID_Matrix_concat_brucells_melitensis.csv'
# # # filename = 
# # # species_file = 
# # print('Brucella_melitensis')
# # print(calculate_c(L, mu, generations, GC_prop, kappa, phi, iterations, path = path))

# # # calculate_c(species = species_file)


# real_L = 1779675
# L = 1000
# mu = 1/L
# generations = 300
# GC_prop = 0.50
# kappa = 2.20303254463
# phi = 0.50
# iterations = 1000
# path = 'C:/Users/Owner/Documents/UNCG REU/Project/Recombination-Rates/Model_and_Simulations/NEW/Genomic ID Matrices/ID_Matrix_concat_brucella_suis.csv'
# # filename = 
# # species_file = 
# print('Brucella_suis')
# print(calculate_c(real_L, L, mu, generations, GC_prop, kappa, phi, iterations, path = path))

# # calculate_c(species = species_file)

