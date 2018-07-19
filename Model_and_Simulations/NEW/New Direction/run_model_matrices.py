
from get_matrices import get_SCAR_matrices
from get_matrices import print_matrices


species_alignment = '/mnt/c/Users/Owner/Documents/UNCG/Project/BIGG_DATA/Useful_Data/Concatenates,Trees,Homoplasies/A_Strains_with_mu/Vibrio_cholerae/concat_universal.fa'
# raxml_path = '/mnt/c/Users/Owner/Documents/UNCG/Project/standard-RAxML'
# tree_file = 'RAxML_bestTree.tree'
# rooted_tree_file = 'RAxML_rootedTree.root'
ancestral_alignment = 'RAxML_marginalAncestralStates.anc'
# ancestral_tree_file = 'RAxML_nodeLabelledRootedTree.anc'
kappa_file = '/mnt/c/Users/Owner/Documents/UNCG/Project/BIGG_DATA/Useful_Data/Concatenates,Trees,Homoplasies/A_Strains_with_mu/Vibrio_cholerae/kappa.txt'
output_file = 'SCAR_Vibrio_cholerae_universal'
mu = 0.0000000001150
# 0.0000000001330
 # 0.0000000003350


SCAR = get_SCAR_matrices(species_alignment, ancestral_alignment, kappa_file, mu)
strain_names,S,C,A,R,RATES,average_rate = SCAR['strain_names'],SCAR['Shared'], SCAR['Convergent'], SCAR['Ancestral'], SCAR['Recombinant'], SCAR['Rates'], SCAR['average']

print_matrices(output_file, S,C,A,R,RATES, strain_names, average_rate)