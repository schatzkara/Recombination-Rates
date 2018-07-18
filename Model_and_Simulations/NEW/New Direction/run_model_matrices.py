
from get_matrices import get_SCAR_matrices
from get_matrices import print_matrices


species_alignment = '/mnt/c/Users/Owner/Documents/UNCG/Project/BIGG_DATA/Useful_Data/Concatenates,Trees,Homoplasies/Aayyy_Clonal/Bacillus_anthracis/concat_universal.fa'
# raxml_path = '/mnt/c/Users/Owner/Documents/UNCG/Project/standard-RAxML'
# tree_file = 'RAxML_bestTree.tree'
# rooted_tree_file = 'RAxML_rootedTree.root'
ancestral_alignment = 'RAxML_marginalAncestralStates.anc'
# ancestral_tree_file = 'RAxML_nodeLabelledRootedTree.anc'
kappa_file = '/mnt/c/Users/Owner/Documents/UNCG/Project/BIGG_DATA/Useful_Data/Concatenates,Trees,Homoplasies/Aayyy_Clonal/Bacillus_anthracis/kappa.txt'
output_file = 'SCAR_Bacillus_anthracis'


SCAR = get_SCAR_matrices(species_alignment, ancestral_alignment, kappa_file)
strain_names,S,C,A,R = SCAR['strain_names'],SCAR['Shared'], SCAR['Convergent'], SCAR['Ancestral'], SCAR['Recombinant']

print_matrices(output_file, S,C,A,R, strain_names)