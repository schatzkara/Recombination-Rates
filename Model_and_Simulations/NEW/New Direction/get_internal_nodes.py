file_path = 'C:/Users/Owner/Documents/UNCG/Project/standard-RAxML/RAxML_marginalAncestralStates.anc' 


def get_internal_nodes(file_path):
	f = open(file_path, 'r')
	f = list(f)

	internal_nodes = {}

	for x in range(len(f)):
		# print(f[x][:10])
		node_name, genome = f[x].strip('\n').split(' ')
		# print(node_name)
		# print(genome[:10])
		internal_nodes[node_name] = genome

	# print(internal_nodes)
	return internal_nodes

# get_internal_nodes(file_path)