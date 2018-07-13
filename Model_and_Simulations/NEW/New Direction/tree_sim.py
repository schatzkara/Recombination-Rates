

def make_tree(tree_file):

	tree_file = open(tree_file, 'r')
	tree_file = list(tree_file)

	root = None
	parents = {}
	children = {}

	for line in tree_file:
		current = line.strip('\n').split('\t')
		parent = current[0]
		child1 = current[1]
		child2 = current[2]

		parents[parent] = [child1,child2]
		children[child1] = parent
		children[child2] = parent

	for parent in parents:
		if parent not in children.keys():
			root = parent

	tree = []
	tree.append(root)
	# tree.append(parents[root][0])
	# tree.append(parents[root][1])

	for i in range(len(tree_file)):
		parent = tree[i]
		tree.append(parents[parent][0])
		tree.append(parents[parent][1])

	print_tree(tree)

def print_tree(tree_list):
	root = tree_list[0]
	print(root)
	# parents = tree_list[0:int(len(tree_list)/2)]
	current_parents = [root]
	while current_parents != []:
		for j in range(len(current_parents)):
			i = tree_list.index(current_parents[j])
			c1 = 2*i + 1
			c2 = 2*i + 2
			current_parents.remove(current_parents[j])
			if c1 < len(tree_list) and c2 < len(tree_list):
				child1 = tree_list[c1]
				child2 = tree_list[c2]
				current_parents.append(child1)
				current_parents.append(child2)
			print(current_parents)
			# print(child1 + '\t' + child2)


make_tree('test_tree.txt')
