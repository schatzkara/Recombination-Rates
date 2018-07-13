
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

def get_newick_format(tree_list):

	string = ';'

	root = tree_list[0]

	string = root + string

	# ((s1,s2)n1,(s3,s4)n2)n3;




# (,,(,));                               no nodes are named
# (A,B,(C,D));                           leaf nodes are named
# (A,B,(C,D)E)F;                         all nodes are named
# (:0.1,:0.2,(:0.3,:0.4):0.5);           all but root node have a distance to parent
# (:0.1,:0.2,(:0.3,:0.4):0.5):0.0;       all have a distance to parent
# (A:0.1,B:0.2,(C:0.3,D:0.4):0.5);       distances and leaf names (popular)
# (A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;     distances and all names
# ((B:0.2,(C:0.3,D:0.4)E:0.5)A:0.1)F;    a tree rooted on a leaf node (rare)




# Tree: The full input Newick Format for a single tree
#    Subtree: an internal node (and its descendants) or a leaf node
#    Leaf: a node with no descendants
#    Internal: a node and its one or more descendants
#    BranchSet: a set of one or more Branches
#    Branch: a tree edge and its descendant subtree.
#    Name: the name of a node
#    Length: the length of a tree edge.

   
# The grammar rules

#    Tree → Subtree ";" | Branch ";"
#    Subtree → Leaf | Internal
#    Leaf → Name
#    Internal → "(" BranchSet ")" Name
#    BranchSet → Branch | Branch "," BranchSet
#    Branch → Subtree Length
#    Name → empty | string
#    Length → empty | ":" number