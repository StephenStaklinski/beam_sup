
import ete3


def get_num_migrations_from_nwk_and_labeling(nwk, labeling_tsv):
    # Read labeling file into a dict: node_name -> tissue
    node_to_tissue = {}
    with open(labeling_tsv, 'r') as f:
        for line in f:
            node, tissue = line.strip().split(' ')
            node_to_tissue[str(node)] = str(tissue)

    # Read the tree
    tree = ete3.Tree(nwk, format=1)

    # Count migrations by tree traversal
    migrations = 0
    for node in tree.traverse("preorder"):
        if node.is_root():
            continue
        parent = node.up
        parent_tissue = node_to_tissue.get(str(parent.name))
        node_tissue = node_to_tissue.get(str(node.name))
        if parent_tissue != node_tissue:
            migrations += 1

    return migrations

