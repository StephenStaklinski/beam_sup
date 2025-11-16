
import sys
from Bio import Phylo
from ete3 import Tree
import dendropy


def convert_newick_to_nexus(newick_file: str, outfile: str) -> None:
    """
    Converts a Newick tree file to Nexus format.

    Args:
        newick_file (str): Path to the input Newick file.

    Returns:
        str: Path to the output Nexus file.
    """
    tree = Phylo.read(newick_file, "newick")
    Phylo.write(tree, outfile, "nexus")


def convert_nexus_to_newick(nexus_file: str, outfile: str) -> None:
    """
    Converts a Nexus tree file to Newick format, removing annotations and features.

    Args:
        nexus_file (str): Path to the input Nexus file.
        outfile (str): Path to the output Newick file.
    """

    def remove_annotations_and_features(clade):
        if hasattr(clade, "comment"):
            del clade.comment
        if hasattr(clade, "branch_length"):
            del clade.branch_length

    tree = Phylo.read(nexus_file, "nexus")
    tree.rooted = True  # ensure the tree is rooted
    tree.format = "newick"
    for clade in tree.find_clades():
        remove_annotations_and_features(clade)
    Phylo.write(tree, outfile, "newick", plain=True)


def get_num_migrations_from_nwk_and_labeling(nwk: str, labeling_tsv: str) -> int:
    """
    Counts the number of migration events in a phylogenetic tree based on node tissue labeling provided
    seperately in a tsv file. A migration event is defined as a transition between different tissue
    types along the branches of the tree.

    Args:
        nwk (str): Newick string representing the phylogenetic tree.
        labeling_tsv (str): Path to a TSV file mapping node names to tissue types. Each line should contain 'node_name tissue' separated by a space.

    Returns:
        int: The number of migration events detected in the tree.
    """

    # Read labeling file into a dict: node_name -> tissue
    node_to_tissue = {}
    with open(labeling_tsv, "r") as f:
        for line in f:
            node, tissue = line.strip().split(" ")
            node_to_tissue[str(node)] = str(tissue)

    # Read the tree
    tree = Tree(nwk, format=1)

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


def get_mig_comig_counts_and_topologies_from_nwk(
    newick: str, primary_tissue: str
) -> tuple[int, int, int, bool, bool, str]:
    """
    Analyzes a phylogenetic tree in Newick format to count migration and co-migration events
    between tissues, and to determine clonality and other migration-related properties. Requires
    the newick string to have node names formatted as 'nodeID_tissue'.

    Args:
        newick (str): The Newick string representing the phylogenetic tree.
        primary_tissue (str): The name of the primary tissue to use as the root context.

    Returns:
        tuple[int, int, int, bool, bool, str]: A tuple containing:
            - migration_count (int): Total number of migration events (edges where tissue changes).
            - comigration_count (int): Number of unique migration edges (co-migration events).
            - num_multiedges (int): Number of repeated migration edges (multi-edges).
            - met_to_met (bool): True if there is a migration between two metastatic tissues.
            - reseeding (bool): True if there is a migration back to the primary tissue.
            - clonality (str): "Monoclonal" if no multi-edges, "Polyclonal" otherwise.
    """

    tree = Tree(newick, format=8)

    root = tree.get_tree_root()
    if not root.name:
        tree.get_tree_root().name = f"root_{primary_tissue}"

    migration_count = 0
    comigration_count = 0
    num_multiedges = 0
    met_to_met = False
    reseeding = False
    edges = []
    multiedges_already_checked = []

    for node in tree.traverse():
        if node.is_root():
            continue
        if node.up.is_root():
            parent_tissue = primary_tissue
        else:
            parent_tissue = node.up.name.split("_")[-1]
        child_tissue = node.name.split("_")[-1]
        if parent_tissue != child_tissue:
            migration_count += 1
            edge = f"{parent_tissue}_{child_tissue}"
            if edge not in edges:  # All unique edges are considered a co-migration
                comigration_count += 1
            if edge in edges and edge not in multiedges_already_checked:
                num_multiedges += 1
                multiedges_already_checked.append(edge)
            if (
                not met_to_met
                and parent_tissue != primary_tissue
                and child_tissue != primary_tissue
            ):
                met_to_met = True
            if not reseeding and child_tissue == primary_tissue:
                reseeding = True
            edges.append(edge)

    if num_multiedges != 0:
        clonality = "Polyclonal"
    else:
        clonality = "Monoclonal"

    return (
        migration_count,
        comigration_count,
        num_multiedges,
        met_to_met,
        reseeding,
        clonality,
    )


def annotate_tree_with_tissues(nwk_file: str, tissues_file: str, out_file: str) -> None:
    """
    Annotates each node in a Newick tree with its tissue type from a mapping file and writes the result.

    Args:
        nwk_file (str): Path to the input Newick tree file.
        tissues_file (str): Path to the TSV file mapping node names to tissue types.
        out_file (str): Path to the output Newick file with annotated node names.
    """

    tree = Tree(nwk_file, format=3)

    tissues = {}
    with open(tissues_file, "r") as file:
        for line in file:
            fields = line.strip().split(" ")
            tissues[fields[0]] = fields[1]

    for node in tree.traverse():
        name = node.name
        node.name = name + "_" + tissues[name]

    tree.write(outfile=out_file, format=8)
    

def remove_tissues_from_tree(nwk_file: str, out_file: str) -> None:
    """
    Removes appended tissue types from each node name in a Newick tree and writes the result.

    Args:
        nwk_file (str): Path to the input Newick tree file.
        out_file (str): Path to the output Newick file with cleaned node names.
    """

    # Read the tree (DendroPy preserves underscores by default)
    tree = dendropy.Tree.get(
        path=nwk_file,
        schema="newick",
        preserve_underscores=True
    )

    # Modify leaf and internal node labels
    for node in tree:
        if node.label:
            node.label = "_".join(node.label.split("_")[0:-1])
        if node.taxon and node.taxon.label:
            node.taxon.label = "_".join(node.taxon.label.split("_")[0:-1])

    # Write out to file
    tree.write(
        path=out_file,
        schema="newick",
        suppress_rooting=True,
        suppress_internal_node_labels=False,
        unquoted_underscores=True
        
    )


def format_edge_list_to_newick(
    edge_list_file: str,
    outfile: str,
) -> None:
    """
    Convert an edge list to a newick tree and write to outfile.

    Args:
        edge_list_file (str): Input edge list file.
        outfile (str): Ouput newick file
    """
    edges = []
    with open(edge_list_file, "r") as file:
        for line in file:
            fields = line.strip().split()
            edge = (fields[0], fields[1])
            edges.append(edge)

    tree = Tree.from_parent_child_table(edges)
    tree.write(outfile=outfile, format=8)
