import sys
import os
import random
from collections import Counter
from ete3 import Tree
import cassiopeia as cas
import pandas as pd
import networkx as nx


def run_parsimony_migration_cassiopeia_fitchcount(
    newick_file: str, tissues_file: str, outdir: str, run_fitchcount: bool = True
) -> None:
    """
    Runs the Cassiopeia Fitch-Hartigan parsimony migration analysis and optionally run FitchCount.

    Args:
        newick_file: Path to the input Newick tree file.
        tissues_file: Path to the CSV file mapping cells to tissues.
        outdir: Output directory to save results.
    """

    # Read in newick to ete3 tree
    ete_tree = Tree(newick_file, format=3)
    ete_tree.name = "root"

    # Load the tissues
    tissues_df = pd.read_csv(
        tissues_file, header=None, index_col=0, names=["cell", "tissue"], dtype=str
    )
    tissues_df.index = tissues_df.index.astype(str)

    # Load the tree to cassiopeia object
    tree = cas.data.CassiopeiaTree(tree=ete_tree, cell_meta=tissues_df)

    # Run fitch-hartigan to get a randomly selected parsimonious tissue labeling on the tree
    fh_tree = cas.tl.fitch_hartigan(cassiopeia_tree=tree, meta_item="tissue", copy=True)

    # Show and save results (keep in mind that only the origin is known above the root, which the parsimony here does not consider but it should not matter)
    name_map = {}
    for node in fh_tree.depth_first_traverse_nodes(postorder=False):
        label = fh_tree.get_attribute(node, "label")
        new_name = f"{node}_{label}"
        name_map[node] = new_name
    fh_tree.relabel_nodes(name_map)

    with open(f"{outdir}/cassiopeia_fitch_hartigan_result.nwk", "w") as f:
        f.write(cas.data.to_newick(fh_tree._CassiopeiaTree__network))

    # Optional: run fitch-count to get a transition matrix of tissue change frequencies across all parsimonious labelings
    if run_fitchcount:
        fc_matrix = cas.tl.fitch_count(cassiopeia_tree=tree, meta_item="tissue")
        fc_matrix.to_csv(f"{outdir}/cassiopeia_fitch_count_result.csv")


def label_tissues_fitch_parsimony(
    tree: Tree,
    tissues_df: pd.DataFrame,
    primary_tissue: str,
    threshold_num_solutions: int = 1000,
) -> tuple[Tree, list[Tree]]:
    """
    Infers ancestral tissue states for internal nodes of a phylogenetic tree using the Fitch parsimony algorithm.

    Given a tree structure and known tissue types for leaf nodes, this function assigns tissue labels to all nodes in
    the tree to minimize the number of tissue transitions (parsimony score).
    It supports enumeration of all possible equally parsimonious solutions if their number is below a specified threshold.

    Args:
        tree (Tree): ete3.Tree object representing the phylogenetic tree.
        tissues_df (pd.DataFrame): DataFrame mapping leaf node names to tissue types. Must have columns 'cell' and 'tissue'.
        primary_tissue (str): The tissue type to assign to the root node.
        threshold_num_solutions (int, optional): Maximum number of equally parsimonious solutions to enumerate. If the number of solutions exceeds this threshold, only a single solution is returned. Defaults to 1000.

    Returns:
        tuple[Tree, list[Tree]]:
            - The tree with one set of inferred tissue labels and parsimony scores assigned to each node.
            - A list of all possible equally parsimonious trees (with different tissue assignments), if their number is less than the threshold; otherwise, an empty list or a list with a single tree.

    """

    def postorder(node, tissues_df):
        if node.is_leaf():
            # Assign known tissue type from the tissues_df to the leaf node
            node.final_tissue = tissues_df.loc[
                tissues_df["cell"] == str(node.name), "tissue"
            ].values[0]
            node.tissue_set = {node.final_tissue}
            node.name = f"{node.name}_{node.final_tissue}"
            node.decision = "leaf"
        else:
            # Process all children
            children_tissue_sets = [postorder(child, tissues_df) for child in node.children]

            # Compute the possible tissues for internal nodes based on the children's tissue sets
            intersection = set.intersection(*children_tissue_sets)
            if intersection:
                node.tissue_set = intersection
            else:
                node.tissue_set = set.union(*children_tissue_sets)
        return node.tissue_set

    def preorder(node, total_solutions, primary_tissue, parent_tissue=None):
        # Leaf tissues are already known so skip them for tissue assignment
        if not node.is_leaf():
            if node.is_root():
                # The root tissue is known
                node.final_tissue = f"{primary_tissue}"
                node.decision = "root"
            elif parent_tissue and parent_tissue in node.tissue_set:
                # If parent tissue is in the node's set, choose it
                node.final_tissue = parent_tissue
                node.decision = "parent"
            else:
                # If not then make an arbitrary choice from those available and increment the parsimony score
                num_tissues = len(list(node.tissue_set))
                if num_tissues == 1:
                    node.final_tissue = list(node.tissue_set)[0]
                    node.decision = "one_option"
                else:
                    node.final_tissue = random.choice(list(node.tissue_set))
                    node.decision = "random"
                total_solutions = total_solutions * num_tissues
                node.parsimony_score += 1
            node.name = f"{node.name}_{node.final_tissue}"
            # Recursively process children
            for child in node.children:
                total_solutions = preorder(child, total_solutions, node.final_tissue)
        else:
            # Check if leaf nodes are different tissues than their parents
            if parent_tissue != node.final_tissue:
                node.parsimony_score += 1
        return total_solutions

    def traverse_all_solutions(root):
        tree = root.copy()
        all_solutions = [tree]

        for node in tree.traverse():
            if node.decision == "random":
                print("random")
                tissues = list(node.tissue_set)
                split = node.name.split("_")
                label = split[0]
                first_tissue = split[1]
                tissues_subset = [tis for tis in tissues if tis != first_tissue]
                print(tissues_subset)
                count = len(all_solutions)
                for tissue in tissues_subset:
                    print(tissue)
                    for i in range(count):
                        print(i)
                        tree_copy = all_solutions[i].copy()
                        node_copy = tree_copy.search_nodes(name=node.name)[0]
                        node_copy.name = f"{label}_{tissue}"
                        all_solutions.append(tree_copy)

        return all_solutions

    # Run the postorder to get candidate tissues at each node
    postorder(tree, tissues_df)

    # Initialize the parsimony scores
    for node in tree.traverse():
        node.parsimony_score = 0

    # Assign the ancestral tissues for each node and update the parsimony score
    num_solutions = preorder(tree, total_solutions=1, primary_tissue=primary_tissue)

    # Obtain the total parsimony score for the tree with random node selections
    total_parsimony_score = sum(node.parsimony_score for node in tree.traverse())

    # If the total number of solutions is less than the specified threshold, then re-run the preorder and enumerate all solutions to be returned as a list of trees
    # print(f"Num solutions: {num_solutions}")
    if num_solutions < threshold_num_solutions and num_solutions != 1:
        all_solutions = traverse_all_solutions(tree)
    elif num_solutions == 1:
        all_solutions = [tree]
    else:
        all_solutions = []

    return tree, all_solutions


def label_tissues_consensus(
    tree: Tree, tissues_df: pd.DataFrame, primary_tissue: str
) -> Tree:
    """
    Assigns consensus tissue labels to the nodes of a phylogenetic tree.
    This function traverses a given tree and appends tissue labels to each node's name.
    For leaf nodes, the tissue label is taken directly from the provided DataFrame.
    For internal nodes, the tissue label is determined by the most common tissue among its descendant leaves.
    In the event of a tie, the specified primary tissue is preferred as a tiebreaker; if the primary tissue is not involved in the tie, a random choice is made among the tied tissues.

    Args:
        tree (Tree): The phylogenetic tree whose nodes will be labeled.
        tissues_df (pd.DataFrame): DataFrame containing 'cell' and 'tissue' columns mapping cell names to tissue types.
        primary_tissue (str): The tissue type to prefer in the event of a tie.

    Returns:
        Tree: A copy of the input tree with node names appended with consensus tissue labels.
    """

    copy_tree = tree.copy()
    for node in copy_tree.traverse():
        # leaves are assumed to be known labels
        if node.is_leaf():
            tissue = tissues_df.loc[
                tissues_df["cell"] == str(node.name), "tissue"
            ].values[0]
            node.name = f"{node.name}_{tissue}"
        # internal nodes are chosen by consensus of leaf tissues
        else:
            node_name = node.name
            children = [
                tissues_df.loc[tissues_df["cell"] == str(leaf.name), "tissue"].values[0]
                for leaf in node.get_leaves()
            ]
            counts = Counter(children)
            most_common_elements = counts.most_common(2)
            # if there is a clear winner, choose that tissue
            if (
                len(most_common_elements) == 1
                or most_common_elements[0][1] != most_common_elements[1][1]
            ):
                consensus_tissue = most_common_elements[0][0]
            else:
                # tie breaker goes to primary
                if any(primary_tissue in elem for elem in most_common_elements):
                    consensus_tissue = primary_tissue
                # if tie doesn't involve primary, choose randomly
                else:
                    tied_elements = [elem[0] for elem in most_common_elements]
                    consensus_tissue = random.choice(tied_elements)
            # rename with consensus tissue label
            node.name = f"{node_name}_{consensus_tissue}"
    return copy_tree


def label_tissues_random(tree: Tree, tissues_df: pd.DataFrame) -> Tree:
    """
    Assigns random tissue labels to the nodes of a phylogenetic tree.
    This function traverses a given tree and appends tissue labels to each node's name.
    For leaf nodes, the tissue label is taken directly from the provided DataFrame.
    For internal nodes, a tissue label is randomly selected from the list of unique tissues present in the DataFrame.

    Args:
        tree (Tree): A phylogenetic tree whose nodes will be labeled.
        tissues_df (pd.DataFrame): A DataFrame with columns "cell" (node names) and "tissue" (tissue labels).

    Returns:
        Tree: A copy of the input tree with node names appended with tissue labels.
    """
    copy_tree = tree.copy()
    tissues = tissues_df["tissue"].unique().tolist()
    for node in copy_tree.traverse():
        # leaves are assumed to be known labels
        if node.is_leaf():
            print(node.name)
            tissue = tissues_df.loc[
                tissues_df["cell"] == str(node.name), "tissue"
            ].values[0]
            node.name = f"{node.name}_{tissue}"
        # internal nodes are chosen randomly
        else:
            random_tissue = random.choice(tissues)
            node.name = f"{node.name}_{random_tissue}"
    return copy_tree
