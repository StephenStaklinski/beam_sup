import pandas as pd
import cassiopeia as cas
from ete3 import Tree


def infer_parsimony_tree_cassiopeia_greedy(
    character_matrix_tsv: str, outdir: str
) -> None:
    """
    Infers a parsimony tree from a character matrix TSV file using Cassiopeia's greedy solver,
    and writes the resulting tree in Newick format with internal node labels.

    Args:
        character_matrix_tsv (str): Path to the input character matrix TSV file.
        outdir (str): Output directory prefix for saving the inferred tree.
    """
    # read in final matrix
    final_matrix = pd.read_csv(character_matrix_tsv, sep="\t", index_col=0)
    final_matrix.index = final_matrix.index.astype(str)

    # solve cassiopeia greedy
    reconstructed_tree = cas.data.CassiopeiaTree(
        character_matrix=final_matrix, missing_state_indicator=-1
    )
    greedy_solver = cas.solver.VanillaGreedySolver()
    greedy_solver.solve(reconstructed_tree)

    # make ete3 tree to write newick with internal node labels
    connections = reconstructed_tree.edges
    tree = Tree.from_parent_child_table(connections)

    # rename internal nodes from cassiopeia defaults
    i = 0
    for node in tree.traverse():
        if not node.is_leaf():
            node.name = f"node{i}"
            i += 1

    out_tree_infer = outdir + "/cassiopeia_greedy_inferred.nwk"
    with open(out_tree_infer, "w") as it:
        it.write(tree.write(format=8))
