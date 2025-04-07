"""Functions for statistical analysis of migration patterns."""

import numpy as np
from typing import Tuple, List, Optional
import pandas as pd
from concurrent.futures import ThreadPoolExecutor
from collections import defaultdict


def compute_posterior_mutual_info(
    trees,
    origin_tissue: str,
    threads: int = 1,
    output_file_matrix: Optional[str] = None,
    output_file_information: Optional[str] = None
) -> Tuple[float, np.ndarray, List[str]]:
    """Calculate mutual information from migration patterns in a set of trees.
    
    Args:
        trees: Dendropy TreeList object
        origin_tissue: Name of the origin tissue
        threads: Number of threads to use for parallel processing
        output_file_matrix: Optional path to save the count matrix
        output_file_information: Optional path to save the mutual information score
        
    Returns:
        Tuple containing:
            - Mutual information score
            - Count matrix as numpy array
            - List of tissue names in order
    """
    # Process trees in parallel to get migration counts
    migration_counts = defaultdict(lambda: defaultdict(int))
    tissue_types = set()

    def process_tree(tree):
        counts = defaultdict(lambda: defaultdict(int))
        types = set()
        
        for node in tree.preorder_node_iter():
            if node.parent_node:
                source = node.parent_node.annotations.get_value("location")
                target = node.annotations.get_value("location")
                types.add(source)
                types.add(target)
                counts[source][target] += 1
            else:
                target = node.annotations.get_value("location")
                types.add(origin_tissue)
                types.add(target)
                counts[origin_tissue][target] += 1
                
        return types, counts

    with ThreadPoolExecutor(max_workers=threads) as executor:
        results = executor.map(process_tree, trees)

    for t_types, m_counts in results:
        tissue_types.update(t_types)
        for source, targets in m_counts.items():
            for target, count in targets.items():
                migration_counts[source][target] += count

    # Create count matrix
    tissue_list = sorted(tissue_types - {origin_tissue})
    tissue_list.insert(0, origin_tissue)
    n = len(tissue_list)
    count_matrix = np.zeros((n, n))
    for i, source in enumerate(tissue_list):
        for j, target in enumerate(tissue_list):
            count_matrix[i, j] = migration_counts[source][target]

    # Calculate mutual information
    P = count_matrix / np.sum(count_matrix)
    p_x = np.sum(P, axis=1)
    p_y = np.sum(P, axis=0)

    MI = 0
    for i in range(P.shape[0]):
        for j in range(P.shape[1]):
            if P[i, j] > 0:
                MI += P[i, j] * np.log2(P[i, j] / (p_x[i] * p_y[j]))

    H_x = -np.sum(p_x * np.log2(p_x, where=p_x>0))
    H_y = -np.sum(p_y * np.log2(p_y, where=p_y>0))
    mutual_info = (2 * MI) / (H_x + H_y)

    # Save results if output files are provided
    if output_file_matrix:
        df = pd.DataFrame(count_matrix, index=tissue_list, columns=tissue_list)
        df.to_csv(output_file_matrix)
        
    if output_file_information:
        with open(output_file_information, 'w') as f:
            f.write(str(mutual_info))

    return mutual_info, count_matrix, tissue_list 