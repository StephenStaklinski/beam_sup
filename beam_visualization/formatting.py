from typing import Dict, Optional, List
from copy import deepcopy
from multiprocessing import Pool
import dendropy
from ete3 import Tree
import random

def _process_tree_for_consensus(tree, primary_tissue):
    """Process a single tree to get migration counts."""
    # Convert to ete3 tree
    tree_copy = deepcopy(tree)
    for node in tree_copy.preorder_node_iter():
        try:
            prediction = node.taxon.label + "_" + node.annotations.get_value("location")
            node.taxon.label = prediction
        except:
            prediction = f"node_{node.annotations.get_value('location')}"
        node.label = prediction
    
    ete_tree = Tree(tree_copy.as_string(schema="newick").replace("'", ""), format=3)
    
    # Get migration counts
    counts = {}
    for node in ete_tree.traverse():
        if node.is_root():
            continue
        node_tissue = node.name.split("_")[-1]
        parent_tissue = node.up.name.split("_")[-1]
        if node_tissue != parent_tissue:
            migration = f"{parent_tissue}_{node_tissue}"
            counts[migration] = counts.get(migration, 0) + 1
    
    # Add origin migration if needed
    if primary_tissue:
        root_tissue = ete_tree.get_tree_root().name.split("_")[-1]
        if root_tissue != primary_tissue:
            migration = f"{primary_tissue}_{root_tissue}"
            counts[migration] = counts.get(migration, 0) + 1
    
    return counts

def _process_tree_wrapper(args):
    """Wrapper function for parallel processing of trees."""
    tree, primary_tissue = args
    return _process_tree_for_consensus(tree, primary_tissue)

def get_consensus_graph(trees: dendropy.TreeList, primary_tissue: str, burnin_percent: float = 0.0, cores: int = 1) -> Dict[str, float]:
    """
    Calculate consensus graph from migration counts in trees.
    
    Args:
        trees (dendropy.TreeList): The trees to analyze
        primary_tissue (str): Primary tissue label for migration analysis
        burnin_percent (float): Percentage of trees to discard as burnin. Default is 0.0 (no burnin).
        cores (int): Number of CPU cores to use for parallel processing. Default is 1 (single core).
        
    Returns:
        Dict[str, float]: Dictionary mapping migration patterns to their probabilities.
        Each key is in the format "source_tissue_dest_tissue_count" and the value is the probability
        of that migration pattern occurring in the posterior distribution.
    """
    print("\nCalculating consensus graph...")
    
    # Process posterior trees
    num_trees = len(trees)
    num_discard = round(num_trees * burnin_percent)
    trees_to_analyze = trees[num_discard:]
    
    print(f"  Analyzing {len(trees_to_analyze)} trees (after {num_discard} burnin)")
    
    # Process trees in parallel
    with Pool(processes=cores) as pool:
        all_counts = pool.map(
            _process_tree_wrapper,
            [(tree, primary_tissue) for tree in trees_to_analyze]
        )
    
    # Calculate consensus graph
    prob = 1 / len(all_counts)
    consensus_graph = {}
    for counts in all_counts:
        for migration, count in counts.items():
            for i in range(1, count + 1):
                edge = f"{migration}_{i}"
                consensus_graph[edge] = consensus_graph.get(edge, 0) + prob
    
    print("  Consensus graph calculation complete")
    return consensus_graph

def sample_trees(trees: dendropy.TreeList, n: int = 1, burnin_percent: float = 0.1, output_prefix: Optional[str] = None) -> List[str]:
    """
    Sample trees from the posterior distribution.
    
    Args:
        trees: Dendropy TreeList object
        n: Number of trees to sample
        burnin_percent: Percentage of trees to discard as burnin
        output_prefix: Optional prefix for output files. If None, trees are not written to file.
        
    Returns:
        List[str]: List of newick strings for sampled trees with location annotations preserved
    """
    # Calculate burnin
    num_discard = round(len(trees) * burnin_percent)
    available_trees = range(num_discard, len(trees))
    
    # Sample tree indices
    sampled_indices = random.sample(list(available_trees), k=min(n, len(available_trees)))
    
    # Process only the sampled trees
    sampled_trees = []
    for idx in sampled_indices:
        # Create a copy to modify
        tree_copy = trees[idx].clone()
        
        # Process each node to add names and preserve location annotations
        i = 1
        for node in tree_copy.preorder_node_iter():
            location = node.annotations.get_value("location")
            if node.taxon:
                if not node.taxon.label:
                    node.taxon.label = f"node{i}"
                    i += 1
                # Only add location annotation if it's not already there
                if "[&location=" not in node.taxon.label:
                    node.taxon.label = f"{node.taxon.label}[&location=\"{location}\"]"
            else:
                if not node.label:
                    node.label = f"node{i}"
                    i += 1
                # Only add location annotation if it's not already there
                if "[&location=" not in node.label:
                    node.label = f"{node.label}[&location=\"{location}\"]"
        
        # Convert to newick string
        sampled_trees.append(tree_copy.as_string(schema="newick").strip())
    
    # Write to file if output_prefix is specified
    if output_prefix:
        output_file = f"{output_prefix}_sampled_trees.txt"
        with open(output_file, 'w') as f:
            for tree in sampled_trees:
                f.write(tree + '\n')
    
    return sampled_trees 