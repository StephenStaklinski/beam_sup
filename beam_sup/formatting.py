
from typing import Optional, Dict, List, Union, Tuple
from copy import deepcopy
from multiprocessing import Pool
import dendropy
from ete3 import Tree
import random
import pickle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict

from . import config


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


def get_consensus_graph(
    trees: dendropy.TreeList,
    primary_tissue: str,
    burnin_percent: float = config.DEFAULT_BURNIN_PERCENT,
    cores: int = config.DEFAULT_CORES,
) -> Dict[str, float]:
    """
    Calculate consensus graph from migration counts in trees.

    Args:
        trees (dendropy.TreeList): The trees to analyze
        primary_tissue (str): Primary tissue label for migration analysis
        burnin_percent (float): Percentage of trees to discard as burnin. config.DEFAULT is 0.0 (no burnin).
        cores (int): Number of CPU cores to use for parallel processing. config.DEFAULT is 1 (single core).

    Returns:
        Dict[str, float]: Dictionary mapping migration patterns to their probabilities.
        Each key is in the format "source_target_count" and the value is the probability
        of that migration pattern occurring in the posterior distribution.

    Raises:
        ValueError: If trees is None or empty
    """
    if not isinstance(trees, dendropy.TreeList):
        raise ValueError("Trees must be a dendropy.TreeList object")
        
    if len(trees) == 0:
        raise ValueError("No trees to analyze")

    print("\nCalculating consensus graph...")

    # Process posterior trees
    num_trees = len(trees)
    num_discard = round(num_trees * burnin_percent)
    trees_to_analyze = trees[num_discard:]

    print(f"  Analyzing {len(trees_to_analyze)} trees (after {num_discard} burnin)")

    # Process trees in parallel
    with Pool(processes=cores) as pool:
        all_counts = pool.map(
            _process_tree_wrapper, [(tree, primary_tissue) for tree in trees_to_analyze]
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


def sample_trees(
    trees: dendropy.TreeList,
    n: int = config.DEFAULT_NUM_SAMPLES,
    burnin_percent: float = config.DEFAULT_BURNIN_PERCENT,
    output_prefix: Optional[str] = None,
) -> List[str]:
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
    sampled_indices = random.sample(
        list(available_trees), k=min(n, len(available_trees))
    )

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
                    node.taxon.label = f'{node.taxon.label}[&location="{location}"]'
            else:
                if not node.label:
                    node.label = f"node{i}"
                    i += 1
                # Only add location annotation if it's not already there
                if "[&location=" not in node.label:
                    node.label = f'{node.label}[&location="{location}"]'

        # Convert to newick string
        sampled_trees.append(tree_copy.as_string(schema="newick").strip())

    # Write to file if output_prefix is specified
    if output_prefix:
        output_file = f"{output_prefix}_sampled_trees.txt"
        with open(output_file, "w") as f:
            for tree in sampled_trees:
                f.write(tree + "\n")

    return sampled_trees


def get_all_posterior_metastasis_times(
    trees: dendropy.TreeList,
    total_time: float,
    primary_tissue: Optional[str] = None,
    burnin_percent: float = config.DEFAULT_BURNIN_PERCENT,
    output_prefix: Optional[str] = None,
) -> Dict[str, Dict[str, Tuple[float, float]]]:
    """
    Calculate metastasis times for all trees in the posterior distribution.

    Args:
        trees: Dendropy TreeList object containing the posterior trees
        total_time: Total time of the experiment
        primary_tissue: Optional tissue label for the primary site. If not provided, must be set globally.
        burnin_percent: Percentage of trees to discard as burnin
        output_prefix: Optional prefix for output files. If provided, results will be written to {output_prefix}.pkl

    Returns:
        Dict[str, Dict[str, Tuple[float, float]]]: Dictionary mapping tree labels to dictionaries of metastasis events.
        Each metastasis event is a tuple of (start_time, end_time) for the migration.
    """

    # Apply burnin
    num_discard = round(len(trees) * burnin_percent)
    trees_to_analyze = trees[num_discard:]

    all_met_events = {}

    for tree in trees_to_analyze:
        # Create a copy to modify
        tree_copy = tree.clone()

        # Process each node to add names and preserve location annotations
        i = 0
        for node in tree_copy.preorder_node_iter():
            try:
                prediction = (
                    node.taxon.label + "_" + node.annotations.get_value("location")
                )
                node.taxon.label = prediction
            except Exception:
                prediction = f"node{i}" + "_" + node.annotations.get_value("location")
                i += 1
            node.label = prediction

        # Convert to newick string and create ete3 tree
        newick = tree_copy.as_string(
            schema="newick",
            suppress_edge_lengths=False,
            node_label_element_separator=",",
        )
        newick = newick.replace("'", "")  # Remove quoted node names
        ete_tree = Tree(newick, format=3)

        # Verify tree is ultrametric
        root = ete_tree.get_tree_root()
        leaf_distances = set()
        for leaf in ete_tree.iter_leaves():
            leaf_distances.add(round(root.get_distance(leaf.name), 3))
        if len(leaf_distances) != 1:
            raise ValueError("Tree is not ultrametric")

        # Get tree height and calculate origin to root height
        tree_height = ete_tree.get_farthest_leaf()[1]
        origin_to_root_height = total_time - tree_height

        # Calculate metastasis times
        met_times = {}
        migrations = set()

        for node in ete_tree.traverse("levelorder"):
            if node.is_root():
                parent_tissue = primary_tissue
                parent_time = 0  # origin is at start of experiment
                node_time = origin_to_root_height
            else:
                parent_node = node.up
                parent_tissue = parent_node.name.split("_")[-1]
                root = ete_tree.get_tree_root()
                parent_time = origin_to_root_height + root.get_distance(
                    parent_node.name
                )
                node_time = origin_to_root_height + root.get_distance(node.name)

            node_tissue = node.name.split("_")[-1]

            if node_tissue != parent_tissue:
                migration = f"{parent_tissue}_{node_tissue}"
                migration_time = (parent_time, node_time)
                if migration not in migrations:
                    migrations.add(migration)
                    migration = migration + "_1"
                    met_times[migration] = migration_time
                else:
                    existing_migrations = [
                        key for key in met_times.keys() if migration in key
                    ]
                    i = (
                        max([int(key.split("_")[-1]) for key in existing_migrations])
                        + 1
                    )
                    migration = migration + "_" + str(i)
                    met_times[migration] = migration_time

        # Store results with tree label
        all_met_events[tree.label] = met_times

    # Save results if output_prefix is provided
    if output_prefix:
        output_file = f"{output_prefix}.pkl"
        with open(output_file, "wb") as f:
            pickle.dump(all_met_events, f)

    return all_met_events
