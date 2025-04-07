import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Optional, Dict, List, Union
import dendropy
from scipy.stats import gaussian_kde
from ete3 import Tree
from copy import deepcopy
from multiprocessing import Pool

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

class BeamResults:
    """A class to handle BEAM output visualization and analysis."""
    
    def __init__(self, trees_file: str, log_file: str, primary_tissue: Optional[str] = None):
        """
        Initialize BeamResults with BEAM output files.
        
        Args:
            trees_file (str): Path to the .trees file
            log_file (str): Path to the .log file
            primary_tissue (str, optional): Primary tissue label for migration analysis
        """
        print(f"Initializing BeamResults with files:")
        print(f"  Trees file: {trees_file}")
        print(f"  Log file: {log_file}")
        self.trees_file = trees_file
        self.log_file = log_file
        self.primary_tissue = primary_tissue
        self.trees = None
        self.log_data = None
        self.consensus_graph = None
        self._load_data()
    
    def _load_data(self):
        """Load and parse the trees and log files."""
        print("\nLoading trees file...")
        # Load trees (assuming NEXUS format)
        try:
            print("  Attempting to load as NEXUS format...")
            self.trees = dendropy.TreeList.get(path=self.trees_file, schema="nexus")
            print(f"  Successfully loaded {len(self.trees)} trees")
        except Exception as e:
            print(f"  Initial NEXUS parsing failed: {e}")
            print("  Attempting to fix case sensitivity...")
            # If NEXUS parsing fails, try to read the file and fix case sensitivity
            with open(self.trees_file, 'r') as f:
                content = f.read()
            # Replace case-sensitive NEXUS keywords
            content = content.replace('Begin taxa', 'BEGIN TAXA')
            content = content.replace('Dimensions', 'DIMENSIONS')
            content = content.replace('Taxlabels', 'TAXLABELS')
            content = content.replace('Begin trees', 'BEGIN TREES')
            # Write to temporary file
            temp_file = self.trees_file + '.temp'
            with open(temp_file, 'w') as f:
                f.write(content)
            print("  Created temporary file with fixed case sensitivity")
            # Try to load the fixed file
            self.trees = dendropy.TreeList.get(path=temp_file, schema="nexus")
            print(f"  Successfully loaded {len(self.trees)} trees after case sensitivity fix")
            # Clean up temporary file
            os.remove(temp_file)
            print("  Cleaned up temporary file")
        
        print("\nLoading log file...")
        # Load log data
        try:
            print("  Reading log file...")
            self.log_data = pd.read_csv(self.log_file, sep="\t", comment="#")
            print(f"  Successfully loaded log file with {len(self.log_data)} rows and {len(self.log_data.columns)} columns")
            print("  Available parameters:", ", ".join(self.log_data.columns))
        except Exception as e:
            print(f"  Error loading log file: {e}")
            raise
    
    
    def info(self) -> None:
        """Print information about the loaded data."""
        print("\nBeamResults Object Information:")
        print("===============================")
        print(f"Trees file: {self.trees_file}")
        print(f"Log file: {self.log_file}")
        print(f"Primary tissue: {self.primary_tissue}")
        print(f"\nTrees:")
        print(f"  Number of trees: {len(self.trees)}")
        print(f"  Taxa: {', '.join(self.trees.taxon_namespace.labels())}")
        print(f"\nLog Data:")
        print(f"  Number of samples: {len(self.log_data)}")
        print(f"  Parameters: {', '.join(self.log_data.columns)}")
        if self.consensus_graph is not None:
            print(f"\nMigration Analysis:")
            print(f"  Number of migrations: {len(self.consensus_graph)}")
            top_migrations = sorted(
                self.consensus_graph.items(),
                key=lambda x: x[1],
                reverse=True
            )[:5]
            print("  Top migrations:")
            for migration, prob in top_migrations:
                print(f"    {migration}: {prob:.4f}")
        print("\nAvailable Methods:")
        print("  - info(): Show this information")
        print("  - get_parameters(): List available parameters")
        print("  - get_parameter_stats(parameter): Get statistics for a parameter")
        print("  - get_trees(): Get the loaded trees")
        print("  - plot_parameters(): Plot parameter distributions")
        print("  - get_consensus_graph(): Calculate consensus graph from migration counts")
    
    def get_parameters(self) -> List[str]:
        """
        Get a list of available parameters from the log file.
        
        Returns:
            List[str]: List of parameter names
        """
        return list(self.log_data.columns)
    
    def get_parameter_stats(self, parameter: str) -> Dict:
        """
        Get statistics for a specific parameter.
        
        Args:
            parameter (str): Name of the parameter
            
        Returns:
            Dict: Dictionary containing parameter statistics
        """
        if parameter not in self.log_data.columns:
            raise ValueError(f"Parameter '{parameter}' not found in log file")
        
        stats = {
            'mean': self.log_data[parameter].mean(),
            'std': self.log_data[parameter].std(),
            'min': self.log_data[parameter].min(),
            'max': self.log_data[parameter].max(),
            'median': self.log_data[parameter].median(),
            'q1': self.log_data[parameter].quantile(0.25),
            'q3': self.log_data[parameter].quantile(0.75)
        }
        return stats
    
    def get_trees(self) -> dendropy.TreeList:
        """
        Get the loaded trees.
        
        Returns:
            dendropy.TreeList: The loaded trees
        """
        return self.trees
    
    def plot_parameters(self, parameter: Optional[str] = None, output_file: Optional[str] = None) -> None:
        """
        Plot parameter distributions from the log file.
        
        Args:
            parameter (str, optional): Specific parameter to plot. If None, plots all parameters.
            output_file (str, optional): Path to save the plot. If None, displays the plot.
        """
        if parameter:
            if parameter not in self.log_data.columns:
                raise ValueError(f"Parameter '{parameter}' not found in log file")
            data = self.log_data[parameter]
        else:
            data = self.log_data
        
        plt.figure(figsize=(12, 6))
        sns.histplot(data=data, kde=False)
        plt.title(f"Distribution of {parameter if parameter else 'Parameters'}")
        plt.tight_layout()
        
        if output_file:
            plt.savefig(output_file)
            plt.close()
        else:
            plt.show()
    
    def get_consensus_graph(self, primary_tissue: Optional[str] = None, burnin_percent: float = 0.0, cores: int = 1) -> Dict[str, float]:
        """
        Calculate consensus graph from migration counts in trees.
        
        Args:
            primary_tissue (str, optional): Primary tissue label for migration analysis. If None, uses the one set during initialization.
            burnin_percent (float): Percentage of trees to discard as burnin. Default is 0.0 (no burnin).
            cores (int): Number of CPU cores to use for parallel processing. Default is 1 (single core).
            
        Returns:
            Dict[str, float]: Dictionary mapping migration patterns to their probabilities.
            Each key is in the format "source_tissue_dest_tissue_count" and the value is the probability
            of that migration pattern occurring in the posterior distribution.
            
        Raises:
            ValueError: If primary_tissue is not specified either during initialization or as an argument.
        """
        # Use provided primary_tissue or the one from initialization
        if primary_tissue is None:
            if self.primary_tissue is None:
                raise ValueError("Primary tissue must be specified either during initialization or as an argument")
            primary_tissue = self.primary_tissue
        
        print("\nCalculating consensus graph...")
        
        # Process posterior trees
        num_trees = len(self.trees)
        num_discard = round(num_trees * burnin_percent)
        trees_to_analyze = self.trees[num_discard:]
        
        print(f"  Analyzing {len(trees_to_analyze)} trees (after {num_discard} burnin)")
        
        # Process trees in parallel
        with Pool(processes=cores) as pool:
            all_counts = pool.map(
                _process_tree_wrapper,
                [(tree, primary_tissue) for tree in trees_to_analyze]
            )
        
        # Calculate consensus graph
        prob = 1 / len(all_counts)
        self.consensus_graph = {}
        for counts in all_counts:
            for migration, count in counts.items():
                for i in range(1, count + 1):
                    edge = f"{migration}_{i}"
                    self.consensus_graph[edge] = self.consensus_graph.get(edge, 0) + prob
        
        print("  Consensus graph calculation complete")
        return self.consensus_graph
    

    