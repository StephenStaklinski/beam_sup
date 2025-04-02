import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Optional, Dict, List, Union
import dendropy

class BeamResults:
    """A class to handle BEAM output visualization and analysis."""
    
    def __init__(self, trees_file: str, log_file: str):
        """
        Initialize BeamResults with BEAM output files.
        
        Args:
            trees_file (str): Path to the .trees file
            log_file (str): Path to the .log file
        """
        print(f"Initializing BeamResults with files:")
        print(f"  Trees file: {trees_file}")
        print(f"  Log file: {log_file}")
        self.trees_file = trees_file
        self.log_file = log_file
        self.trees = None
        self.log_data = None
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
        print(f"\nTrees (self.trees):")
        print(f"  Number of trees: {len(self.trees)}")
        print(f"  Taxa: {', '.join(self.trees.taxon_namespace.labels())}")
        print(f"\nLog Data (self.log_data):")
        print(f"  Number of samples: {len(self.log_data)}")
        print(f"  Parameters: {', '.join(self.log_data.columns)}")
        print("\nAvailable Methods:")
        print("  - info(): Show this information")
        print("  - get_parameters(): List available parameters")
        print("  - get_parameter_stats(parameter): Get statistics for a parameter")
        print("  - plot_parameters(): Plot parameter distributions")
        print("  - plot_migration_graphs(): Plot migration graphs")
        print("  - plot_rate_matrix(): Plot rate matrix")
    
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
    

    