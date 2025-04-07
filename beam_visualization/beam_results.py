"""Main class for BEAM results analysis."""

from typing import Optional, Dict, List
import dendropy
import pandas as pd
from .data_loader import load_beam_files
from .plotting import plot_parameters
from .formatting import get_consensus_graph
from .config import DEFAULT_BURNIN_PERCENT, DEFAULT_CORES

class BeamResults:
    """A class to handle BEAM output visualization and analysis."""
    
    def __init__(self, trees_file: str, log_file: str, primary_tissue: Optional[str] = None):
        """
        Initialize BeamResults with BEAM output files.
        
        Args:
            trees_file: Path to the .trees file
            log_file: Path to the .log file
            primary_tissue: Primary tissue label for migration analysis
        """
        self.trees_file = trees_file
        self.log_file = log_file
        self.primary_tissue = primary_tissue
        self.trees, self.log_data = load_beam_files(trees_file, log_file)
        self.consensus_graph = None
        
        # Validate data
        if len(self.trees) == 0:
            raise ValueError("No trees loaded")
        if len(self.log_data) == 0:
            raise ValueError("No log data loaded")
    
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
            parameter: Name of the parameter
            
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
            parameter: Specific parameter to plot. If None, plots all parameters.
            output_file: Path to save the plot. If None, displays the plot.
        """
        plot_parameters(self.log_data, parameter, output_file)
    
    def get_consensus_graph(
        self,
        primary_tissue: Optional[str] = None,
        burnin_percent: float = DEFAULT_BURNIN_PERCENT,
        cores: int = DEFAULT_CORES
    ) -> Dict[str, float]:
        """
        Calculate consensus graph from migration counts in trees.
        
        Args:
            primary_tissue: Primary tissue label for migration analysis
            burnin_percent: Percentage of trees to discard as burnin
            cores: Number of CPU cores to use for parallel processing
            
        Returns:
            Dict[str, float]: Dictionary mapping migration patterns to their probabilities
        """
        self.consensus_graph = get_consensus_graph(
            self.trees,
            primary_tissue or self.primary_tissue,
            burnin_percent,
            cores
        )
        return self.consensus_graph 