"""Main class for BEAM results analysis."""

from typing import Optional, Dict, List, Union
import dendropy
import pandas as pd
from . import data_loader
from . import plotting
from . import formatting
from . import config

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
        self.trees, self.log_data = data_loader.load_beam_files(trees_file, log_file)
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
        plotting.plot_parameters(self.log_data, parameter, output_file)
    
    def get_consensus_graph(
        self,
        primary_tissue: Optional[str] = None,
        burnin_percent: float = config.DEFAULT_BURNIN_PERCENT,
        cores: int = config.DEFAULT_CORES,
        output_file: Optional[str] = None
    ) -> Dict[str, float]:
        """
        Calculate consensus graph from migration counts in trees.
        
        Args:
            primary_tissue: Primary tissue label for migration analysis
            burnin_percent: Percentage of trees to discard as burnin
            cores: Number of CPU cores to use for parallel processing
            output_file: Optional path to write the consensus graph to a file
            
        Returns:
            Dict[str, float]: Dictionary mapping migration patterns to their probabilities
        """
        self.consensus_graph = formatting.get_consensus_graph(
            self.trees,
            primary_tissue or self.primary_tissue,
            burnin_percent,
            cores
        )
        
        if output_file:
            # Sort the graph by probability in descending order
            sorted_graph = dict(
                sorted(self.consensus_graph.items(), key=lambda x: x[1], reverse=True)
            )
            # Write to file
            with open(output_file, "w") as file:
                for key, value in sorted_graph.items():
                    file.write(f"{key},{value}\n")
        
        return self.consensus_graph
    
    def plot_probability_graph(
        self,
        output_file: Optional[str] = None,
        primary_tissue: Optional[str] = None
    ) -> None:
        """
        Plot the consensus migration graph with edge thicknesses proportional to probability.
        
        Args:
            output_file: Optional path to save the plot. If None, displays the plot.
            primary_tissue: Primary tissue label for migration analysis. If None, uses the primary_tissue from initialization.
            
        Raises:
            ValueError: If primary_tissue is not provided either here or during initialization
        """
        if primary_tissue is None and self.primary_tissue is None:
            raise ValueError("Primary tissue must be provided either during initialization or when calling plot_probability_graph")
        
        if self.consensus_graph is None:
            # Compute consensus graph if not already available
            self.get_consensus_graph(primary_tissue=primary_tissue or self.primary_tissue)
        
        plotting.plot_probability_graph(
            data=self.log_data,
            primary_tissue=primary_tissue or self.primary_tissue,
            output_file=output_file,
            consensus_graph=self.consensus_graph
        )
    
    def plot_thresholded_graph(
        self,
        threshold: Union[float, List[float]] = 0.5,
        output_file_prefix: Optional[str] = None,
        primary_tissue: Optional[str] = None
    ) -> None:
        """
        Plot the thresholded consensus migration graph with collapsed multiedges.
        
        Args:
            threshold: Single threshold value or list of thresholds for filtering edges
            output_file_prefix: Optional prefix for output files. If None, displays the plot(s).
                              For multiple thresholds, each file will be named {prefix}_{threshold}.pdf
            primary_tissue: Primary tissue label for migration analysis. If None, uses the primary_tissue from initialization.
            
        Raises:
            ValueError: If primary_tissue is not provided either here or during initialization
        """
        if primary_tissue is None and self.primary_tissue is None:
            raise ValueError("Primary tissue must be provided either during initialization or when calling plot_thresholded_graph")
        
        if self.consensus_graph is None:
            # Compute consensus graph if not already available
            self.get_consensus_graph(primary_tissue=primary_tissue or self.primary_tissue)
        
        plotting.plot_thresholded_graph(
            data=self.log_data,
            primary_tissue=primary_tissue or self.primary_tissue,
            threshold=threshold,
            output_file_prefix=output_file_prefix,
            consensus_graph=self.consensus_graph
        )

