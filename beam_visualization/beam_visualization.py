"""Main class for BEAM results analysis."""

from typing import Optional, Dict, List, Union, Tuple
import dendropy
import numpy as np
import pandas as pd
import os

from . import data_loader
from . import plotting
from . import formatting
from . import config
from . import statistics


class BeamResults:
    """A class to handle BEAM output visualization and analysis."""

    def _ensure_output_dir(self, filepath: Optional[str]) -> None:
        """
        Ensure the directory for the given filepath exists.
        
        Args:
            filepath: Path to the file to be written
        """
        if filepath:
            output_dir = os.path.dirname(filepath)
            if output_dir and not os.path.exists(output_dir):
                os.makedirs(output_dir, exist_ok=True)

    def __init__(
        self, 
        trees_file: str, 
        log_file: str, 
        primary_tissue: str,
        total_time: float
    ):
        """
        Initialize BeamResults with BEAM output files.

        Args:
            trees_file: Path to the .trees file
            log_file: Path to the .log file
            primary_tissue: Primary tissue label for migration analysis
            total_time: Total time of the experiment
        """
        self.trees_file = trees_file
        self.log_file = log_file
        self.primary_tissue = primary_tissue
        self.total_time = total_time
        self.trees, self.log_data = data_loader.load_beam_files(trees_file, log_file)
        self.consensus_graph = None
        self.metastasis_times = None

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
        print(f"Total time: {self.total_time}")
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
                self.consensus_graph.items(), key=lambda x: x[1], reverse=True
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
            "mean": self.log_data[parameter].mean(),
            "std": self.log_data[parameter].std(),
            "min": self.log_data[parameter].min(),
            "max": self.log_data[parameter].max(),
            "median": self.log_data[parameter].median(),
            "q1": self.log_data[parameter].quantile(0.25),
            "q3": self.log_data[parameter].quantile(0.75),
        }
        return stats

    def get_trees(self) -> dendropy.TreeList:
        """
        Get the loaded trees.

        Returns:
            dendropy.TreeList: The loaded trees
        """
        return self.trees

    def plot_parameters(
        self, parameter: Optional[str] = None, output_file: Optional[str] = None
    ) -> None:
        """
        Plot parameter distributions from the log file.

        Args:
            parameter: Specific parameter to plot. If None, plots all parameters.
            output_file: Path to save the plot. If None, displays the plot.
        """

        if output_file:
            self._ensure_output_dir(output_file)

        plotting.plot_parameters(self.log_data, parameter, output_file)

    def get_consensus_graph(
        self,
        burnin_percent: float = config.DEFAULT_BURNIN_PERCENT,
        cores: int = config.DEFAULT_CORES,
        output_file: Optional[str] = None,
        force_recompute: bool = False,
    ) -> Dict[str, float]:
        """
        Calculate consensus graph from migration counts in trees.

        Args:
            burnin_percent: Percentage of trees to discard as burnin
            cores: Number of CPU cores to use for parallel processing
            output_file: Optional path to write the consensus graph to a file
            force_recompute: If True, forces recalculation of consensus graph even if already computed

        Returns:
            Dict[str, float]: Dictionary mapping migration patterns to their probabilities
        """
        # Only compute consensus graph if not already computed or if force_recompute is True
        if self.consensus_graph is None or force_recompute:
            self.consensus_graph = formatting.get_consensus_graph(
                self.trees, self.primary_tissue, burnin_percent, cores
            )

        if output_file:
            self._ensure_output_dir(output_file)
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
        output_file: Optional[str] = None
    ) -> None:
        """
        Plot the consensus migration graph with edge thicknesses proportional to probability.

        Args:
            output_file: Optional path to save the plot. If None, displays the plot.
        """
        if self.consensus_graph is None:
            # Compute consensus graph if not already available
            self.get_consensus_graph()

        if output_file:
            self._ensure_output_dir(output_file)

        plotting.plot_probability_graph(
            data=self.log_data,
            primary_tissue=self.primary_tissue,
            output_file=output_file,
            consensus_graph=self.consensus_graph,
        )

    def plot_thresholded_graph(
        self,
        threshold: Union[float, List[float]] = 0.5,
        output_file_prefix: Optional[str] = None,
    ) -> None:
        """
        Plot the thresholded consensus migration graph with collapsed multiedges.

        Args:
            threshold: Single threshold value or list of thresholds for filtering edges
            output_file_prefix: Optional prefix for output files. If None, displays the plot(s).
                              For multiple thresholds, each file will be named {prefix}_{threshold}.pdf
        """
        if self.consensus_graph is None:
            # Compute consensus graph if not already available
            self.get_consensus_graph()

        if output_file_prefix:
            self._ensure_output_dir(output_file_prefix)

        plotting.plot_thresholded_graph(
            data=self.log_data,
            primary_tissue=self.primary_tissue,
            threshold=threshold,
            output_file_prefix=output_file_prefix,
            consensus_graph=self.consensus_graph,
        )

    def compute_posterior_mutual_info(
        self,
        output_file_matrix: Optional[str] = None,
        output_file_information: Optional[str] = None,
        threads: int = 1,
    ) -> Tuple[float, np.ndarray, List[str]]:
        """Calculate mutual information from migration patterns in the trees.

        Args:
            output_file_matrix: Optional path to save the count matrix
            output_file_information: Optional path to save the mutual information score
            threads: Number of threads to use for parallel processing

        Returns:
            Tuple containing:
                - Mutual information score
                - Count matrix as numpy array
                - List of tissue names in order
        """
        if self.trees is None:
            raise ValueError("Trees must be loaded before computing mutual information")

        if output_file_matrix:
            self._ensure_output_dir(output_file_matrix)
        if output_file_information:
            self._ensure_output_dir(output_file_information)

        return statistics.compute_posterior_mutual_info(
            trees=self.trees,
            origin_tissue=self.primary_tissue,
            threads=threads,
            output_file_matrix=output_file_matrix,
            output_file_information=output_file_information,
        )

    def sample_and_plot_trees(
        self,
        n: int = 1,
        output_prefix: Optional[str] = None,
        burnin_percent: float = config.DEFAULT_BURNIN_PERCENT,
    ) -> None:
        """
        Sample trees from the posterior and plot them with their migration graphs.

        Args:
            n: Number of trees to sample
            output_prefix: Optional prefix for output files. If None, plots are displayed instead of saved.
            burnin_percent: Percentage of trees to discard as burnin
        """
        if self.trees is None:
            raise ValueError("Trees must be loaded before plotting sampled trees")

        if output_prefix:
            self._ensure_output_dir(output_prefix)

        # Sample trees
        sampled_trees = formatting.sample_trees(
            self.trees, n=n, burnin_percent=burnin_percent, output_prefix=output_prefix
        )

        if output_prefix:
            # Plot each tree
            for i, tree in enumerate(sampled_trees, 1):
                plotting.plot_sampled_tree(
                    newick_str=tree,
                    primary_tissue=self.primary_tissue,
                    total_time=self.total_time,
                    output_prefix=output_prefix,
                    tree_num=i,
                )

    def get_metastasis_times(
        self,
        burnin_percent: float = config.DEFAULT_BURNIN_PERCENT,
        min_prob_threshold: float = 0.5,
        output_prefix: Optional[str] = None,
        force_recompute: bool = False,
    ) -> Dict[str, Dict[str, Tuple[float, float]]]:
        """
        Calculate metastasis times for all trees in the posterior distribution.

        Args:
            burnin_percent: Percentage of trees to discard as burnin
            min_prob_threshold: Minimum probability threshold for migrations to include in plots
            output_prefix: Optional prefix for output files. If provided:
                - Metastasis times will be written to {output_prefix}.pkl
                - Plots will be saved with the given prefix
            force_recompute: If True, forces recalculation of metastasis times even if already computed

        Returns:
            Dict[str, Dict[str, Tuple[float, float]]]: Dictionary mapping tree labels to dictionaries of metastasis events.
            Each metastasis event is a tuple of (start_time, end_time) for the migration.
        """
        if output_prefix:
            self._ensure_output_dir(output_prefix)

        # Only compute metastasis times if not already computed or if force_recompute is True
        if self.metastasis_times is None or force_recompute:
            # Get metastasis times
            self.metastasis_times = formatting.get_all_posterior_metastasis_times(
                self.trees,
                total_time=self.total_time,
                primary_tissue=self.primary_tissue,
                burnin_percent=burnin_percent,
                output_prefix=output_prefix,
            )

        # Get consensus graph if not already computed
        if self.consensus_graph is None:
            self.get_consensus_graph()

        # Plot metastasis timing
        plotting.plot_metastasis_timing(
            met_times=self.metastasis_times,
            consensus_graph=self.consensus_graph,
            origin_time=self.total_time,
            origin_tissue=self.primary_tissue,
            min_prob_threshold=min_prob_threshold,
            output_prefix=output_prefix,
        )

        return self.metastasis_times
