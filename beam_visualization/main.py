"""Main class for BEAM results analysis."""

from typing import Optional, Dict, List, Union, Tuple
import dendropy
import numpy as np
import pandas as pd

from . import data_loader
from . import plotting
from . import formatting
from . import config
from . import statistics


class BeamResults:
    """A class to handle BEAM output visualization and analysis."""

    def __init__(
        self, trees_file: str, log_file: str, primary_tissue: Optional[str] = None
    ):
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

    def get_parameter_stats(self, parameter: str) -> Dict[str, float]:
        """
        Get summary statistics for a parameter.

        Args:
            parameter: Name of the parameter to analyze

        Returns:
            Dict[str, float]: Dictionary containing summary statistics
        """
        if parameter not in self.log_data.columns:
            raise ValueError(f"Parameter '{parameter}' not found in log file")

        return {
            "mean": float(self.log_data[parameter].mean()),
            "std": float(self.log_data[parameter].std()),
            "min": float(self.log_data[parameter].min()),
            "max": float(self.log_data[parameter].max()),
            "median": float(self.log_data[parameter].median()),
            "q1": float(self.log_data[parameter].quantile(0.25)),
            "q3": float(self.log_data[parameter].quantile(0.75)),
        }

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
        plotting.plot_parameters(self.log_data, parameter, output_file)

    def get_consensus_graph(
        self,
        primary_tissue: Optional[str] = None,
        burnin_percent: float = config.DEFAULT_BURNIN_PERCENT,
        cores: int = config.DEFAULT_CORES,
        output_file: Optional[str] = None,
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
            self.trees, primary_tissue or self.primary_tissue, burnin_percent, cores
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
        self, output_file: Optional[str] = None, primary_tissue: Optional[str] = None
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
            raise ValueError(
                "Primary tissue must be provided either during initialization or when calling plot_probability_graph"
            )

        if self.consensus_graph is None:
            # Compute consensus graph if not already available
            self.get_consensus_graph(
                primary_tissue=primary_tissue or self.primary_tissue
            )

        plotting.plot_probability_graph(
            data=self.log_data,
            primary_tissue=primary_tissue or self.primary_tissue,
            output_file=output_file,
            consensus_graph=self.consensus_graph,
        )

    def plot_thresholded_graph(
        self,
        threshold: Union[float, List[float]] = 0.5,
        output_file_prefix: Optional[str] = None,
        primary_tissue: Optional[str] = None,
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
            raise ValueError(
                "Primary tissue must be provided either during initialization or when calling plot_thresholded_graph"
            )

        if self.consensus_graph is None:
            # Compute consensus graph if not already available
            self.get_consensus_graph(
                primary_tissue=primary_tissue or self.primary_tissue
            )

        plotting.plot_thresholded_graph(
            data=self.log_data,
            primary_tissue=primary_tissue or self.primary_tissue,
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

        Raises:
            ValueError: If primary_tissue is not set
        """
        if self.primary_tissue is None:
            raise ValueError(
                "Primary tissue must be set before computing mutual information"
            )

        if self.trees is None:
            raise ValueError("Trees must be loaded before computing mutual information")

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
        total_time: float = 1.0,
        output_prefix: Optional[str] = None,
        burnin_percent: float = 0.1,
    ) -> None:
        """
        Sample trees from the posterior and plot them with their migration graphs.

        Args:
            n: Number of trees to sample
            total_time: Total time of the tree
            output_prefix: Optional prefix for output files. If None, plots are displayed instead of saved.
            burnin_percent: Percentage of trees to discard as burnin

        Raises:
            ValueError: If primary_tissue is not set
        """
        if self.primary_tissue is None:
            raise ValueError("Primary tissue must be set before plotting sampled trees")

        if self.trees is None:
            raise ValueError("Trees must be loaded before plotting sampled trees")

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
                    total_time=total_time,
                    output_prefix=output_prefix,
                    tree_num=i,
                )

    def get_metastasis_times(
        self,
        total_time: float,
        primary_tissue: Optional[str] = None,
        burnin_percent: float = 0.1,
        min_prob_threshold: float = 0.5,
        output_prefix: Optional[str] = None,
    ) -> Dict[str, Dict[str, Tuple[float, float]]]:
        """
        Calculate metastasis times for all trees in the posterior distribution.

        Args:
            total_time: Total time of the experiment
            primary_tissue: Optional tissue label for the primary site. If not provided, uses the global primary_tissue.
            burnin_percent: Percentage of trees to discard as burnin
            min_prob_threshold: Minimum probability threshold for migrations to include in plots
            output_prefix: Optional prefix for output files. If provided:
                - Metastasis times will be written to {output_prefix}.pkl
                - Plots will be saved with the given prefix

        Returns:
            Dict[str, Dict[str, Tuple[float, float]]]: Dictionary mapping tree labels to dictionaries of metastasis events.
            Each metastasis event is a tuple of (start_time, end_time) for the migration.
        """

        if primary_tissue is None:
            if self.primary_tissue is None:
                raise ValueError("Primary tissue must be provided")
            else:
                primary_tissue = self.primary_tissue

        # Get metastasis times
        met_times = formatting.get_all_posterior_metastasis_times(
            self.trees,
            total_time=total_time,
            primary_tissue=primary_tissue,
            burnin_percent=burnin_percent,
            output_prefix=output_prefix,
        )

        # Get consensus graph if not already computed
        if self.consensus_graph is None:
            self.get_consensus_graph()

        # Plot metastasis timing
        plotting.plot_metastasis_timing(
            met_times=met_times,
            consensus_graph=self.consensus_graph,
            origin_time=total_time,
            origin_tissue=(
                primary_tissue if primary_tissue is not None else self.primary_tissue
            ),
            min_prob_threshold=min_prob_threshold,
            output_prefix=output_prefix,
        )

        return met_times
