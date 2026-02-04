import argparse
from typing import Optional, Dict, List, Union, Tuple
from dendropy import TreeList
import os


DEFAULT_BURNIN = 0.0
DEFAULT_CORES = 1
DEFAULT_MIN_PROB_THRESHOLD = 0.5


class BeamResults:
    """A class to handle migration graph posterior distribution analysis."""

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
        self, trees_file: str, primary_tissue: str, total_time: float, log_file: str = None, cores: int = DEFAULT_CORES, state_key: str = "location", burnin: float = DEFAULT_BURNIN
    ) -> None:
        """
        Initialize BeamResults with BEAM output files.

        Args:
            trees_file: Path to the .trees file
            primary_tissue: Primary tissue label for migration analysis
            total_time: Total time of the experiment
            log_file: Path to the .log file (optional)
            state_key: Key in the node annotations for tissue states. Default is "location".
            cores: Number of CPU threads to use for parallel processing
        """
        from .data_loader import load_beam_files
        
        self.trees_file = trees_file
        self.log_file = log_file
        self.primary_tissue = primary_tissue
        self.total_time = total_time
        self.trees, self.log_data = load_beam_files(trees_file, log_file)
        self.consensus_graph = None
        self.metastasis_times = None
        self.cores = cores
        self.state_key = state_key
        self.burnin = burnin

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
        print(f"  Number of taxa: {len(self.trees.taxon_namespace.labels())}")
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

    def get_trees(self) -> TreeList:
        """
        Get the loaded trees.

        Returns:
            TreeList: The loaded trees
        """
        return self.trees

    def plot_parameters(
        self, output_file: str, parameter: Optional[str] = None
    ) -> None:
        """
        Plot parameter distributions from the log file.

        Args:
            output_file: Path to save the plot.
            parameter: Specific parameter to plot. If None, plots all parameters.
        """
        from .plotting import plot_parameters
        self._ensure_output_dir(output_file)
        plot_parameters(self.log_data, output_file, parameter)

    def get_consensus_graph(
        self,
        output_file: Optional[str] = None,
        force_recompute: bool = False,
    ) -> Dict[str, float]:
        """
        Calculate consensus graph from migration counts in trees.

        Args:
            output_file: Optional path to write the consensus graph to a file
            force_recompute: If True, forces recalculation of consensus graph even if already computed

        Returns:
            Dict[str, float]: Dictionary mapping migration patterns to their probabilities
        """
        from .posterior_processing import get_consensus_graph

        if self.trees is None or len(self.trees) == 0:
            raise ValueError("No trees to analyze")

        # Only compute consensus graph if not already computed or if force_recompute is True
        if self.consensus_graph is None or force_recompute:
            self.consensus_graph = get_consensus_graph(
                self.trees, self.primary_tissue, self.burnin, self.cores, self.state_key
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

    def plot_probability_graph(self, output_file: str) -> None:
        """
        Plot the consensus migration graph with edge thicknesses proportional to probability.

        Args:
            output_file: Path to save the plot.
        """
        from .plotting import plot_probability_graph
        if self.consensus_graph is None:
            # Compute consensus graph if not already available
            self.get_consensus_graph()

        self._ensure_output_dir(output_file)

        plot_probability_graph(
            data=self.log_data,
            output_file=output_file,
            primary_tissue=self.primary_tissue,
            consensus_graph=self.consensus_graph,
        )

    def plot_thresholded_graph(
        self,
        output_file_prefix: str,
        threshold: Union[float, List[float]] = 0.5,
    ) -> None:
        """
        Plot the thresholded consensus migration graph with collapsed multiedges.

        Args:
            output_file_prefix: Prefix for output files. Each file will be named {prefix}_{threshold}.pdf
            threshold: Single threshold value or list of thresholds for filtering edges
        """
        from .plotting import plot_thresholded_graph
        if self.consensus_graph is None:
            # Compute consensus graph if not already available
            self.get_consensus_graph()

        self._ensure_output_dir(output_file_prefix)

        plot_thresholded_graph(
            data=self.log_data,
            output_file_prefix=output_file_prefix,
            primary_tissue=self.primary_tissue,
            threshold=threshold,
            consensus_graph=self.consensus_graph,
        )

    def compute_posterior_mutual_info(
        self,
        output_file_matrix: Optional[str] = None,
        output_file_information: Optional[str] = None,
    ) -> Tuple[float, 'np.ndarray', List[str]]:
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
        from .posterior_processing import compute_posterior_mutual_info
        
        if self.trees is None:
            raise ValueError("Trees must be loaded before computing mutual information")

        if output_file_matrix:
            self._ensure_output_dir(output_file_matrix)
        if output_file_information:
            self._ensure_output_dir(output_file_information)

        return compute_posterior_mutual_info(
            trees=self.trees,
            primary_tissue=self.primary_tissue,
            cores=self.cores,
            output_file_matrix=output_file_matrix,
            output_file_information=output_file_information,
            state_key=self.state_key,
        )

    def sample_and_plot_trees(
        self,
        n: int = 1,
        output_prefix: str = None,
    ) -> None:
        """
        Sample trees from the posterior and plot them with their migration graphs.

        Args:
            n: Number of trees to sample
            output_prefix: Prefix for output files.
        """
        from .posterior_processing import sample_trees
        from .plotting import plot_sampled_tree
        
        self._ensure_output_dir(output_prefix)

        # Sample trees
        sampled_trees = sample_trees(
            self.trees, n=n, burnin_percent=self.burnin, output_prefix=output_prefix, state_key=self.state_key
        )

        # Plot each tree
        for i, tree in enumerate(sampled_trees, 1):
            plot_sampled_tree(
                newick_str=tree,
                primary_tissue=self.primary_tissue,
                total_time=self.total_time,
                output_prefix=output_prefix,
                tree_num=i,
                state_key=self.state_key
            )

    def get_metastasis_times(
        self,
        min_prob_threshold: float = DEFAULT_MIN_PROB_THRESHOLD,
        output_prefix: str = None,
        force_recompute: bool = False,
    ) -> Dict[str, Dict[str, Tuple[float, float]]]:
        """
        Calculate metastasis times for all trees in the posterior distribution.

        Args:
            min_prob_threshold: Minimum probability threshold for migrations to include in plots
            output_prefix: Prefix for output files:
                - Metastasis times will be written to {output_prefix}.pkl
                - Plots will be saved with the given prefix
            force_recompute: If True, forces recalculation of metastasis times even if already computed

        Returns:
            Dict[str, Dict[str, Tuple[float, float]]]: Dictionary mapping tree labels to dictionaries of metastasis events.
            Each metastasis event is a tuple of (start_time, end_time) for the migration.
        """
        from .posterior_processing import get_all_posterior_metastasis_times
        from .plotting import plot_metastasis_timing
        
        self._ensure_output_dir(output_prefix)

        # Only compute metastasis times if not already computed or if force_recompute is True
        if self.metastasis_times is None or force_recompute:
            # Get metastasis times
            self.metastasis_times = get_all_posterior_metastasis_times(
                self.trees,
                total_time=self.total_time,
                primary_tissue=self.primary_tissue,
                burnin_percent=self.burnin,
                output_prefix=output_prefix,
                cores=self.cores,
                state_key=self.state_key
            )

        # Get consensus graph if not already computed
        if self.consensus_graph is None:
            self.get_consensus_graph()

        # Plot metastasis timing
        plot_metastasis_timing(
            met_times=self.metastasis_times,
            consensus_graph=self.consensus_graph,
            origin_time=self.total_time,
            primary_tissue=self.primary_tissue,
            min_prob_threshold=min_prob_threshold,
            output_prefix=output_prefix,
        )

        return self.metastasis_times

    def plot_rate_matrix(self, output_file: str) -> None:
        """
        Plot the tissue substitution rate matrix as a heatmap.

        Args:
            output_file: Path to save the plot
        """
        from .plotting import plot_rate_matrix
        self._ensure_output_dir(output_file)
        plot_rate_matrix(
            data=self.log_data,
            output_file=output_file,
            primary_tissue=self.primary_tissue,
        )


def run_full_simulation():
    from .simulate import simulate_metastatic_cancer_population, overlay_simulated_crispr_barcode_data

    parser = argparse.ArgumentParser(
        description="Run metastatic cancer population simulation and overlay CRISPR barcode data."
    )
    parser.add_argument(
        "--outdir", default="./", help="Output directory for simulation results. Must be unique since reruns will remove existing files in this directory."
    )
    parser.add_argument(
        "--outprefix", default="test", help="Output prefix for simulation results. Will be written in the outdir specified."
    )
    parser.add_argument(
        "--num_generations",
        type=int,
        default=250,
        help="Number of generations to simulate.",
    )
    parser.add_argument(
        "--migration_rate",
        type=str,
        default="1e-6",
        help="Migration rate between anatomical sites.",
    )
    parser.add_argument(
        "--num_cells_downsample",
        type=int,
        default=50,
        help="Number of cells to downsample.",
    )
    parser.add_argument(
        "--num_possible_tissues",
        type=int,
        default=10,
        help="Number of tissues to simulate up to.",
    )
    parser.add_argument(
        "--max_anatomical_sites",
        type=int,
        default=-1,
        help="Maximum number of anatomical sites.",
    )
    parser.add_argument(
        "--migration_start_generation",
        type=int,
        default=0,
        help="Generation at which migrations can start occurring.",
    )
    parser.add_argument(
        "--migration_end_generation",
        type=int,
        default=-1,
        help="Generation after which migrations stop occurring.",
    )
    parser.add_argument(
        "--seed", type=int, default=None, help="Random seed for reproducibility."
    )
    parser.add_argument(
        "--num_sites",
        type=int,
        default=50,
        help="Number of CRISPR cassettes (cuts) to simulate.",
    )
    parser.add_argument(
        "--mutationrate", type=float, default=0.1, help="Mutation rate per cassette."
    )
    parser.add_argument(
        "--heritable_silencing_rate",
        type=float,
        default=0.0001,
        help="Rate of heritable silencing of CRISPR cassettes.",
    )
    parser.add_argument(
        "--stochastic_silencing_rate",
        type=float,
        default=0.01,
        help="Rate of stochastic silencing of CRISPR cassettes.",
    )
    parser.add_argument(
        "--transition_matrix_filepath",
        type=str,
        default=None,
        help="File path to a custom transition matrix CSV file.",
    )

    args = parser.parse_args()

    # Run the cancer population simulation
    seed = simulate_metastatic_cancer_population(
        outdir=args.outdir,
        outprefix=args.outprefix,
        num_generations=args.num_generations,
        migration_rate=args.migration_rate,
        num_cells_downsample=args.num_cells_downsample,
        num_possible_tissues=args.num_possible_tissues,
        max_anatomical_sites=args.max_anatomical_sites,
        migration_start_generation=args.migration_start_generation,
        migration_end_generation=args.migration_end_generation,
        transition_matrix_filepath=args.transition_matrix_filepath,
        seed=args.seed,
    )

    # Find the ground truth tree file
    outputdir = args.outdir
    ground_truth_tree = None
    if os.path.isdir(outputdir):
        for fname in os.listdir(outputdir):
            if fname.endswith(".nwk") and fname.startswith(args.outprefix):
                ground_truth_tree = os.path.join(outputdir, fname)
                break
    if ground_truth_tree is None:
        raise FileNotFoundError(f"Ground truth tree file (.nwk) not found in {outputdir}.")

    # Overlay simulated CRISPR barcode data
    overlay_simulated_crispr_barcode_data(
        ground_truth_tree_filepath=ground_truth_tree,
        outprefix=args.outprefix,
        num_sites=args.num_sites,
        mutationrate=args.mutationrate,
        heritable_silencing_rate=args.heritable_silencing_rate,
        stochastic_silencing_rate=args.stochastic_silencing_rate
    )


def simulate_barcode_overlay_only():
    parser = argparse.ArgumentParser(
        description="Overlay simulated CRISPR barcode data onto a given ground truth tree."
    )
    parser.add_argument(
        "--ground_truth_tree",
        required=True,
        help="File path to the ground truth tree in Newick format.",
    )
    parser.add_argument(
        "--outprefix", required=True, help="Output prefix for the simulated data."
    )
    parser.add_argument(
        "--num_sites",
        type=int,
        default=50,
        help="Number of CRISPR cassettes (cuts) to simulate.",
    )
    parser.add_argument(
        "--mutationrate", type=float, default=0.1, help="Mutation rate per cassette."
    )
    parser.add_argument(
        "--heritable_silencing_rate",
        type=float,
        default=0.0001,
        help="Rate of heritable silencing of CRISPR cassettes.",
    )
    parser.add_argument(
        "--stochastic_silencing_rate",
        type=float,
        default=0.01,
        help="Rate of stochastic silencing of CRISPR cassettes.",
    )

    args = parser.parse_args()

    overlay_simulated_crispr_barcode_data(
        ground_truth_tree_filepath=args.ground_truth_tree,
        outprefix=args.outprefix,
        num_sites=args.num_sites,
        mutationrate=args.mutationrate,
        heritable_silencing_rate=args.heritable_silencing_rate,
        stochastic_silencing_rate=args.stochastic_silencing_rate
    )