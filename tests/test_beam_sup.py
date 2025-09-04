import pytest
import os
import numpy as np
from beam_sup.beam_sup import BeamResults
from beam_sup import data_loader
import pandas as pd
import dendropy
import random

# Sample test data paths - these should be replaced with actual test data
TEST_TREES_FILE = "tests/data/test.trees"
TEST_LOG_FILE = "tests/data/test.log"

# Test data for error cases
NONEXISTENT_FILE = "tests/data/nonexistent.trees"
INVALID_TREES_FILE = "tests/data/invalid.trees"
INVALID_LOG_FILE = "tests/data/invalid.log"

@pytest.fixture
def beam_results():
    """Fixture to create a BeamResults instance for testing."""
    return BeamResults(
        trees_file=TEST_TREES_FILE,
        log_file=TEST_LOG_FILE,
        primary_tissue="LL",
        total_time=54.0
    )

def test_initialization(beam_results):
    """Test that BeamResults initializes correctly."""
    assert beam_results.trees_file == TEST_TREES_FILE
    assert beam_results.log_file == TEST_LOG_FILE
    assert beam_results.primary_tissue == "LL"
    assert beam_results.total_time == 54.0
    assert beam_results.trees is not None
    assert beam_results.log_data is not None

def test_get_parameters(beam_results):
    """Test that get_parameters returns a list of strings."""
    parameters = beam_results.get_parameters()
    assert isinstance(parameters, list)
    assert all(isinstance(param, str) for param in parameters)

def test_get_parameter_stats(beam_results):
    """Test that get_parameter_stats returns a dictionary with expected keys."""
    # Get the first parameter to test with
    parameters = beam_results.get_parameters()
    if parameters:
        stats = beam_results.get_parameter_stats(parameters[0])
        assert isinstance(stats, dict)
        assert "mean" in stats
        assert "std" in stats
        assert "min" in stats
        assert "max" in stats

def test_get_parameter_stats_empty_parameters(beam_results):
    """Test that get_parameter_stats handles empty parameters list."""
    parameters = beam_results.get_parameters()
    if not parameters:
        stats = beam_results.get_parameter_stats("nonexistent_parameter")
        assert stats == {}

def test_get_parameter_stats_invalid_parameter(beam_results):
    """Test that get_parameter_stats raises ValueError for invalid parameter."""
    with pytest.raises(ValueError, match="Parameter 'nonexistent_parameter' not found in log file"):
        beam_results.get_parameter_stats("nonexistent_parameter")

def test_get_trees(beam_results):
    """Test that get_trees returns a TreeList object."""
    trees = beam_results.get_trees()
    assert trees is not None
    # Add more specific assertions based on the dendropy.TreeList structure

def test_plot_parameters(beam_results, tmp_path):
    """Test that plot_parameters creates an output file."""
    output_file = os.path.join(tmp_path, "test_plot.png")
    beam_results.plot_parameters(output_file=output_file, parameter=beam_results.get_parameters()[10])
    assert os.path.exists(output_file)

def test_plot_parameters_invalid_parameter(beam_results, tmp_path):
    """Test that plot_parameters handles invalid parameter name."""
    output_file = os.path.join(tmp_path, "test_plot.png")
    with pytest.raises(ValueError):
        beam_results.plot_parameters(parameter="nonexistent_parameter", output_file=output_file)

def test_get_consensus_graph(beam_results):
    """Test that get_consensus_graph returns a dictionary with probabilities."""
    graph = beam_results.get_consensus_graph()
    assert isinstance(graph, dict)
    assert all(0 <= prob <= 1 for prob in graph.values())
    assert beam_results.consensus_graph is not None

def test_get_consensus_graph_empty_trees(beam_results):
    """Test that get_consensus_graph handles empty trees list."""
    # Create a copy of beam_results with empty trees
    beam_results.trees = dendropy.TreeList()
    with pytest.raises(ValueError, match="No trees to analyze"):
        beam_results.get_consensus_graph()

def test_get_consensus_graph_invalid_trees(beam_results):
    """Test that get_consensus_graph handles invalid trees."""
    # Create a copy of beam_results with invalid trees
    beam_results.trees = None
    with pytest.raises(ValueError):
        beam_results.get_consensus_graph()

def test_plot_probability_graph(beam_results, tmp_path):
    """Test that plot_probability_graph creates an output file."""
    output_file = os.path.join(tmp_path, "test_plot.png")
    beam_results.plot_probability_graph(output_file=output_file)
    assert os.path.exists(output_file)

def test_plot_thresholded_graph(beam_results, tmp_path):
    """Test that plot_thresholded_graph creates an output file."""
    output_file = os.path.join(tmp_path, "test_plot")
    beam_results.plot_thresholded_graph(output_file_prefix=output_file)
    assert os.path.exists(output_file + "_50.pdf")

def test_plot_thresholded_graph_single_edge(beam_results, tmp_path):
    """Test plot_thresholded_graph with a single edge case."""
    from beam_sup.plotting import plot_thresholded_graph
    
    # Create a simple consensus graph with a single edge that has count 1
    consensus_graph = {"LL_LR_1": 1.0}  # The "_1" suffix means count=1
    
    # Plot thresholded graph
    output_prefix = os.path.join(tmp_path, "test_plot")
    plot_thresholded_graph(
        beam_results.log_data,
        output_prefix,
        primary_tissue="LL",
        threshold=0.5,
        consensus_graph=consensus_graph
    )
    
    # Verify output file was created
    output_file = f"{output_prefix}_50.pdf"
    assert os.path.exists(output_file)

def test_compute_posterior_mutual_info(beam_results, tmp_path):
    """Test that compute_posterior_mutual_info returns expected types."""
    output_file_matrix = os.path.join(tmp_path, "test_plot_matrix.csv")
    output_file_information = os.path.join(tmp_path, "test_plot_information.txt")
    mi, matrix, labels = beam_results.compute_posterior_mutual_info(output_file_matrix=output_file_matrix, output_file_information=output_file_information)
    assert isinstance(mi, float)
    assert isinstance(matrix, np.ndarray)
    assert isinstance(labels, list)
    assert all(isinstance(label, str) for label in labels)
    assert os.path.exists(output_file_matrix)
    assert os.path.exists(output_file_information)

def test_compute_posterior_mutual_info_empty_trees(beam_results):
    """Test that compute_posterior_mutual_info handles empty trees list."""
    # Create a copy of beam_results with empty trees
    beam_results.trees = dendropy.TreeList()
    with pytest.raises(ValueError, match="No trees to analyze"):
        beam_results.compute_posterior_mutual_info()

def test_compute_posterior_mutual_info_invalid_trees(beam_results):
    """Test that compute_posterior_mutual_info handles invalid trees."""
    # Create a copy of beam_results with invalid trees
    beam_results.trees = None
    with pytest.raises(ValueError):
        beam_results.compute_posterior_mutual_info()

def test_sample_and_plot_trees(beam_results, tmp_path):
    """Test that sample_and_plot_trees creates an output file."""
    output_prefix = os.path.join(tmp_path, "test_plot")
    beam_results.sample_and_plot_trees(n=1, output_prefix=output_prefix)
    assert os.path.exists(output_prefix + "_tree_1.pdf")

def test_get_metastasis_times(beam_results, tmp_path):
    """Test that get_metastasis_times returns expected types."""
    output_prefix = os.path.join(tmp_path, "test_plot")
    metastasis_times = beam_results.get_metastasis_times(output_prefix=output_prefix)
    assert isinstance(metastasis_times, dict)
    assert all(isinstance(key, str) for key in metastasis_times)
    assert all(isinstance(value, dict) for value in metastasis_times.values())
    assert os.path.exists(output_prefix + ".pkl")

def test_plot_rate_matrix(beam_results, tmp_path):
    """Test that plot_rate_matrix creates an output file."""
    output_file = os.path.join(tmp_path, "test_plot.png")
    beam_results.plot_rate_matrix(output_file=output_file)
    assert os.path.exists(output_file)

def test_plot_rate_matrix_no_data(beam_results, tmp_path):
    """Test that plot_rate_matrix handles case with no rate data."""
    output_file = os.path.join(tmp_path, "test_plot.png")
    # Create a copy of beam_results with empty log data
    beam_results.log_data = pd.DataFrame()
    with pytest.raises(ValueError):
        beam_results.plot_rate_matrix(output_file=output_file)

def test_load_nonexistent_trees_file():
    """Test that loading a nonexistent trees file raises FileNotFoundError."""
    with pytest.raises(FileNotFoundError):
        data_loader.load_trees_file(NONEXISTENT_FILE)

def test_load_nonexistent_log_file():
    """Test that loading a nonexistent log file raises FileNotFoundError."""
    with pytest.raises(FileNotFoundError):
        data_loader.load_log_file(NONEXISTENT_FILE)

def test_beam_results_invalid_input():
    """Test that BeamResults initialization with invalid input raises appropriate error."""
    with pytest.raises(ValueError):
        BeamResults(
            trees_file=INVALID_TREES_FILE,
            log_file=TEST_LOG_FILE,
            primary_tissue="LL",
            total_time=54.0
        )
    with pytest.raises(ValueError):
        BeamResults(
            trees_file=TEST_TREES_FILE,
            log_file=INVALID_LOG_FILE,
            primary_tissue="LL",
            total_time=54.0
        )

def test_info_method(beam_results, capsys):
    """Test that the info method prints correct information."""
    # First compute consensus graph to ensure it's available
    beam_results.get_consensus_graph()
    
    # Call info method
    beam_results.info()
    
    # Capture the output
    captured = capsys.readouterr()
    output = captured.out
    
    # Verify the output contains expected information
    assert "BeamResults Object Information" in output
    assert f"Trees file: {TEST_TREES_FILE}" in output
    assert f"Log file: {TEST_LOG_FILE}" in output
    assert "Primary tissue: LL" in output
    assert "Total time: 54.0" in output
    assert "Number of trees:" in output
    assert "Number of samples:" in output
    assert "Parameters:" in output
    assert "Migration Analysis:" in output
    assert "Top migrations:" in output

def test_get_consensus_graph_force_recompute(beam_results, tmp_path):
    """Test get_consensus_graph with force_recompute and output_file."""
    # First compute consensus graph
    initial_graph = beam_results.get_consensus_graph()
    
    # Modify the consensus graph
    beam_results.consensus_graph = {"test": 0.5}
    
    # Recompute with force_recompute=True and output_file
    output_file = os.path.join(tmp_path, "consensus_graph.csv")
    recomputed_graph = beam_results.get_consensus_graph(
        force_recompute=True,
        output_file=output_file
    )
    
    # Verify the graph was recomputed (should be different from our modification)
    assert recomputed_graph != {"test": 0.5}
    assert recomputed_graph == initial_graph
    
    # Verify the output file was created
    assert os.path.exists(output_file)
    
    # Verify the output file contents
    with open(output_file, "r") as f:
        content = f.read()
        assert "," in content  # Should contain CSV data

def test_ensure_output_dir(beam_results, tmp_path):
    """Test that _ensure_output_dir creates directories as needed."""
    # Create a path in a non-existent directory
    test_dir = os.path.join(tmp_path, "test_dir")
    test_file = os.path.join(test_dir, "test.txt")
    
    # Call _ensure_output_dir
    beam_results._ensure_output_dir(test_file)
    
    # Verify the directory was created
    assert os.path.exists(test_dir)
    assert os.path.isdir(test_dir)

def test_process_tree_for_consensus(beam_results):
    """Test the _process_tree_for_consensus function with a single tree."""
    from beam_sup.posterior_processing import _process_tree_for_consensus
    
    # Get a single tree from the beam_results
    tree = beam_results.trees[0]
    
    # Process the tree
    counts = _process_tree_for_consensus(tree, "LL")
    
    # Verify the output
    assert isinstance(counts, dict)
    # Verify that counts are non-negative integers
    assert all(isinstance(v, int) and v >= 0 for v in counts.values())
    # Verify that migration keys are in the format "source_target"
    assert all("_" in k for k in counts.keys())

def test_process_tree_wrapper(beam_results):
    """Test the _process_tree_wrapper function."""
    from beam_sup.posterior_processing import _process_tree_wrapper
    
    # Get a single tree from the beam_results
    tree = beam_results.trees[0]
    
    # Process the tree using the wrapper
    counts = _process_tree_wrapper((tree, "LL"))
    
    # Verify the output
    assert isinstance(counts, dict)
    # Verify that counts are non-negative integers
    assert all(isinstance(v, int) and v >= 0 for v in counts.values())

def test_get_consensus_graph_with_burnin(beam_results):
    """Test get_consensus_graph with non-zero burnin."""
    # Set a non-zero burnin percentage
    burnin_percent = 0.1
    
    # Get consensus graph
    graph = beam_results.get_consensus_graph(burnin_percent=burnin_percent)
    
    # Verify the output
    assert isinstance(graph, dict)
    assert all(0 <= prob <= 1 for prob in graph.values())
    assert beam_results.consensus_graph is not None

def test_sample_trees_with_output(beam_results, tmp_path):
    """Test sample_trees with output file generation."""
    from beam_sup.posterior_processing import sample_trees
    
    # Set output prefix
    output_prefix = os.path.join(tmp_path, "test_sampled_trees")
    
    # Sample trees
    sampled_trees = sample_trees(
        beam_results.trees,
        n=2,
        burnin_percent=0.1,
        output_prefix=output_prefix
    )
    
    # Verify the output
    assert isinstance(sampled_trees, list)
    assert len(sampled_trees) == 2
    assert all(isinstance(tree, str) for tree in sampled_trees)
    
    # Verify output file was created
    output_file = f"{output_prefix}_sampled_trees.txt"
    assert os.path.exists(output_file)
    
    # Verify file contents
    with open(output_file, "r") as f:
        content = f.read()
        assert len(content.strip().split("\n")) == 2

def test_get_all_posterior_metastasis_times_with_output(beam_results, tmp_path):
    """Test get_all_posterior_metastasis_times with output file generation."""
    from beam_sup.posterior_processing import get_all_posterior_metastasis_times
    
    # Set output prefix
    output_prefix = os.path.join(tmp_path, "test_metastasis_times")
    
    # Get metastasis times
    met_times = get_all_posterior_metastasis_times(
        beam_results.trees,
        total_time=54.0,
        primary_tissue="LL",
        burnin_percent=0.1,
        output_prefix=output_prefix
    )
    
    # Verify the output
    assert isinstance(met_times, dict)
    assert all(isinstance(key, str) for key in met_times)
    assert all(isinstance(value, dict) for value in met_times.values())
    
    # Verify output file was created
    output_file = f"{output_prefix}.pkl"
    assert os.path.exists(output_file)

def test_get_consensus_graph_invalid_input(beam_results):
    """Test get_consensus_graph with invalid input types."""
    from beam_sup.posterior_processing import get_consensus_graph
    
    # Test with non-TreeList input
    with pytest.raises(ValueError, match="Trees must be a dendropy.TreeList object"):
        get_consensus_graph(None, "LL")
    
    # Test with empty TreeList
    empty_trees = dendropy.TreeList()
    with pytest.raises(ValueError, match="No trees to analyze"):
        get_consensus_graph(empty_trees, "LL")

def test_sample_trees_node_labeling(beam_results):
    """Test sample_trees with nodes that have no labels."""
    from beam_sup.posterior_processing import sample_trees
    
    # Create a copy of the first tree
    tree = beam_results.trees[0].clone()
    
    # Remove labels from some nodes
    for node in tree.preorder_node_iter():
        if node.taxon:
            node.taxon.label = None
        else:
            node.label = None
    
    # Create a new TreeList with the modified tree
    trees = dendropy.TreeList()
    trees.append(tree)
    
    # Sample trees
    sampled_trees = sample_trees(trees, n=1)
    
    # Verify that all nodes got labels
    assert len(sampled_trees) == 1
    tree_str = sampled_trees[0]
    assert "node" in tree_str  # Should contain generated node labels

def test_get_all_posterior_metastasis_times_non_ultrametric(beam_results):
    """Test get_all_posterior_metastasis_times with non-ultrametric tree."""
    from beam_sup.posterior_processing import get_all_posterior_metastasis_times
    
    # Create a copy of the first tree
    tree = beam_results.trees[0].clone()
    
    # Modify edge lengths to make tree non-ultrametric
    for node in tree.preorder_node_iter():
        if node.edge_length is not None:
            node.edge_length *= random.uniform(0.5, 1.5)
    
    # Create a new TreeList with the modified tree
    trees = dendropy.TreeList()
    trees.append(tree)
    
    # Try to get metastasis times
    with pytest.raises(ValueError, match="Tree is not ultrametric"):
        get_all_posterior_metastasis_times(trees, total_time=54.0, primary_tissue="LL")



