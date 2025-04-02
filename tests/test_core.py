import os
import pytest
import pandas as pd
import dendropy
from beam_visualization import BeamResults

@pytest.fixture
def sample_trees_file(tmp_path):
    """Create a sample trees file."""
    trees_file = tmp_path / "sample.trees"
    with open(trees_file, "w") as f:
        f.write("""#NEXUS
BEGIN TREES;
[!R] TREE tree_1 = ((A:1,B:1):1,C:2);
[!R] TREE tree_2 = ((A:1,C:1):1,B:2);
END;
""")
    return str(trees_file)

@pytest.fixture
def sample_log_file(tmp_path):
    """Create a sample log file."""
    log_file = tmp_path / "sample.log"
    df = pd.DataFrame({
        "Sample": range(100),
        "posterior": [0.1] * 100,
        "likelihood": [0.2] * 100,
        "prior": [0.3] * 100,
        "rate": [0.5] * 100
    })
    df.to_csv(log_file, sep="\t", index=False)
    return str(log_file)

def test_beam_results_initialization(sample_trees_file, sample_log_file):
    """Test that BeamResults can be initialized with valid files."""
    results = BeamResults(sample_trees_file, sample_log_file)
    assert results.trees_file == sample_trees_file
    assert results.log_file == sample_log_file
    assert results.trees is not None
    assert results.log_data is not None

def test_beam_results_invalid_files():
    """Test that BeamResults raises appropriate errors for invalid files."""
    with pytest.raises(FileNotFoundError):
        BeamResults("nonexistent.trees", "nonexistent.log")

def test_plot_parameters(sample_trees_file, sample_log_file):
    """Test the plot_parameters method."""
    results = BeamResults(sample_trees_file, sample_log_file)
    
    # Test plotting all parameters
    results.plot_parameters()
    
    # Test plotting specific parameter
    results.plot_parameters("rate")
    
    # Test invalid parameter
    with pytest.raises(ValueError):
        results.plot_parameters("nonexistent_parameter")
