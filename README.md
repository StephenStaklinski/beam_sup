# BEAM Visualization

A Python package for visualizing and analyzing outputs from the BEAM phylogenetic tree analysis package.

## Installation

You can install the package using pip:

```bash
pip install beam_visualization
```

For development installation (including test dependencies):

```bash
git clone https://github.com/StephenStaklinski/beam_visualization.git
cd beam_visualization
pip install -e ".[test]"
```

Or similarly, but with conda:

```bash
git clone https://github.com/StephenStaklinski/beam_visualization.git
conda create -n beam_visualization python
conda activate beam_visualization
pip install -e ".[test]"

```

## Usage

The package provides functionality to load and visualize BEAM outputs, including:
- Parameter distributions from log files
- Migration graphs
- Rate matrices

### Basic Usage

```python
from beam_visualization import BeamResults

data = "./examples/data/"
trees = data + "example.trees" 
log = data + "example.log"

# Load BEAM results
results = BeamResults(trees, log)
results.info()

# Plot parameter distributions
results.plot_parameters(output_file="./all_parameters.pdf")  # Plot all parameters

param = results.get_parameters()[14]    # View the names of parameters that were logged and pick one
results.plot_parameters(param, output_file="./rate.pdf")  # Plot specific parameter
results.get_parameter_stats(param)  # Calculate basic statistics about specific parameter
```

## Features

- Load and parse BEAM output files (.trees and .log)
- Plot parameter distributions and statistics
- Generate migration graphs
- Analyze rate matrices

## Development

### Running Tests

To run the test suite:

```bash
# Install development dependencies
pip install -e ".[test]"

# Run tests with coverage report
pytest

# Run tests without coverage report
pytest --no-cov
```

The test suite includes:
- Unit tests for all major functionality
- Test fixtures for sample data
- Coverage reporting
- Error handling tests

## License

This project is licensed under the MIT License - see the LICENSE file for details.
