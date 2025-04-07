"""BEAM visualization package."""

from .main import BeamResults
from .plotting import plot_parameters
from .formatting import get_consensus_graph
from .data_loader import load_beam_files

__version__ = "0.1.0"
__all__ = [
    'BeamResults',
    'plot_parameters',
    'get_consensus_graph',
    'load_beam_files'
] 