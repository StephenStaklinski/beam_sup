import pytest
import pandas as pd
import os

@pytest.fixture(scope="session")
def test_data_dir():
    """Get the path to the test data directory."""
    return os.path.join(os.path.dirname(__file__), "data")

@pytest.fixture(scope="session")
def ensure_test_data_dir(test_data_dir):
    """Ensure the test data directory exists."""
    os.makedirs(test_data_dir, exist_ok=True)
    return test_data_dir 