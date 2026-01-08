import os
import pandas as pd
from dendropy import TreeList
from typing import Tuple


def load_trees_file(trees_file: str) -> TreeList:
    """
    Load and parse a BEAM trees file.

    Args:
        trees_file: Path to the .trees file

    Returns:
        TreeList: The loaded trees
    """
    if not os.path.exists(trees_file):
        raise FileNotFoundError(f"Trees file not found: {trees_file}")

    try:
        # Try loading as NEXUS format
        trees = TreeList.get(path=trees_file, schema="nexus")
    except Exception as e:
        raise ValueError(f"Failed to load trees file: {e}")

    return trees


def load_log_file(log_file: str) -> pd.DataFrame:
    """
    Load and parse a BEAM log file.

    Args:
        log_file: Path to the .log file

    Returns:
        pd.DataFrame: The loaded log data
    """
    if not os.path.exists(log_file):
        raise FileNotFoundError(f"Log file not found: {log_file}")

    try:
        log_data = pd.read_csv(log_file, sep="\t", comment="#")

        required_columns = ["Sample", "posterior", "likelihood", "prior"]
        if not all(col in log_data.columns for col in required_columns):
            raise ValueError("Log file missing required columns")
    except Exception as e:
        raise ValueError(f"Failed to load log file: {e}")

    return log_data


def load_beam_files(
    trees_file: str, log_file: str
) -> Tuple[TreeList, pd.DataFrame]:
    """
    Load both BEAM output files.

    Args:
        trees_file: Path to the .trees file
        log_file: Path to the .log file

    Returns:
        Tuple containing:
            - dendropy.TreeList: The loaded trees
            - pd.DataFrame: The loaded log data
    """
    trees = load_trees_file(trees_file)
    
    log_data = None
    if log_file is not None:
        log_data = load_log_file(log_file)
        
    return trees, log_data
