import matplotlib.pyplot as plt
import seaborn as sns
from typing import Optional
import pandas as pd

def plot_parameters(data: pd.DataFrame, parameter: Optional[str] = None, output_file: Optional[str] = None) -> None:
    """
    Plot parameter distributions from the log file.
    
    Args:
        data (pd.DataFrame): The log data
        parameter (str, optional): Specific parameter to plot. If None, plots all parameters.
        output_file (str, optional): Path to save the plot. If None, displays the plot.
    """
    if parameter:
        if parameter not in data.columns:
            raise ValueError(f"Parameter '{parameter}' not found in log file")
        data = data[parameter]
    
    plt.figure(figsize=(12, 6))
    sns.histplot(data=data, kde=False)
    plt.title(f"Distribution of {parameter if parameter else 'Parameters'}")
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file)
        plt.close()
    else:
        plt.show() 