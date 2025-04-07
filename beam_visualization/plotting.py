"""Plotting functions for BEAM visualization."""

import matplotlib.pyplot as plt
import seaborn as sns
from typing import Optional, Dict, Union, List
import pandas as pd
import networkx as nx
from .config import DEFAULT_COLORS

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

def plot_probability_graph(
    data: pd.DataFrame,
    primary_tissue: Optional[str] = None,
    output_file: Optional[str] = None,
    consensus_graph: Optional[Dict[str, float]] = None
) -> None:
    """
    Plot the consensus migration graph with edge thicknesses proportional to probability.
    
    Args:
        data: DataFrame containing the log data
        primary_tissue: Primary tissue label for migration analysis
        output_file: Optional path to save the plot. If None, displays the plot.
        consensus_graph: Pre-computed consensus graph
    """
    # Find all tissues to set the node colors
    all_tissues = sorted(
        list(
            set([value for node in consensus_graph.keys() for value in node.split("_")[0:2]])
            - {primary_tissue}
        )
    )
    
    # Assign colors to tissues
    custom_colors = {
        node: color
        for node, color in zip(all_tissues, DEFAULT_COLORS[0:len(all_tissues)])
        if node != primary_tissue
    }
    custom_colors[primary_tissue] = "black"
    all_tissues = [primary_tissue] + all_tissues
    
    # Plot probability graph with edge thicknesses proportional to probability
    G = nx.MultiDiGraph()
    for node in all_tissues:
        G.add_node(
            node,
            color=custom_colors[node],
            shape="box",
            fillcolor="white",
            penwidth=3.0,
            fontsize=32,
        )
    
    for edge, probability in consensus_graph.items():
        source, target, num = edge.split("_")
        G.add_edge(
            source,
            target,
            color=f'"{custom_colors[source]};0.5:{custom_colors[target]}"',
            penwidth=probability * 3,
            fontsize=24,
        )
    
    dot = nx.nx_pydot.to_pydot(G)
    if output_file:
        dot.write_pdf(output_file)
    else:
        dot.write_pdf("temp_plot.pdf")
        plt.figure(figsize=(12, 8))
        plt.imshow(plt.imread("temp_plot.pdf"))
        plt.axis('off')
        plt.tight_layout()
        plt.show()
        import os
        os.remove("temp_plot.pdf")

def plot_thresholded_graph(
    data: pd.DataFrame,
    primary_tissue: Optional[str] = None,
    threshold: Union[float, List[float]] = 0.5,
    output_file_prefix: Optional[str] = None,
    consensus_graph: Optional[Dict[str, float]] = None
) -> None:
    """
    Plot the thresholded consensus migration graph with collapsed multiedges.
    
    Args:
        data: DataFrame containing the log data
        primary_tissue: Primary tissue label for migration analysis
        threshold: Single threshold value or list of thresholds for filtering edges
        output_file_prefix: Optional prefix for output files. If None, displays the plot.
                           For multiple thresholds, each file will be named {prefix}_{threshold}.pdf
        consensus_graph: Pre-computed consensus graph
    """
    # Convert single threshold to list for uniform handling
    thresholds = [threshold] if isinstance(threshold, (int, float)) else threshold
    
    # Find all tissues to set the node colors
    all_tissues = sorted(
        list(
            set([value for node in consensus_graph.keys() for value in node.split("_")[0:2]])
            - {primary_tissue}
        )
    )
    
    # Assign colors to tissues
    custom_colors = {
        node: color
        for node, color in zip(all_tissues, DEFAULT_COLORS[0:len(all_tissues)])
        if node != primary_tissue
    }
    custom_colors[primary_tissue] = "black"
    all_tissues = [primary_tissue] + all_tissues
    
    for current_threshold in thresholds:
        # Plot thresholded graph with collapsed multiedges
        G = nx.MultiDiGraph()
        for node in all_tissues:
            G.add_node(
                node,
                color=custom_colors[node],
                shape="box",
                fillcolor="white",
                penwidth=3.0,
                fontsize=32,
            )
        
        for edge, probability in consensus_graph.items():
            if probability > current_threshold:
                source, target, num = edge.split("_")
                if G.has_edge(source, target):
                    G[source][target][0]["label"] = str(int(G[source][target][0]["label"]) + 1)
                else:
                    G.add_edge(
                        source,
                        target,
                        color=f'"{custom_colors[source]};0.5:{custom_colors[target]}"',
                        penwidth=3,
                        label="1",
                        fontsize=24,
                    )
        
        # Set empty label for single edges
        for source, target, data in G.edges(data=True):
            if data.get("label") == "1":
                data["label"] = ""
        
        dot = nx.nx_pydot.to_pydot(G)
        if output_file_prefix:
            out_threshold = str(int(current_threshold * 100))
            output_file = f"{output_file_prefix}_{out_threshold}.pdf"
            dot.write_pdf(output_file)
        else:
            dot.write_pdf("temp_plot.pdf")
            plt.figure(figsize=(12, 8))
            plt.imshow(plt.imread("temp_plot.pdf"))
            plt.axis('off')
            plt.tight_layout()
            plt.show()
            import os
            os.remove("temp_plot.pdf")

