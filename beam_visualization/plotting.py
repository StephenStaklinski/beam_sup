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

def plot_sampled_tree(
    newick_str: str,
    primary_tissue: str,
    total_time: float,
    output_prefix: str,
    tree_num: int
) -> None:
    """
    Plot a sampled tree with its migration graph and timing.
    
    Args:
        newick_str: Newick string of the tree with location annotations
        primary_tissue: Primary tissue label
        total_time: Total time of the tree
        output_prefix: Prefix for output files
        tree_num: Number of this tree in the sample
    """
    import os
    from ete3 import Tree, TreeStyle, NodeStyle, CircleFace, TextFace
    import re
    
    # Parse tree and get tissues
    tree = Tree(newick_str, format=1)
    all_tissues = {primary_tissue}
    
    # Process nodes and get tissues
    i = 1
    for node in tree.traverse():
        # Extract location from annotation if present
        if '&location=' in node.name:
            match = re.search(r'&location="([^"]+)"', node.name)
            if match:
                node.tissue = match.group(1)
                # Keep original name without annotation
                node.name = node.name.split('[')[0]
        else:
            # For nodes without location annotation, set tissue to none
            node.tissue = None
                
        all_tissues.add(node.tissue)
    
    # Check ultrametricity
    dists = {round(node.get_distance(tree), 5) for node in tree.traverse() if node.is_leaf()}
    if len(dists) != 1:
        print(f"WARNING: Tree {tree_num} is not ultrametric")
    
    # Add origin node
    tree_height = dists.pop()
    origin = Tree(name="origin", dist=0)
    origin.tissue = primary_tissue
    root = tree.get_tree_root()
    root.dist = total_time - tree_height
    origin.add_child(root)
    
    # Get migration counts
    migration_counts = {}
    for node in origin.traverse():
        if node.is_root():
            continue
        if node.up.tissue != node.tissue:
            migration = f"{node.up.tissue}_{node.tissue}"
            migration_counts[migration] = migration_counts.get(migration, 0) + 1
    
    # Get colors
    all_tissues = sorted(list(all_tissues - {primary_tissue}))
    custom_colors = {
        node: color
        for node, color in zip(all_tissues, DEFAULT_COLORS[0:len(all_tissues)])
        if node != primary_tissue
    }
    all_tissues = [primary_tissue] + all_tissues
    custom_colors[primary_tissue] = "black"
    custom_colors[None] = "white"
    
    # Plot migration graph
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
    
    for edge, count in migration_counts.items():
        source, target = edge.split("_")
        G.add_edge(
            source,
            target,
            color=f'"{custom_colors[source]};0.5:{custom_colors[target]}"',
            penwidth=3,
            label="" if count == 1 else str(count),
            fontsize=24,
        )
    
    dot = nx.nx_pydot.to_pydot(G)
    dot.write_pdf(f"{output_prefix}_migration_graph_{tree_num}.pdf")
    
    # Plot migration timing
    metastasis_times = {}
    for node in origin.traverse():
        if node.is_root():
            continue
        if node.up.tissue != node.tissue:
            migration = f"{node.up.tissue}_{node.tissue}"
            time = node.up.get_distance(origin) + (node.dist / 2)
            metastasis_times.setdefault(migration, []).append(time)
    
    if metastasis_times:  # Only plot if we have migration events
        fig, ax = plt.subplots(figsize=(12, 3))
        for migration, times in metastasis_times.items():
            source, target = migration.split("_")
            for time in times:
                ax.plot([time, time], [0.5, 1], color=custom_colors[source], linewidth=3)
                ax.plot([time, time], [0, 0.5], color=custom_colors[target], linewidth=3)
        
        ax.set_xlim(0, total_time)
        ax.set_ylim(0, 1)
        ax.set_xlabel("Time", fontsize=18)
        ax.set_yticks([0.25, 0.75])
        ax.set_yticklabels(["Target", "Source"], fontsize=18)
        ax.tick_params(axis="x", labelsize=18)
        
        handles = [plt.Line2D([0], [0], color=color, lw=3) for color in custom_colors.values()]
        ax.legend(
            handles,
            custom_colors.keys(),
            bbox_to_anchor=(1.05, 1),
            loc="upper left",
            borderaxespad=0.0,
            frameon=False,
            fontsize=14,
        )
        
        plt.tight_layout()
        plt.savefig(f"{output_prefix}_migration_timing_{tree_num}.pdf")
        plt.close()

    # Plot tree
    os.environ['QT_QPA_PLATFORM'] = 'offscreen'
    ts = TreeStyle()
    ts.rotation = 90
    ts.scale = 1
    ts.show_leaf_name = True
    ts.show_branch_length = False
    ts.show_border = False
    ts.show_scale = False
    ts.mode = "r"
    
    for tissue in all_tissues:
        ts.legend.add_face(CircleFace(10, custom_colors[tissue]), column=0)
        ts.legend.add_face(TextFace(tissue, fsize=12), column=1)
    
    for node in origin.traverse():
        nstyle = NodeStyle()
        nstyle["shape"] = "circle"
        nstyle["size"] = 10
        nstyle["hz_line_color"] = custom_colors[node.tissue]
        nstyle["vt_line_color"] = custom_colors[node.tissue]
        nstyle["hz_line_width"] = 3
        nstyle["vt_line_width"] = 3
        nstyle["fgcolor"] = custom_colors[node.tissue]
        node.set_style(nstyle)
    
    origin.render(f"{output_prefix}_tree_{tree_num}.pdf", tree_style=ts)

