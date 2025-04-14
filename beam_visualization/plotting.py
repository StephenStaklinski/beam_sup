"""Plotting functions for BEAM visualization."""

import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import seaborn as sns
from typing import Optional, Dict, List, Union, Tuple
import pandas as pd
import networkx as nx
import numpy as np
from collections import defaultdict
from ete3 import Tree, TreeStyle, NodeStyle, TextFace, CircleFace
import os
import re
from .config import (
    DEFAULT_COLORS,
    DEFAULT_NODE_STYLES,
    DEFAULT_PLOT_STYLES,
    DEFAULT_TREE_STYLE,
    DEFAULT_FIGURE_SIZE,
    DEFAULT_FONT_SIZE,
    DEFAULT_MIN_PROB_THRESHOLD,
)


def plot_parameters(
    data: pd.DataFrame,
    parameter: Optional[str] = None,
    output_file: Optional[str] = None,
) -> None:
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

    plt.figure(figsize=DEFAULT_FIGURE_SIZE)
    sns.histplot(data=data, kde=False)
    plt.title(
        f"Distribution of {parameter if parameter else 'Parameters'}",
        fontsize=DEFAULT_FONT_SIZE,
    )
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
    consensus_graph: Optional[Dict[str, float]] = None,
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
            set(
                [
                    value
                    for node in consensus_graph.keys()
                    for value in node.split("_")[0:2]
                ]
            )
            - {primary_tissue}
        )
    )

    # Assign colors to tissues
    custom_colors = {
        node: color
        for node, color in zip(all_tissues, DEFAULT_COLORS[0 : len(all_tissues)])
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
            shape=DEFAULT_NODE_STYLES["shape"],
            fillcolor="white",
            penwidth=DEFAULT_NODE_STYLES["line_width"],
            fontsize=DEFAULT_FONT_SIZE,
        )

    for edge, probability in consensus_graph.items():
        source, target, num = edge.split("_")
        G.add_edge(
            source,
            target,
            color=f'"{custom_colors[source]};0.5:{custom_colors[target]}"',
            penwidth=probability * DEFAULT_NODE_STYLES["line_width"],
            fontsize=DEFAULT_FONT_SIZE - 8,
        )

    dot = nx.nx_pydot.to_pydot(G)
    if output_file:
        dot.write_pdf(output_file)
    else:
        dot.write_pdf("temp_plot.pdf")
        plt.figure(figsize=DEFAULT_FIGURE_SIZE)
        plt.imshow(plt.imread("temp_plot.pdf"))
        plt.axis("off")
        plt.tight_layout()
        plt.show()
        os.remove("temp_plot.pdf")


def plot_thresholded_graph(
    data: pd.DataFrame,
    primary_tissue: Optional[str] = None,
    threshold: Union[float, List[float]] = DEFAULT_MIN_PROB_THRESHOLD,
    output_file_prefix: Optional[str] = None,
    consensus_graph: Optional[Dict[str, float]] = None,
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
            set(
                [
                    value
                    for node in consensus_graph.keys()
                    for value in node.split("_")[0:2]
                ]
            )
            - {primary_tissue}
        )
    )

    # Assign colors to tissues
    custom_colors = {
        node: color
        for node, color in zip(all_tissues, DEFAULT_COLORS[0 : len(all_tissues)])
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
                shape=DEFAULT_NODE_STYLES["shape"],
                fillcolor="white",
                penwidth=DEFAULT_NODE_STYLES["line_width"],
                fontsize=DEFAULT_FONT_SIZE,
            )

        for edge, probability in consensus_graph.items():
            if probability > current_threshold:
                source, target, num = edge.split("_")
                if G.has_edge(source, target):
                    G[source][target][0]["label"] = str(
                        int(G[source][target][0]["label"]) + 1
                    )
                else:
                    G.add_edge(
                        source,
                        target,
                        color=f'"{custom_colors[source]};0.5:{custom_colors[target]}"',
                        penwidth=DEFAULT_NODE_STYLES["line_width"],
                        label="1",
                        fontsize=DEFAULT_FONT_SIZE - 8,
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
            plt.figure(figsize=DEFAULT_FIGURE_SIZE)
            plt.imshow(plt.imread("temp_plot.pdf"))
            plt.axis("off")
            plt.tight_layout()
            plt.show()
            os.remove("temp_plot.pdf")


def plot_sampled_tree(
    newick_str: str,
    primary_tissue: str,
    total_time: float,
    output_prefix: str,
    tree_num: int,
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

    # Parse tree and get tissues
    tree = Tree(newick_str, format=1)
    all_tissues = {primary_tissue}

    # Process nodes and get tissues
    i = 1
    for node in tree.traverse():
        # Extract location from annotation if present
        if "&location=" in node.name:
            match = re.search(r'&location="([^"]+)"', node.name)
            if match:
                node.tissue = match.group(1)
                # Keep original name without annotation
                node.name = node.name.split("[")[0]
        else:
            # For nodes without location annotation, set tissue to none
            node.tissue = None

        all_tissues.add(node.tissue)

    # Check ultrametricity
    dists = {
        round(node.get_distance(tree), 5) for node in tree.traverse() if node.is_leaf()
    }
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
        for node, color in zip(all_tissues, DEFAULT_COLORS[0 : len(all_tissues)])
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
            shape=DEFAULT_NODE_STYLES["shape"],
            fillcolor="white",
            penwidth=DEFAULT_NODE_STYLES["line_width"],
            fontsize=DEFAULT_FONT_SIZE,
        )

    for edge, count in migration_counts.items():
        source, target = edge.split("_")
        G.add_edge(
            source,
            target,
            color=f'"{custom_colors[source]};0.5:{custom_colors[target]}"',
            penwidth=DEFAULT_NODE_STYLES["line_width"],
            label="" if count == 1 else str(count),
            fontsize=DEFAULT_PLOT_STYLES["legend_fontsize"] - 8,
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
        fig, ax = plt.subplots(figsize=DEFAULT_FIGURE_SIZE)
        for migration, times in metastasis_times.items():
            source, target = migration.split("_")
            for time in times:
                ax.plot(
                    [time, time],
                    [0.5, 1],
                    color=custom_colors[source],
                    linewidth=DEFAULT_PLOT_STYLES["linewidth"],
                )
                ax.plot(
                    [time, time],
                    [0, 0.5],
                    color=custom_colors[target],
                    linewidth=DEFAULT_PLOT_STYLES["linewidth"],
                )

        ax.set_xlim(0, total_time)
        ax.set_ylim(0, 1)
        ax.set_xlabel("Time", fontsize=DEFAULT_FONT_SIZE)
        ax.set_yticks([0.25, 0.75])
        ax.set_yticklabels(["Target", "Source"], fontsize=DEFAULT_FONT_SIZE)
        ax.tick_params(axis="x", labelsize=DEFAULT_FONT_SIZE)

        handles = [
            plt.Line2D([0], [0], color=color, lw=DEFAULT_PLOT_STYLES["linewidth"])
            for color in custom_colors.values()
        ]
        ax.legend(
            handles,
            custom_colors.keys(),
            bbox_to_anchor=DEFAULT_PLOT_STYLES["legend_position"],
            loc="upper left",
            borderaxespad=0.0,
            frameon=False,
            fontsize=DEFAULT_PLOT_STYLES["legend_fontsize"] - 8,
        )

        plt.tight_layout()
        plt.savefig(f"{output_prefix}_migration_timing_{tree_num}.pdf")
        plt.close()

    # Plot tree
    os.environ["QT_QPA_PLATFORM"] = "offscreen"
    ts = TreeStyle()
    ts.rotation = DEFAULT_TREE_STYLE["rotation"]
    ts.scale = DEFAULT_TREE_STYLE["scale"]
    ts.show_leaf_name = DEFAULT_TREE_STYLE["show_leaf_name"]
    ts.show_branch_length = DEFAULT_TREE_STYLE["show_branch_length"]
    ts.show_border = DEFAULT_TREE_STYLE["show_border"]
    ts.show_scale = DEFAULT_TREE_STYLE["show_scale"]
    ts.mode = DEFAULT_TREE_STYLE["mode"]

    # Add legend
    for tissue in all_tissues:
        if tissue is not None:  # Skip None tissue type
            ts.legend.add_face(
                CircleFace(DEFAULT_NODE_STYLES["size"], custom_colors[tissue]), column=0
            )
            ts.legend.add_face(TextFace(tissue, fsize=DEFAULT_FONT_SIZE), column=1)

    # Set node styles
    for node in origin.traverse():
        nstyle = NodeStyle()
        nstyle["shape"] = "circle"
        nstyle["size"] = DEFAULT_NODE_STYLES["size"]
        nstyle["hz_line_color"] = custom_colors[node.tissue]
        nstyle["vt_line_color"] = custom_colors[node.tissue]
        nstyle["hz_line_width"] = DEFAULT_NODE_STYLES["line_width"]
        nstyle["vt_line_width"] = DEFAULT_NODE_STYLES["line_width"]
        nstyle["fgcolor"] = custom_colors[node.tissue]
        node.set_style(nstyle)

    origin.render(f"{output_prefix}_tree_{tree_num}.pdf", tree_style=ts)


def plot_metastasis_timing(
    met_times: Dict[str, Dict[str, Tuple[float, float]]],
    consensus_graph: Dict[str, float],
    origin_time: float,
    origin_tissue: str,
    min_prob_threshold: float = DEFAULT_MIN_PROB_THRESHOLD,
    output_prefix: Optional[str] = None,
) -> None:
    """
    Plot metastasis timing distributions.

    Args:
        met_times: Dictionary mapping tree labels to dictionaries of metastasis events
        consensus_graph: Dictionary mapping migration patterns to their probabilities
        origin_time: Time of origin
        origin_tissue: Primary tissue label
        min_prob_threshold: Minimum probability threshold for migrations to include in plots
        output_prefix: Optional prefix for output files
    """
    
    allowable_migrations = set()
    for migration, prob in consensus_graph.items():
        if float(prob) >= min_prob_threshold:
            allowable_migrations.add(migration)
    migration_counts = defaultdict(lambda: np.zeros(origin_time + 1))
    migration_counts_mid_points = defaultdict(list)

    for graph in met_times.values():
        for migration, time in graph.items():
            if migration not in allowable_migrations:
                continue
            # range is from origin at 0 to the end of experiment at origin_time
            start_range = round(time[0])
            end_range = round(time[1])
            num_intervals = end_range - start_range + 1
            prob = 1 / (len(met_times) * num_intervals)
            migration_counts[migration][start_range : end_range + 1] += prob
            migration_counts_mid_points[migration].append((start_range + end_range) / 2)

    df = pd.DataFrame(migration_counts, index=np.arange(0, origin_time + 1)).T
    
    # Split the migration strings into source and target tissues
    df.index = pd.MultiIndex.from_tuples(
        [
            tuple([migration.split("_")[0], "_".join(migration.split("_")[1:])])
            for migration in df.index
        ],
        names=["source", "target"],
    )

    remaining_sources = sorted(
        set(df.index.get_level_values("source").unique()) - {origin_tissue}
    )
    remaining_sources.insert(0, origin_tissue)

    # Get unique tissues
    target_tissues = sorted(set(df.index.get_level_values("target")))
    target_tissues_no_num = set([tissue.split("_")[0] for tissue in target_tissues])
    tissues = sorted(
        set(df.index.get_level_values("source")).union(target_tissues_no_num)
        - {origin_tissue}
    )
    tissues.insert(0, origin_tissue)

    # Create a color palette for the tissues
    all_tissues = sorted(list(set(tissues) - {origin_tissue}))
    custom_colors = DEFAULT_COLORS
    custom_colors = {
        node: color for node, color in zip(all_tissues, custom_colors[0 : len(all_tissues)])
    }
    custom_colors[origin_tissue] = "black"

    # Create a grid of subplots with one row per source tissue
    fig, axes = plt.subplots(
        len(remaining_sources),
        1,
        figsize=(15, 3 * len(remaining_sources)),
        sharex=True,
        sharey=True,
    )

    for i, source in enumerate(remaining_sources):
        if len(remaining_sources) == 1:
            ax = axes
        else:
            ax = axes[i]
        y = 0.1
        for target in target_tissues:
            if (source, target) in df.index:
                y += 1
                if "_1" in target:
                    target_reformatted = target.split("_")[0]
                else:
                    target_reformatted = target
                target_name = target.split("_")[0]
                sns.lineplot(
                    x=df.columns,
                    y=df.loc[(source, target)],
                    ax=ax,
                    color=custom_colors[target_name],
                )
                ax.fill_between(
                    df.columns,
                    df.loc[(source, target)],
                    alpha=0.3,
                    color=custom_colors[target_name],
                )

        ax.set_ylabel(source, fontsize=DEFAULT_FONT_SIZE)
        ax.tick_params(axis="both", which="major", labelsize=DEFAULT_FONT_SIZE)
        # ax.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)
        if i == len(remaining_sources) - 1:
            ax.set_xlabel("Time", fontsize=DEFAULT_FONT_SIZE)

    # Create a single legend for all axes
    handles = [
        plt.Line2D([0], [0], color=color, lw=2) for tissue, color in custom_colors.items()
    ]
    labels = list(custom_colors.keys())
    fig.legend(
        handles,
        labels,
        bbox_to_anchor=(1.05, 0.75),
        title="Target Tissue",
        frameon=False,
        fontsize=DEFAULT_FONT_SIZE,
        title_fontsize=DEFAULT_FONT_SIZE,
    )

    fig.text(
        0.001,
        0.5,
        "Source tissue",
        va="center",
        ha="center",
        rotation="vertical",
        fontsize=DEFAULT_FONT_SIZE,
    )

    plt.tight_layout(rect=[0.02, 0, 0.88, 1])
    if output_prefix:
        plt.savefig(f"{output_prefix}_metastasis_timing_prob.pdf", bbox_inches='tight')
    else:
        plt.show()
    plt.close()


    # plot midpoints
    fig, axes = plt.subplots(
        len(remaining_sources),
        1,
        figsize=(15, 3 * len(remaining_sources)),
        sharex=True,
        sharey=True,
    )
    for i, source in enumerate(remaining_sources):
        if len(remaining_sources) == 1:
            ax = axes
        else:
            ax = axes[i]
        for target in target_tissues:
            migration = f"{source}_{target}"
            if migration in migration_counts_mid_points:
                if "_1" in target:
                    target_reformatted = target.split("_")[0]
                else:
                    target_reformatted = target
                target_name = target.split("_")[0]
                ax.hist(
                    migration_counts_mid_points[migration],
                    bins=100,
                    color=custom_colors[target_name],
                    alpha=0.6,
                    label=target_reformatted,
                )
        ax.set_ylabel(source, fontsize=DEFAULT_FONT_SIZE)
        ax.tick_params(axis="both", which="major", labelsize=DEFAULT_FONT_SIZE)
        if i == len(remaining_sources) - 1:
            ax.set_xlabel("Time", fontsize=DEFAULT_FONT_SIZE)

    fig.legend(
        handles,
        labels,
        bbox_to_anchor=(1.05, 0.75),
        title="Target Tissue",
        frameon=False,
        fontsize=DEFAULT_FONT_SIZE,
        title_fontsize=DEFAULT_FONT_SIZE,
    )
    fig.text(
        0.001,
        0.5,
        "Source tissue",
        va="center",
        ha="center",
        rotation="vertical",
        fontsize=DEFAULT_FONT_SIZE,
    )

    plt.tight_layout()

    if output_prefix:
        plt.savefig(f"{output_prefix}_metastasis_timing_midpoints.pdf", bbox_inches='tight')
    else:
        plt.show()
    plt.close()
