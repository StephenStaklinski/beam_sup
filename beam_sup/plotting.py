"""Plotting functions"""

import matplotlib.pyplot as plt
import matplotlib

matplotlib.use("Agg")
import seaborn as sns
from typing import Optional, Dict, List, Union, Tuple
import pandas as pd
import networkx as nx
import numpy as np
from collections import defaultdict
from ete3 import Tree, TreeStyle, NodeStyle, TextFace, CircleFace
import os
import re

from . import config


def plot_parameters(
    data: pd.DataFrame,
    output_file: str,
    parameter: Optional[str] = None,
) -> None:
    """
    Plot parameter distributions from the log file.

    Args:
        data (pd.DataFrame): The log data
        output_file (str): Path to save the plot.
        parameter (str, optional): Specific parameter to plot. If None, plots all parameters.
    """
    if parameter:
        if parameter not in data.columns:
            raise ValueError(f"Parameter '{parameter}' not found in log file")
        data = data[parameter]

    plt.figure(figsize=config.DEFAULT_FIGURE_SIZE)
    sns.histplot(data=data, kde=False)
    plt.title(
        f"Distribution of {parameter if parameter else 'Parameters'}",
        fontsize=config.DEFAULT_FONT_SIZE,
    )
    plt.tight_layout()

    plt.savefig(output_file)
    plt.close()


def plot_probability_graph(
    data: pd.DataFrame,
    output_file: str,
    primary_tissue: Optional[str] = None,
    consensus_graph: Optional[Dict[str, float]] = None,
) -> None:
    """
    Plot the consensus migration graph with edge thicknesses proportional to probability.

    Args:
        data: DataFrame containing the log data
        output_file: Path to save the plot
        primary_tissue: Primary tissue label for migration analysis
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
        for node, color in zip(all_tissues, config.DEFAULT_COLORS[0 : len(all_tissues)])
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
            shape=config.DEFAULT_NODE_STYLES["shape"],
            fillcolor="white",
            penwidth=config.DEFAULT_NODE_STYLES["line_width"],
            fontsize=config.DEFAULT_FONT_SIZE,
        )

    for edge, probability in consensus_graph.items():
        source, target, num = edge.split("_")
        G.add_edge(
            source,
            target,
            color=f'"{custom_colors[source]};0.5:{custom_colors[target]}"',
            penwidth=probability * config.DEFAULT_NODE_STYLES["line_width"],
            fontsize=config.DEFAULT_FONT_SIZE - 8,
        )

    dot = nx.nx_pydot.to_pydot(G)
    dot.write_pdf(output_file)


def plot_thresholded_graph(
    data: pd.DataFrame,
    output_file_prefix: str,
    primary_tissue: Optional[str] = None,
    threshold: Union[float, List[float]] = config.DEFAULT_MIN_PROB_THRESHOLD,
    consensus_graph: Optional[Dict[str, float]] = None,
) -> None:
    """
    Plot the thresholded consensus migration graph with collapsed multiedges.

    Args:
        data: DataFrame containing the log data
        output_file_prefix: Prefix for output files. Each file will be named {prefix}_{threshold}.pdf
        primary_tissue: Primary tissue label for migration analysis
        threshold: Single threshold value or list of thresholds for filtering edges
        consensus_graph: Pre-computed consensus graph
    """
    # Convert single threshold to list for uniform handling
    thresholds = [threshold] if isinstance(threshold, (int, float)) else threshold

    # Ensure all thresholds are float
    thresholds = [float(threshold) for threshold in thresholds]

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
        for node, color in zip(all_tissues, config.DEFAULT_COLORS[0 : len(all_tissues)])
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
                shape=config.DEFAULT_NODE_STYLES["shape"],
                fillcolor="white",
                penwidth=config.DEFAULT_NODE_STYLES["line_width"],
                fontsize=config.DEFAULT_FONT_SIZE,
            )

        for edge, probability in consensus_graph.items():
            if float(probability) > current_threshold:
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
                        penwidth=config.DEFAULT_NODE_STYLES["line_width"],
                        label="1",
                        fontsize=config.DEFAULT_FONT_SIZE,
                    )

        # Set empty label for single edges
        for source, target, data in G.edges(data=True):
            if data.get("label") == "1":
                data["label"] = ""

        dot = nx.nx_pydot.to_pydot(G)
        out_threshold = str(int(current_threshold * 100))
        output_file = f"{output_file_prefix}_{out_threshold}.pdf"
        dot.write_pdf(output_file)


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
        # Extract location from annotation
        match = re.search(r'&location="([^"]+)"', node.name)
        if match:
            node.tissue = match.group(1)
            # Keep original name without annotation
            node.name = node.name.split("[")[0]

        all_tissues.add(node.tissue)

    dists = {
        round(node.get_distance(tree), 5) for node in tree.traverse() if node.is_leaf()
    }

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
        for node, color in zip(all_tissues, config.DEFAULT_COLORS[0 : len(all_tissues)])
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
            shape=config.DEFAULT_NODE_STYLES["shape"],
            fillcolor="white",
            penwidth=config.DEFAULT_NODE_STYLES["line_width"],
            fontsize=config.DEFAULT_FONT_SIZE,
        )

    for edge, count in migration_counts.items():
        source, target = edge.split("_")
        G.add_edge(
            source,
            target,
            color=f'"{custom_colors[source]};0.5:{custom_colors[target]}"',
            penwidth=config.DEFAULT_NODE_STYLES["line_width"],
            label="" if count == 1 else str(count),
            fontsize=config.DEFAULT_PLOT_STYLES["legend_fontsize"],
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

    if metastasis_times:
        length, height = config.DEFAULT_FIGURE_SIZE
        fig, ax = plt.subplots(figsize=(length, height / 2))
        for migration, times in metastasis_times.items():
            source, target = migration.split("_")
            for time in times:
                ax.plot(
                    [time, time],
                    [0.5, 1],
                    color=custom_colors[source],
                    linewidth=config.DEFAULT_PLOT_STYLES["linewidth"],
                )
                ax.plot(
                    [time, time],
                    [0, 0.5],
                    color=custom_colors[target],
                    linewidth=config.DEFAULT_PLOT_STYLES["linewidth"],
                )

        ax.set_xlim(0, total_time)
        ax.set_ylim(0, 1)
        ax.set_xlabel("Time", fontsize=config.DEFAULT_FONT_SIZE)
        ax.set_yticks([0.25, 0.75])
        ax.set_yticklabels(["Target", "Source"], fontsize=config.DEFAULT_FONT_SIZE)
        ax.tick_params(axis="x", labelsize=config.DEFAULT_FONT_SIZE)

        handles = [
            plt.Line2D(
                [0], [0], color=color, lw=config.DEFAULT_PLOT_STYLES["linewidth"]
            )
            for color in custom_colors.values()
        ]
        ax.legend(
            handles,
            custom_colors.keys(),
            bbox_to_anchor=config.DEFAULT_PLOT_STYLES["legend_position"],
            loc="upper left",
            borderaxespad=0.0,
            frameon=False,
            fontsize=config.DEFAULT_PLOT_STYLES["legend_fontsize"] - 8,
        )

        plt.tight_layout()
        plt.savefig(f"{output_prefix}_migration_timing_{tree_num}.pdf")
        plt.close()

    # Plot tree
    os.environ["QT_QPA_PLATFORM"] = "offscreen"
    ts = TreeStyle()
    ts.rotation = config.DEFAULT_TREE_STYLE["rotation"]
    ts.scale = config.DEFAULT_TREE_STYLE["scale"]
    ts.show_leaf_name = config.DEFAULT_TREE_STYLE["show_leaf_name"]
    ts.show_branch_length = config.DEFAULT_TREE_STYLE["show_branch_length"]
    ts.show_border = config.DEFAULT_TREE_STYLE["show_border"]
    ts.show_scale = config.DEFAULT_TREE_STYLE["show_scale"]
    ts.mode = config.DEFAULT_TREE_STYLE["mode"]

    # Add legend
    for tissue in all_tissues:
        if tissue is not None:  # Skip None tissue type
            ts.legend.add_face(
                CircleFace(config.DEFAULT_NODE_STYLES["size"], custom_colors[tissue]),
                column=0,
            )
            ts.legend.add_face(
                TextFace(tissue, fsize=config.DEFAULT_FONT_SIZE), column=1
            )

    # Set node styles
    for node in origin.traverse():
        nstyle = NodeStyle()
        nstyle["shape"] = "circle"
        nstyle["size"] = config.DEFAULT_NODE_STYLES["size"]
        nstyle["hz_line_color"] = custom_colors[node.tissue]
        nstyle["vt_line_color"] = custom_colors[node.tissue]
        nstyle["hz_line_width"] = config.DEFAULT_NODE_STYLES["line_width"]
        nstyle["vt_line_width"] = config.DEFAULT_NODE_STYLES["line_width"]
        nstyle["fgcolor"] = custom_colors[node.tissue]
        node.set_style(nstyle)

    origin.render(f"{output_prefix}_tree_{tree_num}.pdf", tree_style=ts)


def plot_metastasis_timing(
    met_times: Dict[str, Dict[str, Tuple[float, float]]],
    consensus_graph: Dict[str, float],
    origin_time: float,
    primary_tissue: str,
    output_prefix: str,
    min_prob_threshold: float = config.DEFAULT_MIN_PROB_THRESHOLD,
) -> None:
    """
    Plot metastasis timing distributions.

    Args:
        met_times: Dictionary mapping tree labels to dictionaries of metastasis events
        consensus_graph: Dictionary mapping migration patterns to their probabilities
        origin_time: Time of origin
        primary_tissue: Primary tissue label
        output_prefix: Prefix for output files
        min_prob_threshold: Minimum probability threshold for migrations to include in plots
    """

    allowable_migrations = set()
    for migration, prob in consensus_graph.items():
        if float(prob) >= min_prob_threshold:
            allowable_migrations.add(migration)
    migration_counts = defaultdict(lambda: np.zeros(int(origin_time) + 1))
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

    df = pd.DataFrame(migration_counts, index=np.arange(0, int(origin_time) + 1)).T

    # Split the migration strings into source and target tissues
    df.index = pd.MultiIndex.from_tuples(
        [
            tuple([migration.split("_")[0], "_".join(migration.split("_")[1:])])
            for migration in df.index
        ],
        names=["source", "target"],
    )

    remaining_sources = sorted(
        set(df.index.get_level_values("source").unique()) - {primary_tissue}
    )
    remaining_sources.insert(0, primary_tissue)

    # Get unique tissues
    target_tissues = sorted(set(df.index.get_level_values("target")))
    target_tissues_no_num = set([tissue.split("_")[0] for tissue in target_tissues])
    tissues = sorted(
        set(df.index.get_level_values("source")).union(target_tissues_no_num)
        - {primary_tissue}
    )
    tissues.insert(0, primary_tissue)

    # Create a color palette for the tissues
    all_tissues = sorted(list(set(tissues) - {primary_tissue}))
    custom_colors = config.DEFAULT_COLORS
    custom_colors = {
        node: color
        for node, color in zip(all_tissues, custom_colors[0 : len(all_tissues)])
    }
    custom_colors[primary_tissue] = "black"

    # Create a grid of subplots with one row per source tissue
    if len(remaining_sources) == 1:
        fig, ax = plt.subplots(figsize=(15, 3))
        axes = [ax]  # Wrap single axis in list for consistent handling
    else:
        fig, axes = plt.subplots(
            len(remaining_sources),
            1,
            figsize=(15, 3 * len(remaining_sources)),
            sharex=True,
            sharey=True,
        )

    for i, source in enumerate(remaining_sources):
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
                    ax=axes[i],
                    color=custom_colors[target_name],
                )
                axes[i].fill_between(
                    df.columns,
                    df.loc[(source, target)],
                    alpha=0.3,
                    color=custom_colors[target_name],
                )

        axes[i].set_ylabel(source, fontsize=config.DEFAULT_FONT_SIZE)
        axes[i].tick_params(axis="both", which="major", labelsize=config.DEFAULT_FONT_SIZE)
        # ax.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)
        if i == len(remaining_sources) - 1:
            axes[i].set_xlabel("Time", fontsize=config.DEFAULT_FONT_SIZE)

    # Create a single legend for all axes
    handles = [
        plt.Line2D([0], [0], color=color, lw=2)
        for tissue, color in custom_colors.items()
    ]
    labels = list(custom_colors.keys())
    fig.legend(
        handles,
        labels,
        bbox_to_anchor=(1.05, 0.75),
        title="Target Tissue",
        frameon=False,
        fontsize=config.DEFAULT_FONT_SIZE,
        title_fontsize=config.DEFAULT_FONT_SIZE,
    )

    fig.text(
        0.001,
        0.5,
        "Source tissue",
        va="center",
        ha="center",
        rotation="vertical",
        fontsize=config.DEFAULT_FONT_SIZE,
    )

    plt.tight_layout(rect=[0.02, 0, 0.88, 1])
    plt.savefig(f"{output_prefix}_metastasis_timing_prob.pdf", bbox_inches="tight")
    plt.close()

    # plot midpoints
    if len(remaining_sources) == 1:
        fig, ax = plt.subplots(figsize=(15, 3))
        axes = [ax]  # Wrap single axis in list for consistent handling
    else:
        fig, axes = plt.subplots(
            len(remaining_sources),
            1,
            figsize=(15, 3 * len(remaining_sources)),
            sharex=True,
            sharey=True,
        )

    for i, source in enumerate(remaining_sources):
        for target in target_tissues:
            migration = f"{source}_{target}"
            if migration in migration_counts_mid_points:
                if "_1" in target:
                    target_reformatted = target.split("_")[0]
                else:
                    target_reformatted = target
                target_name = target.split("_")[0]
                axes[i].hist(
                    migration_counts_mid_points[migration],
                    bins=100,
                    color=custom_colors[target_name],
                    alpha=0.6,
                    label=target_reformatted,
                )
        axes[i].set_ylabel(source, fontsize=config.DEFAULT_FONT_SIZE)
        axes[i].tick_params(axis="both", which="major", labelsize=config.DEFAULT_FONT_SIZE)
        if i == len(remaining_sources) - 1:
            axes[i].set_xlabel("Time", fontsize=config.DEFAULT_FONT_SIZE)

    fig.legend(
        handles,
        labels,
        bbox_to_anchor=(1.05, 0.75),
        title="Target Tissue",
        frameon=False,
        fontsize=config.DEFAULT_FONT_SIZE,
        title_fontsize=config.DEFAULT_FONT_SIZE,
    )
    fig.text(
        0.001,
        0.5,
        "Source tissue",
        va="center",
        ha="center",
        rotation="vertical",
        fontsize=config.DEFAULT_FONT_SIZE,
    )

    plt.tight_layout()
    plt.savefig(f"{output_prefix}_metastasis_timing_midpoints.pdf", bbox_inches="tight")
    plt.close()


def plot_rate_matrix(
    data: pd.DataFrame,
    output_file: str,
    primary_tissue: str,
) -> None:
    """
    Plot the tissue substitution rate matrix as a heatmap.

    Args:
        data (pd.DataFrame): The log data containing tissue substitution rates
        output_file (str): Path to save the plot
        primary_tissue (str): Primary tissue label
    """
    # Get tissue rate column names
    tissue_rate_col_names = [
        name for name in data.columns if name.startswith("tissueSubstModelLogger")
    ]

    if len(tissue_rate_col_names) == 0:
        raise ValueError("No tissue substitution rate columns found in log data")
    
    # Get unique tissues
    tissues = list(
        set(
            [
                tissue
                for name in tissue_rate_col_names
                for tissue in name.replace("tissueSubstModelLogger.relGeoRate_", "").split(
                    "_"
                )
            ]
        )
    )
    tissues = [primary_tissue] + sorted([tis for tis in tissues if tis != primary_tissue])

    # Create rate matrix
    num_tissues = len(tissues)
    rate_matrix = np.zeros((num_tissues, num_tissues))

    for source in tissues:
        for recipient in tissues:
            if source == recipient:
                continue
            rate = data[f'tissueSubstModelLogger.relGeoRate_{source}_{recipient}'].mean()
            i = tissues.index(source)
            j = tissues.index(recipient)
            rate_matrix[i, j] = rate

    # Plot the rate matrix
    length, height = config.DEFAULT_FIGURE_SIZE
    plt.figure(figsize=(height, height))
    heatmap = sns.heatmap(
        rate_matrix,
        xticklabels=tissues,
        yticklabels=tissues,
        annot=True,
        cmap="YlOrRd",
        annot_kws={"size": config.DEFAULT_FONT_SIZE},
        cbar_kws={"label": "Rate"},
    )
    heatmap.figure.axes[-1].yaxis.label.set_size(config.DEFAULT_FONT_SIZE)
    heatmap.figure.axes[-1].tick_params(labelsize=config.DEFAULT_FONT_SIZE)
    plt.xlabel("Recipient", fontsize=config.DEFAULT_FONT_SIZE)
    plt.ylabel("Source", fontsize=config.DEFAULT_FONT_SIZE)
    plt.xticks(fontsize=config.DEFAULT_FONT_SIZE)
    plt.yticks(fontsize=config.DEFAULT_FONT_SIZE)
    plt.tight_layout()

    plt.savefig(output_file)
    plt.close()
