"""Configuration settings for beam_visualization."""

# File parsing settings
DEFAULT_BURNIN_PERCENT = 0.0
DEFAULT_CORES = 1

# Plotting settings
DEFAULT_FIGURE_SIZE = (12, 6)
DEFAULT_FONT_SIZE = 18

# Data processing settings
MAX_TREES_TO_LOAD = 10000  # Maximum number of trees to load for memory management

# Visualization colors
DEFAULT_COLORS = [
    "#006400",  # Dark Green
    "#FF0000",  # Red
    "#0000CD",  # Medium Blue
    "#FFA500",  # Orange
    "#800080",  # Purple
    "#808080",  # Gray
    "#FFC0CB",  # Pink
    "#ADD8E6",  # Light Blue
    "#A52A2A",  # Brown
    "#FFFF00",  # Yellow
] * 3  # Repeat colors if needed for more tissues

# Default probability threshold for migrations
DEFAULT_MIN_PROB_THRESHOLD = 0.5

# Default number of trees to sample
DEFAULT_NUM_SAMPLES = 1

# Default plot styles
DEFAULT_PLOT_STYLES = {
    'linewidth': 2,
    'alpha': 0.3,
    'bins': 100,
    'legend_position': (1.05, 0.75),
    'legend_title': "Target Tissue",
    'legend_fontsize': 22,
    'legend_title_fontsize': 22
}

# Default node styles for tree visualization
DEFAULT_NODE_STYLES = {
    'shape': 'box',
    'size': 10,
    'line_width': 3,
    'line_color': 'black'
}

# Default tree style settings
DEFAULT_TREE_STYLE = {
    'rotation': 90,
    'scale': 1,
    'show_leaf_name': True,
    'show_branch_length': False,
    'show_border': False,
    'show_scale': False,
    'mode': 'r'
} 