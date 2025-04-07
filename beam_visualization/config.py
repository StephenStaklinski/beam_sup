"""Configuration settings for beam_visualization."""

# File parsing settings
DEFAULT_BURNIN_PERCENT = 0.0
DEFAULT_CORES = 1

# Plotting settings
DEFAULT_FIGURE_SIZE = (12, 6)
DEFAULT_FONT_SIZE = 12

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