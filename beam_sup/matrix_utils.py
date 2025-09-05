import os
from typing import Sequence
import pandas as pd


def count_informative_characters(site_values: pd.Series) -> int:
    """
    Count the number of phylogenetically informative characters in a site.

    A character is considered informative if:
      1. It is not 0 (unedited) or -1 (missing data).
      2. It appears in more than one cell and less than all cells (i.e., provides phylogenetic information).

    Args:
        site_values (pd.Series): A pandas Series of character values for a site.

    Returns:
        int: The number of unique informative characters that appear more than once and less than the total number of cells.
    """
    # Remove 0s and -1s
    informative_values = site_values[~site_values.isin([0, -1])]
    # Count unique values and their frequencies
    value_counts = informative_values.value_counts()
    # Count how many unique values appear more than once and less than all cells
    informative_count = len(
        value_counts[(value_counts > 1) & (value_counts < len(site_values))]
    )

    return informative_count


def calculate_avg_informative_characters(input_file: str) -> float:
    """
    Reads a character matrix from a file and calculates the average number of informative characters per cell.

    Args:
        input_file (str): Path to the input file (CSV or TSV).

    Returns:
        float: The average number of informative characters per cell.
    """
    ext = os.path.splitext(input_file)[1].lower()
    sep = "\t" if ext == ".tsv" else ","
    df = pd.read_csv(input_file, sep=sep, index_col=0)

    informative_counts = df.apply(count_informative_characters)
    num_informative = sum(informative_counts)
    num_cells = df.shape[0]

    if num_informative == 0:
        return 0.0

    return num_informative / num_cells


def convert_matrix_to_row_successive_matrix(
    character_matrix: pd.DataFrame, mut_dict: dict, indel_priors: pd.DataFrame = None
) -> tuple[pd.DataFrame, dict, dict]:
    """
    Converts a character matrix to a successive character matrix, mapping mutations to unique successive integers.
    Also computes normalized edit rates for each successive mutation.

    Args:
        character_matrix (pd.DataFrame): The input character matrix with clones as rows and sites as columns.
        mut_dict (dict): Dictionary mapping site indices to mutation strings.
        indel_priors (Optional: pd.DataFrame): DataFrame containing mutation prior counts, indexed by mutation string. Default is None.

    Returns:
        tuple: (successive_char_matrix, successive_mut_dict, successive_edit_rates)
            successive_char_matrix (pd.DataFrame): The transformed character matrix.
            successive_mut_dict (dict): Mapping from mutation string to successive integer.
            successive_edit_rates (dict): Normalized edit rates for each successive mutation integer.
    """
    successive_char_matrix = character_matrix.copy()
    successive_mut_dict = {}
    i = 1
    for clone, row in character_matrix.iterrows():
        for site, mut in row.items():
            mut = int(mut)
            # Skip undedited and missing sites
            if mut == 0 or mut == -1:
                continue
            mut_str = mut_dict[int(site[1:]) - 1][mut]
            # Replace the mutation with the successive mutation
            if mut_str not in successive_mut_dict:
                successive_mut_dict[mut_str] = i
                new_mut_value = i
                i += 1
            else:
                new_mut_value = successive_mut_dict[mut_str]
            successive_char_matrix.loc[clone, site] = new_mut_value

    successive_edit_rates = {}

    if indel_priors:
        # Compute normalized edit rates for the successive matrix values
        for mut_str, mut_int in successive_mut_dict.items():
            successive_edit_rates[mut_int] = indel_priors.loc[mut_str]["count"]

        total_count = sum(successive_edit_rates.values())
        for mut_int in successive_edit_rates:
            successive_edit_rates[mut_int] = (
                successive_edit_rates[mut_int] / total_count
            )

    return successive_char_matrix, successive_mut_dict, successive_edit_rates


def expand_clones_with_multiple_tissues(
    matrix_df: pd.DataFrame, tissues_df: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Expands clones that are associated with multiple tissues into separate rows for each tissue.

    Args:
        matrix_df (pd.DataFrame): Character matrix with clones as index.
        tissues_df (pd.DataFrame): DataFrame with columns 'group_name' and 'tissues'.

    Returns:
        tuple: (expanded_matrix_df, expanded_tissues_df)
    """
    # Only keep relevant columns
    tissues_df = tissues_df[["group_name", "tissues"]].copy()

    # Find clones with more than one tissue
    clones_with_multiple_tissues = tissues_df.loc[
        tissues_df["tissues"].str.contains(",", na=False), "group_name"
    ].values.tolist()

    new_matrix_rows = []
    new_tissues_rows = []
    clones_to_drop = []

    for clone in clones_with_multiple_tissues:
        tissues = (
            tissues_df.loc[tissues_df["group_name"] == clone, "tissues"]
            .values[0]
            .split(",")
        )
        for i, tissue in enumerate(tissues):
            new_row = matrix_df.loc[matrix_df.index == clone].values[0]
            new_matrix_rows.append(
                pd.Series(new_row, index=matrix_df.columns, name=f"{clone}_{i}")
            )
            new_tissues_rows.append({"group_name": f"{clone}_{i}", "tissues": tissue})
        clones_to_drop.append(clone)
        tissues_df = tissues_df.loc[tissues_df["group_name"] != clone]

    # Remove the original entries for the clone
    for clone in clones_to_drop:
        matrix_df = matrix_df.drop(index=clone)

    # Concatenate new rows to the original DataFrames
    matrix_df = pd.concat([matrix_df, pd.DataFrame(new_matrix_rows)])
    tissues_df = pd.concat(
        [tissues_df, pd.DataFrame(new_tissues_rows)], ignore_index=True
    )

    return matrix_df, tissues_df
