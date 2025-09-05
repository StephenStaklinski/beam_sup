import os
import pandas as pd
import numpy as np
from Bio import Phylo
from Bio.Phylo.Newick import Tree


def format_beam_inputs(
    indel_matrix_file: str,
    outdir: str,
) -> None:
    """
    Formats the indel matrix into the required input files for BEAM.

    This function reads an indel matrix from a TSV file, converts mutation values to
    sequential integers, computes mutation proportions, and writes the necessary input files
    for BEAM, including the edit rate proportions, mutation dictionary, and FASTA file.

    Args:
        indel_matrix_file (str): Path to the input indel matrix file (TSV format).
        outdir (str): Directory where output files will be written.

    Outputs:
        - <outdir>/<sim_number>_edit_rate_proportions.txt: Proportions of each mutation type.
        - <outdir>/<sim_number>_original_mut_to_edit_model_mut.csv: Mapping of original mutations to new sequential integers.
        - <outdir>/<sim_number>.fasta: FASTA file representing the indel matrix with sequential integers.

    """
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # assumes the output dir is named with the sim number
    outname = os.path.basename(outdir)

    # read in the indel matrix
    indel_matrix = pd.read_csv(indel_matrix_file, sep="\t", index_col=0)

    # replace mutation values with sequential values required for the editing model
    done = []
    mut_dict = {-1: -1}  # keep as a dropout site to replace later
    i = 1
    for vals in indel_matrix.values.tolist():
        for v in vals:
            if v == -1 or v == 0 or v in mut_dict.keys():
                continue
            else:
                mut_dict[v] = i
                i = i + 1

    # replace all entries in the indel_matrix with the mut_dict value
    indel_matrix = indel_matrix.replace(mut_dict)

    # get all indel proportions
    muts = np.array(
        [v for vals in indel_matrix.values.tolist() for v in vals if v != 0 and v != -1]
    )  # unedited and silenced state is not included in the editRatePropostions calculations since they have their own free parameter input in Beam
    ordered_value_counts = np.unique(muts, return_counts=True)[1]
    sum_proportions = sum(ordered_value_counts)
    proportions = [str(count / sum_proportions) for count in ordered_value_counts]

    # replace the -1 with the largest value + 1 for dropout as the last column in tidetree
    max_val = max(mut_dict.values())
    mut_dict[-1] = max_val + 1
    indel_matrix = indel_matrix.replace(-1, mut_dict[-1])

    # write mutation proportions for initial states in the edit model
    outfile_proportions = f"{outdir}/{outname}_edit_rate_proportions.txt"
    with open(outfile_proportions, "w") as file:
        file.write(" ".join(proportions))

    # write mut dict to file
    outfile_mut_dict = f"{outdir}/{outname}_original_mut_to_edit_model_mut.csv"
    with open(outfile_mut_dict, "w") as file:
        file.write(f"original_mut,new_mut\n")
        for key, value in mut_dict.items():
            file.write(f"{key},{value}\n")

    # write fasta file
    outfile_fasta = f"{outdir}/{outname}.fasta"
    with open(outfile_fasta, "w") as file:
        for index, row in indel_matrix.iterrows():
            sequence = ",".join(str(x) for x in row)
            file.write(f">{index}\n{sequence}\n")


def get_fixed_edit_rates_for_beam_sim_matrix(
    original_rates_file: str, reordering_dict_file: str, outfile: str
) -> None:
    """
    Generates a normalized list of edit rates for use in BEAM simulation, reordered according to a mapping file.

    This function reads an original rates file containing code-to-rate mappings and a reordering dictionary file
    mapping new codes to old codes. It creates a new list of rates corresponding to the new codes, normalizes
    them so that their sum is 1.0 (as required by BEAM), and writes the resulting rates to an output file,
    one per line, in ascending order of the new codes.

    Args:
        original_rates_file (str): Path to the CSV file containing original code-to-rate mappings.
                                   The file should have a header and each line should be in the format: code,rate
        reordering_dict_file (str): Path to the CSV file mapping new codes to old codes.
                                    The file should have a header and each line should be in the format: new_code,old_code
        outfile (str): Path to the output file where the normalized, ordered rates will be written, one per line.

    Returns:
        None
    """
    # Get the original code to rate mapping
    original_rates = {}
    with open(original_rates_file) as f:
        # skip the header
        f.readline()
        for line in f:
            code, rate = line.strip().split(",")
            original_rates[code] = rate

    # Get the mapping of new code to old code
    new_rates = {}
    with open(reordering_dict_file) as f:
        # skip the header
        f.readline()
        for line in f:
            new_code, old_code = line.strip().split(",")
            # we dont need to map missing data code to a rate since beam only takes in fixed rates for mutation outcomes, not silenced sites
            if str(new_code) != "-1":
                new_rates[new_code] = original_rates[old_code]

    # Order the new_rates by the new codes in ascending order
    ordered_rates = [
        rate for code, rate in sorted(new_rates.items(), key=lambda item: int(item[0]))
    ]

    # Normalize the rates in case they are not already since beam requires the rates to sum to 1.0
    total_rates = sum([float(rate) for rate in ordered_rates])
    ordered_rates = [str(float(rate) / total_rates) for rate in ordered_rates]

    with open(outfile, "w") as f:
        for rate in ordered_rates:
            f.write(f"{rate}\n")
