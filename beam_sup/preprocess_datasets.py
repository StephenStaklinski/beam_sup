import os
import pandas as pd
import cassiopeia as cas

from .matrix_utils import convert_matrix_to_row_successive_matrix


def convert_quinn_allele_table_to_successive_matrix(
    infile: str,
    lineage: int,
    outdir: str,
) -> None:
    """
    Converts a Quinn-style allele table to a successive character matrix and writes output files.
    This function reads an allele table from a TSV file, filters it for a specified lineage group,
    computes empirical indel priors, and converts the filtered allele table into both an original
    and a successive character matrix. It then writes the resulting matrices, mutation priors, and
    mutation dictionaries to the specified output directory.

    Args:
        infile (str): Path to the input allele table file (TSV format).
        lineage (int): The lineage group identifier to filter the allele table.
        outdir (str): Directory where output files will be written.

    Outputs:
        - <outdir>/<lineage>_mutation_priors.txt: Successive edit rates for each mutation code.
        - <outdir>/<lineage>_original_character_matrix.tsv: Original character matrix for the lineage group.
        - <outdir>/<lineage>_original_chracter_int_to_mutation_string_dict.txt: Mapping of character integers to mutation strings for the original matrix.
        - <outdir>/<lineage>_successive_character_matrix.tsv: Successive character matrix for the lineage group.
        - <outdir>/<lineage>_successive_int_to_mutation_string_dict.txt: Mapping of character integers to mutation strings for the successive matrix.

    """
    # read in the provided allele table
    allele_table = pd.read_csv(
        infile,
        sep="\t",
        usecols=[
            "cellBC",
            "intBC",
            "r1",
            "r2",
            "r3",
            "allele",
            "LineageGroup",
            "sampleID",
            "readCount",
            "UMI",
        ],
    )

    group = allele_table[allele_table["LineageGroup"] == lineage]

    # get indel priors as per Cassiopeia docs
    indel_priors = cas.pp.compute_empirical_indel_priors(group)

    char_matrix_df, priors, mut_dict = cas.pp.convert_alleletable_to_character_matrix(
        group,
        missing_data_state="-1",
        allele_rep_thresh=0.95,
        mutation_priors=indel_priors,
    )

    successive_matrix, new_mut_dict, successive_edit_rates = (
        convert_matrix_to_row_successive_matrix(char_matrix_df, mut_dict, indel_priors)
    )

    os.makedirs(outdir, exist_ok=True)

    # write successive edit rates to a file
    successive_edit_rates = dict(sorted(successive_edit_rates.items()))
    with open(f"{outdir}/{lineage}_mutation_priors.txt", "w") as f:
        f.write(f"mutation_code,rate\n")
        for key, value in successive_edit_rates.items():
            f.write(f"{key},{value}\n")

    # write each lineage group's successive matrix to its own file
    char_matrix_df.to_csv(
        f"{outdir}/{lineage}_original_character_matrix.tsv",
        sep="\t",
        index=True,
        header=True,
    )

    with open(
        f"{outdir}/{lineage}_original_chracter_int_to_mutation_string_dict.txt", "w"
    ) as f:
        f.write(f"site_num,char_int,mut_str\n")
        for key, value in mut_dict.items():
            for k, v in value.items():
                f.write(f"{key},{k},{v}\n")

    successive_matrix.to_csv(
        f"{outdir}/{lineage}_successive_character_matrix.tsv",
        sep="\t",
        index=True,
        header=True,
    )

    with open(
        f"{outdir}/{lineage}_successive_int_to_mutation_string_dict.txt", "w"
    ) as f:
        f.write(f"successive_char_int,mut_str\n")
        for key, value in new_mut_dict.items():
            f.write(f"{value},{key}\n")


def convert_simeonov_barcode_hmid_to_matrix(input_file: str, output_file: str) -> None:
    """
    Converts a Simeonov et al. data hmid barcode file into an indel matrix format.

    This function reads a tab-separated input file containing barcode information, processes the data to
    generate a matrix representation, and writes the results to the specified output file. Additionally,
    it generates two supplementary files: one mapping barcodes to tissues and another mapping mutations.

    Args:
        input_file (str): Path to the input TSV file containing barcode data.
        output_file (str): Path to the output CSV file where the matrix will be saved. Supplementary files
            for tissues and mutation mapping will be generated based on this path.

    Outputs:
        - A matrix CSV file at `output_file`.
        - A tissues mapping text file with suffix `_tissues.txt`.
        - A mutation mapping text file with suffix `_mut_mapping.txt`.
    """
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    output_tissues_file = output_file.replace("_matrix.csv", "_tissues.txt")
    output_mut_mapping_file = output_file.replace("_matrix.csv", "_mut_mapping.txt")

    df = pd.read_csv(input_file, sep="\t", index_col=0)
    df.index = ["barcode" + str(idx) for idx in df.index]

    # Use cells per barcode to determine tissues
    tissues = {}
    for idx, row in df.iterrows():
        cells = [cell.strip() for cell in row["cells"].split(",")]
        tissue_set = set(cell.split("_")[0] for cell in cells)
        tissues[idx] = tissue_set

    df = df.drop(columns=["cells"])

    # Split barcode into sites
    hmid_split = df["hmid"].str.split("-", expand=True)
    num_cols = hmid_split.shape[1]
    hmid_split.columns = [f"r{i+1}" for i in range(num_cols)]
    df = pd.concat([df.drop(columns=["hmid"]), hmid_split], axis=1)

    # Replace UNKNOWN with -1 and NONE with 0
    df = df.replace("UNKNOWN", "-1")
    df = df.replace("NONE", "0")

    # Now replace mut string with integers successively while also assigning multiple site mut str to the first site with other sites as -1
    mut_mapping = {}
    next_mut_id = 1
    for idx, row in df.iterrows():
        prev_val = None
        for col in df.columns:
            val = row[col]
            # Skip unedited or missing data sites
            if val == "-1" or val == "0":
                prev_val = None
            # if val == prev_val, set to -1
            elif val == prev_val:
                df.at[idx, col] = -1
            elif val in mut_mapping:
                df.at[idx, col] = mut_mapping[val]
                prev_val = val
            else:
                mut_mapping[val] = next_mut_id
                df.at[idx, col] = next_mut_id
                next_mut_id += 1
                prev_val = val

    # Write all to output csv files
    with open(output_mut_mapping_file, "w") as f:
        f.write("successive_char_int,mut_str\n")
        for mut_str, mut_id in mut_mapping.items():
            f.write(f"{mut_id},{mut_str}\n")

    with open(output_tissues_file, "w") as f:
        f.write("group_name,tissues\n")
        for barcode, tissue_set in tissues.items():
            f.write(f"{barcode},{';'.join(sorted(tissue_set))}\n")

    df.to_csv(output_file, index=True, header=True)


def preprocess_yang_data(allele_filepath: str, outdir: str) -> None:
    """
    Preprocesses the Yang et al. dataset by converting the allele table into successive character matrices
    for each metastatic lineage. Outputs the successive character matrices, mutation dictionaries,
    and collapsing dictionaries to the specified output directory.
    
    Args:
        allele_filepath (str): Path to the input allele table file (TSV format).
        outdir (str): Directory where output files will be written.
        
    Outputs:
        - <outdir>/<MetFamily>_successive_char_matrix.txt: Successive character matrix for each metastatic lineage.
        - <outdir>/<MetFamily>_mut_dict.txt: Mapping of successive character integers to mutation strings for each lineage.
        - <outdir>/<MetFamily>_collapsed.txt: Collapsed character matrix with unique rows for each lineage.
        - <outdir>/<MetFamily>_collapsing_dict.txt: Dictionary mapping group names to original cell barcodes and tissues for each lineage.
    """
    os.makedirs(outdir, exist_ok=True)

    allele_df = pd.read_csv(allele_filepath, sep="\t", index_col=0)
    # Find metastatic mice
    # Make new column with tissue labels only
    allele_df["tissue"] = [
        "".join(filter(str.isalpha, name[2]))
        for name in allele_df["Tumor"].str.split("_")
    ]

    # Get met lineage names and tissues
    tumor_names = allele_df.groupby("MetFamily")["tissue"].unique()

    # Keep only lineages with both primary and met tissues
    tumor_names = tumor_names[
        tumor_names.apply(
            lambda x: any(name.startswith("T") for name in x)
            and any(not name.startswith("T") for name in x)
        )
    ]

    # keep only subset MetFamily mice
    allele_df = allele_df[allele_df["MetFamily"].isin(tumor_names.index)]

    # Make character matrix
    # Write mouse specific character matrix to file
    for MetFamily, df in allele_df.groupby("MetFamily"):

        # Calculate the threshold to make sure sites are only included if at least one cell has a different mutation than the rest, i.e. remove uninformative sites with 100% of the site as the same mutation
        cell_threshold = 1 - (1 / df["cellBC"].nunique())

        # Built in cassiopeia function to convert allele table to character matrix
        char_matrix = cas.pp.convert_alleletable_to_character_matrix(
            df, allele_rep_thresh=cell_threshold
        )
        char_matrix_df = char_matrix[0]
        mut_dict = char_matrix[2]

        # Convert char matrix to successive char matrix
        successive_char_matrix = char_matrix_df.copy()
        successive_mut_dict = {}
        i = 1
        for clone, row in char_matrix_df.iterrows():
            for site, mut in row.items():
                # skip unedited or missing sites
                mut = int(mut)
                if mut == 0 or mut == -1:
                    continue
                mut_str = mut_dict[int(site[1:]) - 1][mut]
                # replace the mutation with the successive mutation
                if mut_str not in successive_mut_dict:
                    successive_mut_dict[mut_str] = i
                    new_mut_value = i
                    i += 1
                else:
                    new_mut_value = successive_mut_dict[mut_str]
                successive_char_matrix.loc[clone, site] = new_mut_value

        # Rename columns of successive char matrix to be successive themselves
        successive_char_matrix.columns = [
            f"r{i}" for i in range(1, len(successive_char_matrix.columns) + 1)
        ]

        # Output successive char matrix for all cells
        mouse_outfile = f"{outdir}/{MetFamily}_successive_char_matrix.txt"
        successive_char_matrix.index.name = "cellBC"
        # successive_char_matrix.to_csv(mouse_outfile, sep="\t", index=True)

        # Write mutation dictionary to file
        mut_dict_outfile = mouse_outfile.replace(".txt", f"_mut_dict.txt")
        # with open(mut_dict_outfile, "w") as f:
        #     f.write(f"mut_id\tmut_str\n")
        #     for str, id in successive_mut_dict.items():
        #         f.write(f"{id}\t{str}\n")

        # Collapse the cells to only unique rows and output collapsing dict of cellBCs and tissue labels
        all_columns = successive_char_matrix.columns.tolist()
        sorted_char_matrix = successive_char_matrix.sort_values(by=all_columns)
        unique_rows = sorted_char_matrix.drop_duplicates(keep="first")
        group_names = [f"clone{i+1}" for i in range(len(unique_rows))]
        group_to_originals = {}
        group_to_tissues = {}
        for group_name, (_, unique_row) in zip(group_names, unique_rows.iterrows()):
            # Find all rows in sorted_char_matrix that match the unique_row
            original_row_names = sorted_char_matrix[
                sorted_char_matrix.eq(unique_row).all(axis=1)
            ].index.tolist()
            group_to_originals[group_name] = original_row_names
            original_tissues = set(
                df[df["cellBC"].isin(original_row_names)]["tissue"].values.tolist()
            )
            group_to_tissues[group_name] = original_tissues

        # Replace index names in unique_rows with the appropriate group name
        unique_rows.index = group_names

        print(
            MetFamily,
            f"cells: {len(successive_char_matrix)}",
            f"clones: {len(unique_rows)}",
            f"sites: {len(successive_char_matrix.columns)}",
        )

        # Write unique rows to file
        unique_rows_outfile = mouse_outfile.replace(".txt", f"_collapsed.txt")
        # unique_rows.to_csv(unique_rows_outfile, sep="\t", index=True)

        # Write collapsing dict of cellBCs and tissue labels to file
        collapsing_dict_outfile = mouse_outfile.replace(".txt", f"_collapsing_dict.txt")
        # with open(collapsing_dict_outfile, "w") as f:
        #     f.write(f"group_name\tcellBCs\ttissues\n")
        #     for group_name in group_names:
        #         cellBCs = ','.join(list(group_to_originals[group_name]))
        #         tissues = ','.join(list(group_to_tissues[group_name]))
        #         f.write(f"{group_name}\t{cellBCs}\t{tissues}\n")
