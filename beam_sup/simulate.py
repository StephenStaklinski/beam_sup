import os
import numpy as np
from ete3 import Tree
from cassiopeia.data import CassiopeiaTree
from cassiopeia.sim import Cas9LineageTracingDataSimulator
import shutil
import subprocess
import glob

from .tree_utils import get_num_migrations_from_nwk_and_labeling


def simulate_metastatic_cancer_population(
    outdir: str,
    num_generations: int = 250,
    migration_rate: str = "1e-6",
    num_cells_downsample: int = 50,
    max_anatomical_sites: int = -1,
    rerun_no_migrations: bool = True,
    seed: int | None = None,
) -> int:
    """
    Simulates a metastatic cancer population using an external 'simulate' executable.

    Args:
        outdir (str): Output directory for simulation results.
        num_generations (int, optional): Number of generations to simulate. Default is 250.
        migration_rate (str or float, optional): Migration rate between anatomical sites. Default is "1e-6".
        num_cells_downsample (int, optional): Number of cells to downsample in the output. Default is 50.
        max_anatomical_sites (int, optional): Maximum number of anatomical sites. Default is -1 (no limit).
        rerun_no_migrations (bool, optional): If True, rerun the simulation if no migration events occured. Default is True.
        seed (int | None, optional): Random seed for reproducibility. If None, a random seed is generated.

    Raises:
        FileNotFoundError: If the 'simulate' executable is not found in PATH.

    Returns:
        int: The random seed used for the simulation.
    """
    simulate_exe = shutil.which("simulate")

    if simulate_exe is None:
        raise FileNotFoundError(
            "Simulator executable 'simulate' not found in PATH. Please ensure the cpp code is compiled and added to PATH."
        )

    if seed is None:
        seed = np.random.randint(0, 2**31 - 1)

    os.makedirs(outdir, exist_ok=True)

    outprefix = f"{outdir}/{seed}" if str(seed) not in outdir else outdir
    os.makedirs(outprefix, exist_ok=True)

    cmd = [
        simulate_exe,
        "-c",
        str(num_generations),
        "-mig",
        str(migration_rate),
        "-s",
        str(seed),
        "-o",
        outprefix,
        "-m",
        str(max_anatomical_sites),
        "-d",
        str(num_cells_downsample),
    ]

    logfile = os.path.join(outprefix, f"{seed}_terminal.log")
    with open(logfile, "w") as logf:
        subprocess.run(cmd, check=True, stdout=logf, stderr=subprocess.STDOUT)

    if rerun_no_migrations:
        # Check if migrations occurred by inspecting the migration graph files
        # If no migrations, then repeat simulation
        nwk = glob.glob(os.path.join(outprefix, "*.nwk"))[0]
        labeling = glob.glob(os.path.join(outprefix, "*.vertex.labeling"))[0]
        migrations = get_num_migrations_from_nwk_and_labeling(nwk, labeling)
        if int(migrations) == 0:
            print(
                f"No migrations detected in simulation with seed {seed}. Rerunning simulation."
            )
            shutil.rmtree(outprefix)
            # Rerun simulation with a new random seed
            return simulate_metastatic_cancer_population(
                outdir,
                num_generations=num_generations,
                migration_rate=migration_rate,
                num_cells_downsample=num_cells_downsample,
                max_anatomical_sites=max_anatomical_sites,
                rerun_no_migrations=rerun_no_migrations,
                seed=None,
            )

    return seed


def overlay_simulated_crispr_barcode_data(
    ground_truth_tree_filepath: str,
    outprefix: str,
    num_sites: int = 50,
    mutationrate: float = 0.1,
) -> None:
    """
    Simulates Cas9 lineage tracing data and overlays it onto a given ground truth tree.

    Args:
        ground_truth_tree_filepath (str): File path to the ground truth tree in Newick format.
        outprefix (str): Output prefix.
        num_cuts (int): Number of cassettes (cuts) to simulate.
        m (float): Mutation rate per cassette.

    Output:
        A TSV file containing the simulated indel character matrix is saved to the specified output directory.
    """
    mut_rates = [float(mutationrate)] * num_sites

    ete_tree = Tree(ground_truth_tree_filepath, format=3)
    branch_length_dict = {}
    for node in ete_tree.traverse():
        if not node.is_root():
            parent_node = node.up
            branch_length = node.dist
            node_names = (parent_node.name, node.name)
            branch_length_dict[node_names] = branch_length

    ground_truth_tree = CassiopeiaTree(tree=ete_tree)
    ground_truth_tree.set_branch_lengths(branch_length_dict)

    # Set indel priors manually to have actual rates later on, rather than always drawing them from scratch
    state_priors = {
        1: 0.001891385794602738,
        2: 0.0024762237066947427,
        3: 0.006434777250119938,
        4: 0.022626187451291142,
        5: 0.0004519611853084519,
        6: 0.004252307867693467,
        7: 0.008695199487743388,
        8: 0.005287228198826541,
        9: 0.00852515654186136,
        10: 0.0028375427438811205,
        11: 0.01630715871347483,
        12: 0.004283226603781587,
        13: 0.0067464292466787685,
        14: 0.016033678821811725,
        15: 0.01917491749069289,
        16: 0.0036585141765295664,
        17: 0.004891872021313965,
        18: 0.00867199597266662,
        19: 0.011140999534685158,
        20: 0.010581282853573013,
        21: 0.004541688131273875,
        22: 0.0007617862007893749,
        23: 0.036165552581349764,
        24: 0.007900631443893674,
        25: 0.00019973184290800428,
        26: 0.009877241632563307,
        27: 0.005210678254720317,
        28: 0.004181565880363428,
        29: 0.017589109966820247,
        30: 0.0017611954710861511,
        31: 0.00022959343560415528,
        32: 0.002994537359426538,
        33: 0.000801068502385592,
        34: 0.014061233100234731,
        35: 0.015163857946017672,
        36: 0.04639414188685866,
        37: 0.001504774488150319,
        38: 0.014571829241040163,
        39: 0.00011550294775905987,
        40: 0.0017412028611578365,
        41: 0.0019912235849203347,
        42: 0.015435544626458327,
        43: 0.0026174926141904324,
        44: 0.03034937733094538,
        45: 0.005414023081912718,
        46: 0.015567352766629855,
        47: 0.0018764294412664411,
        48: 0.03822133912883424,
        49: 0.005016017821502824,
        50: 0.017923582049993315,
        51: 0.00707709038325807,
        52: 0.01894418099476499,
        53: 0.006552348025361078,
        54: 0.011483251612581855,
        55: 0.00213045546278859,
        56: 0.00044723689376124275,
        57: 0.028216265577151605,
        58: 0.007500191290250902,
        59: 9.64495624736785e-06,
        60: 0.014456941979911852,
        61: 0.0015131231024375432,
        62: 0.04034276177940931,
        63: 0.007743384492752599,
        64: 0.0012231748469220054,
        65: 0.018702920184101884,
        66: 0.0008980883166958691,
        67: 0.005709122901667588,
        68: 0.014609172657832641,
        69: 0.008227991029964507,
        70: 0.00027099882320938255,
        71: 0.006800810951957274,
        72: 0.0016024909826768872,
        73: 0.002702309545993545,
        74: 0.013068054552918141,
        75: 0.021222580867047328,
        76: 0.006778296749981686,
        77: 0.005021657046780748,
        78: 0.03394700190346978,
        79: 0.01064798652834624,
        80: 0.0015709747044370804,
        81: 0.006723215374919598,
        82: 0.00918494570764104,
        83: 0.0007148514517145434,
        84: 0.033386853032898484,
        85: 0.004952562215118869,
        86: 0.0016703080977210618,
        87: 0.008780194621411974,
        88: 0.0021936132072342554,
        89: 0.0037054891861055513,
        90: 0.013720737264473232,
        91: 0.02129789157325353,
        92: 0.013693745961487449,
        93: 0.018491991117775546,
        94: 0.0011056635477868833,
        95: 0.0008658270062053632,
        96: 0.008272143321129521,
        97: 0.001439189221306065,
        98: 0.011779293293448219,
        99: 0.032353702218426296,
        100: 0.021099922150975396,
    }

    lt_sim = Cas9LineageTracingDataSimulator(
        number_of_cassettes=num_sites,
        size_of_cassette=1,
        mutation_rate=mut_rates,
        state_generating_distribution=lambda: np.random.exponential(1e-5),
        number_of_states=100,
        state_priors=state_priors,
        heritable_silencing_rate=0.0001,
        stochastic_silencing_rate=0.01,
        heritable_missing_data_state=-1,
        stochastic_missing_data_state=-1,
    )

    lt_sim.overlay_data(ground_truth_tree)
    final_matrix = ground_truth_tree.character_matrix
    out_matrix = outprefix + "_indel_character_matrix.tsv"
    final_matrix.to_csv(out_matrix, sep="\t")
