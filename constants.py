# TODO: Migrate all constants from other files to this file.
# TODO: Migrate hardcoded values to the config file.
import yaml

# ======================= HELPER FUNCTIONS =======================

# NOTE: This config file is used for both the pulldown pipeline and the individual complex feature prediction, depending on which script is run.
CONFIG_FILE = "configs/complex/PxdA.yaml"

def load_yaml_config(file_path):
    with open(file_path, "r") as file:
        return yaml.safe_load(file)


CONFIG = load_yaml_config(CONFIG_FILE)
CURRENT_PIPELINE = CONFIG["current_pipeline"]

# ======================= CONFIGURATION =======================

# Compute-specific constants:
USE_MULTIPROCESSING = True
PROCESS_COUNT = 6
MAX_ID_LENGTH = 40

# Input & AF3 constants:
DATABASE_FOLDER = "$HOME/public_databases"
MODEL_WEIGHTS_FOLDER = "$HOME/"
NUMBER_OF_SEEDS = CONFIG["number_of_seeds"]
MAX_COMBINED_SEQ_LENGTH = CONFIG.get("max_combined_seq_length", 10000)
TEMPLATE_FILE = CONFIG["template_file"]

# ================== Pulldown pipeline specific, not used in individual complex feature ==================
if CURRENT_PIPELINE == "pulldown":
    BAIT_FILENAME = CONFIG["bait_filename"]
    SUPERFOLDER_TO_FASTA_AND_FOLDER = CONFIG["superfolder_to_fasta_and_folder"]
    SUPERFOLDER_TO_FASTA = {k: [fasta for fasta, _ in fasta_and_folder] for k, fasta_and_folder in SUPERFOLDER_TO_FASTA_AND_FOLDER.items()}
    FASTA_FILES = sum(SUPERFOLDER_TO_FASTA.values(), [])
    SUPERFOLDER_TO_FOLDER = {k: [folder for _, folder in fasta_and_folder] for k, fasta_and_folder in SUPERFOLDER_TO_FASTA_AND_FOLDER.items()}
    FOLDERS = sum(SUPERFOLDER_TO_FOLDER.values(), [])
    # The model file that contains a larger region (or the entire) prey protein that were not included in the input to the bulkalphafold run;
    # Used to calculate clashes of the bait protein with other regions of the prey protein.
    # For example, with LRRK2, predictions were run using just the ROC-COR region, but the clashes_model is a pdb file containing the entire RCKW region.
    CLASHES_MODEL = CONFIG.get("clashes_model") 

    # Further downstream processing constants (thresholds are inclusive):
    DOMAINS_TO_RESIDUES = CONFIG["domains_to_residues"]
    # FULL LRRK2
    # DOMAINS_TO_RESIDUES = {
    #     "ARM": [1, 704],
    #     "ANK": [705, 799],
    #     "LRR": [800, 1331],
    #     "ROC": [1332, 1521],
    #     "COR-A": [1522, 1671],
    #     "COR-B": [1672, 1878],
    #     "KINASE": [1879, 2139],
    #     "WD40": [2140, 2527]
    # }
    ALL_PLDDT_WINDOWS = [-1, 11, 25]
    PLDDT_SLIDING_WINDOW = 11
    PREDICTION_THRESHOLD_METRIC_VALUE = 0.4
    PREDICTION_THRESHOLD_METRIC = "ipTM"
    MAX_CLASHES_THRESHOLD = 1000
    MIN_CONTACTS_THRESHOLDS = [1, 20, 80]

# ================== Indivdiual complex prediction specific, not used in pulldown pipeline ==================
if CURRENT_PIPELINE == "complex":
    FIXED_PROTEINS = CONFIG["fixed_proteins"]
    INPUT_FASTA = CONFIG["input_fasta"]
    OUTPUT_FOLDER = CONFIG["output_folder"]