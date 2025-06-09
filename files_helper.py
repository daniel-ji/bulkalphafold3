import re
import os
import pandas as pd
from modify_mmcif_plddt import get_plddt_sliding_window_mmcif
import multiprocessing as mp
from constants import PROCESS_COUNT, MAX_ID_LENGTH


N_MODEL = -1  # -1 for all models
# N_MODEL = 2


def get_sequences_from_fasta(fasta_file):
    """
    Extract sequences from a FASTA file. Returns a list of protein IDs and a list of tuples (sequence ID, sequence).
    Returns `proteins` are a list of UniProt IDs that would be used to fetch the sequences.
    Returns `sequences` are a list of tuples, where each tuple contains a sequence ID and the corresponding sequence.
    """

    with open(fasta_file) as f:
        proteins = [line.rstrip() for line in f if line.strip()]
        sequences = []
        if proteins[0][0] == ">":
            id = ""
            running_sequence = ""
            for line in proteins:
                if line[0] == ">":
                    if running_sequence != "":
                        assert id != ""
                        print("Found sequence " + id)
                        sequences.append((id, running_sequence))
                    id = re.sub(r"\W", "_", line[1 : MAX_ID_LENGTH + 1])
                    running_sequence = ""
                else:
                    running_sequence += line

            if running_sequence != "":
                assert id != ""
                print("Found sequence " + id)
                sequences.append((id, running_sequence))
            proteins = []

    return proteins, sequences


def get_model_files(folder, residue_sliding_window, n_model=N_MODEL):
    """
    Get model files from the given folder. Returned files vary depending on the provided sliding window length for calculating per-residue PLDDT scores.
    """

    with mp.Pool(processes=PROCESS_COUNT) as pool:
        model_folders = [f for f in os.listdir(folder) if os.path.isdir(os.path.join(folder, f))]
        all_model_files = pool.starmap(process_model_folder, [(model_folder, folder, residue_sliding_window, n_model) for model_folder in model_folders])

    all_model_files = [item for sublist in all_model_files for item in sublist]  # flatten the list of lists

    print(f"Found {len(all_model_files)} model files in {folder}")
    return all_model_files


def process_model_folder(model_folder, folder, residue_sliding_window, n_model=N_MODEL):
    model_files = []
    model_root = os.path.join(folder, model_folder)
    ranking_file = [f for f in os.listdir(model_root) if f.endswith("_ranking_scores.csv")][0]
    ranking_csv = pd.read_csv(os.path.join(model_root, ranking_file))
    ranking_csv.sort_values(by=["ranking_score"], ascending=False, inplace=True)
    ranking_data = ranking_csv.head(n_model) if n_model > 0 else ranking_csv
    for _, row in ranking_data.iterrows():
        seed_folder = f"seed-{int(row['seed'])}_sample-{int(row['sample'])}"
        model_file = [f for f in os.listdir(os.path.join(model_root, seed_folder)) if f.endswith("model.cif")][0]
        model_file_path = os.path.join(model_root, seed_folder, model_file)
        sliding_window_model_file = get_plddt_sliding_window_mmcif(model_file_path, residue_sliding_window=residue_sliding_window)
        model_files.append(sliding_window_model_file)

    return model_files

#TODO: In the future, support this for the complex pipeline as well.
if __name__ == "__main__":
    from constants import CURRENT_PIPELINE

    assert CURRENT_PIPELINE == "pulldown", "This script is only for the pulldown pipeline"
    from constants import FOLDERS

    for folder in FOLDERS:
        print(f"Processing {folder}", flush=True)
        model_files = get_model_files(folder, residue_sliding_window=11, n_model=N_MODEL)
