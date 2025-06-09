# install chimerax for calculating clashes if needed
import os
from multiprocessing import Pool
from calculate_clashes import calculate_clashes
from files_helper import get_model_files
from constants import CURRENT_PIPELINE
assert CURRENT_PIPELINE == "pulldown", "This script is only for the pulldown pipeline"
from constants import FOLDERS, PLDDT_SLIDING_WINDOW, PROCESS_COUNT, CLASHES_MODEL, CONFIG_FILE

OUTPUT_FOLDER = "output_all_clashes/"

def process_folder(folder_name):
    csv_file = OUTPUT_FOLDER + folder_name + f"_clashes_plddt_window_{PLDDT_SLIDING_WINDOW}.csv"

    print(f"Processing {folder_name}", flush=True)

    with open(csv_file, "w") as f:
        f.write("model,input clashes,reference clashes,both clashes,between clashes\n")

    model_files = get_model_files(folder_name, residue_sliding_window=PLDDT_SLIDING_WINDOW)

    with Pool(processes=PROCESS_COUNT) as pool:
        results = pool.starmap(calculate_clashes, zip(model_files, [CLASHES_MODEL] * len(model_files)))

    with open(csv_file, "a") as f:
        for model_file, (input_only_clashes, reference_only_clashes, both_clashes, between_clashes) in zip(model_files, results):
            f.write(f"{model_file},{input_only_clashes},{reference_only_clashes},{both_clashes},{between_clashes}\n")

    print(f"Done processing {folder_name}", flush=True)


if __name__ == "__main__":
    if not os.path.exists(OUTPUT_FOLDER):
        os.makedirs(OUTPUT_FOLDER)

    if not CLASHES_MODEL:
        raise ValueError(f"CLASHES_MODEL is not set. Please set it in your YAML configuration file ({CONFIG_FILE}).")

    for folder in FOLDERS:
        process_folder(folder)
