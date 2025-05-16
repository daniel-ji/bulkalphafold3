# install chimerax for calculating clashes if needed

from multiprocessing import Pool
from calculate_clashes import calculate_clashes
from files_helper import get_model_files
from constants import FOLDERS, PLDDT_SLIDING_WINDOW, PROCESS_COUNT

LRRK2_REFERENCE_PDB = "LRRK2_RCKW.pdb"


def process_folder(folder_name):
    csv_file = "output_all_clashes/" + folder_name + f"_clashes_plddt_window_{PLDDT_SLIDING_WINDOW}.csv"

    print(f"Processing {folder_name}", flush=True)

    with open(csv_file, "w") as f:
        f.write("model,input clashes,reference clashes,both clashes,between clashes\n")

    model_files = get_model_files(folder_name, residue_sliding_window=PLDDT_SLIDING_WINDOW)

    with Pool(processes=PROCESS_COUNT) as pool:
        results = pool.starmap(calculate_clashes, zip(model_files, [LRRK2_REFERENCE_PDB] * len(model_files)))

    with open(csv_file, "a") as f:
        for model_file, (input_only_clashes, reference_only_clashes, both_clashes, between_clashes) in zip(model_files, results):
            f.write(f"{model_file},{input_only_clashes},{reference_only_clashes},{both_clashes},{between_clashes}\n")

    print(f"Done processing {folder_name}", flush=True)


if __name__ == "__main__":
    for folder in FOLDERS:
        process_folder(folder)
