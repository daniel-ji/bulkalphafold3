import os
import pandas as pd
from constants import SUPERFOLDER_TO_FOLDER, PLDDT_SLIDING_WINDOW

output_folder = "output_merged_results/"
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

for superfolder, folders in SUPERFOLDER_TO_FOLDER.items():
    print(f"Processing {superfolder}", flush=True)
    superfolder_df = pd.DataFrame()

    for folder in folders:
        print(f"Processing {folder}", flush=True)
        predictions_file = "output_raw_results/" + folder + "_results.csv"

        if PLDDT_SLIDING_WINDOW > 0:
            clashes_file = "output_all_clashes/" + folder + f"_clashes_plddt_window_{PLDDT_SLIDING_WINDOW}.csv"
            binding_file = "output_binding_domain/" + folder + f"_binding_domain_plddt_window_{PLDDT_SLIDING_WINDOW}.csv"
            output_file = output_folder + folder + f"_merged_plddt_window_{PLDDT_SLIDING_WINDOW}.csv"
        else:
            clashes_file = "output_all_clashes/" + folder + "_clashes.csv"
            binding_file = "output_binding_domain/" + folder + "_binding_domain.csv"
            output_file = output_folder + folder + "_merged.csv"

        predictions_df = pd.read_csv(predictions_file)

        if PLDDT_SLIDING_WINDOW > 0:
            predictions_df["model"] = predictions_df["model"].str.replace(".cif", f"_plddt_window_{PLDDT_SLIDING_WINDOW}.cif")

        clashes_df = pd.read_csv(clashes_file)
        binding_df = pd.read_csv(binding_file)

        assert len(predictions_df) == len(clashes_df) == len(binding_df), "DataFrames have different lengths"

        merged_df = pd.merge(predictions_df, clashes_df, on="model")

        merged_df = pd.merge(merged_df, binding_df, on="model")

        assert len(merged_df) == len(predictions_df), f"Merged DataFrame has different length than original DataFrames {len(merged_df)} vs {len(predictions_df)}"

        merged_df.to_csv(output_file, index=False)

        superfolder_df = pd.concat([superfolder_df, merged_df], ignore_index=True)
    
    if PLDDT_SLIDING_WINDOW < 0:
        superfolder_df.to_csv(output_folder + superfolder + "_merged.csv", index=False)
    else:
        superfolder_df.to_csv(output_folder + superfolder + "_merged_plddt_window_" + str(PLDDT_SLIDING_WINDOW) + ".csv", index=False)