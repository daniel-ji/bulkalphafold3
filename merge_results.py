import os
import pandas as pd
from constants import CURRENT_PIPELINE
assert CURRENT_PIPELINE == "pulldown", "This script is only for the pulldown pipeline"
from constants import SUPERFOLDER_TO_FOLDER, PLDDT_SLIDING_WINDOW, CLASHES_MODEL, MAX_CLASHES_THRESHOLD

OUTPUT_FOLDER = "output_merged_results/"
if not os.path.exists(OUTPUT_FOLDER):
    os.makedirs(OUTPUT_FOLDER)

for superfolder, folders in SUPERFOLDER_TO_FOLDER.items():
    print(f"Processing {superfolder}", flush=True)
    superfolder_df = pd.DataFrame()

    for folder in folders:
        print(f"Processing {folder}", flush=True)
        output_file = OUTPUT_FOLDER + folder + f"_merged_plddt_window_{PLDDT_SLIDING_WINDOW}.csv"
        predictions_file = "output_raw_results/" + folder + "_results.csv"
        predictions_df = pd.read_csv(predictions_file)
        if PLDDT_SLIDING_WINDOW > 0:
            predictions_df["model"] = predictions_df["model"].str.replace(".cif", f"_plddt_window_{PLDDT_SLIDING_WINDOW}.cif")

        binding_file = "output_binding_domain/" + folder + f"_binding_domain_plddt_window_{PLDDT_SLIDING_WINDOW}.csv"
        binding_df = pd.read_csv(binding_file)

        # conditional merging of clashes data, if it exists
        if CLASHES_MODEL and MAX_CLASHES_THRESHOLD is not None:
            clashes_file = "output_all_clashes/" + folder + f"_clashes_plddt_window_{PLDDT_SLIDING_WINDOW}.csv"
            clashes_df = pd.read_csv(clashes_file)
            assert len(predictions_df) == len(clashes_df) == len(binding_df), "DataFrames have different lengths"
            merged_df = pd.merge(predictions_df, clashes_df, on="model")
        else:
            print("No clashes data available, merging only predictions and binding domain data")
            assert len(predictions_df) == len(binding_df), "DataFrames have different lengths"
            merged_df = pd.merge(predictions_df, binding_df, on="model")

        assert len(merged_df) == len(predictions_df), f"Merged DataFrame has different length than original DataFrames {len(merged_df)} vs {len(predictions_df)}"

        merged_df.to_csv(output_file, index=False)

        superfolder_df = pd.concat([superfolder_df, merged_df], ignore_index=True)
    
    superfolder_df.to_csv(OUTPUT_FOLDER + superfolder + "_merged_plddt_window_" + str(PLDDT_SLIDING_WINDOW) + ".csv", index=False)