import os
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from process_results import print_uniprot_details
from constants import SUPERFOLDER_TO_FOLDER, PLDDT_SLIDING_WINDOW, MAX_CLASHES_THRESHOLD, PREDICTION_THRESHOLD_METRIC, PREDICTION_THRESHOLD_METRIC_VALUE, DOMAINS_TO_RESIDUES, MIN_CONTACTS_THRESHOLDS
from get_binding_domain_combinations import generate_binding_domain_combinations

INPUT_FOLDER = "output_merged_results/"
OUTPUT_FOLDER = "output_analyze_results/"


def plot_binding_domain_combinations(folder, merged_df, filtered_df, threshold):
    # plot bar chart of all binding domain combinations frequency
    binding_domain_combinations = generate_binding_domain_combinations(DOMAINS_TO_RESIDUES.keys())
    binding_domain_combinations = ["_".join(combination) + f"_MIN_{threshold}" for combination in binding_domain_combinations]
    freq_merged = merged_df[binding_domain_combinations].sum()
    freq_filtered = filtered_df[binding_domain_combinations].sum()

    plt.figure(figsize=(10, 8))
    sns.barplot(x=freq_merged.index, y=freq_merged.values, alpha=0.5, label="All Data")
    sns.barplot(x=freq_filtered.index, y=freq_filtered.values, alpha=0.8, label="Filtered")
    plt.xlabel("Binding domain combinations")
    plt.xticks(rotation=45)
    plt.ylabel("Count")
    plt.yscale("log")
    ticks = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000]
    plt.yticks(ticks, [str(tick) for tick in ticks])
    plt.title(f"Binding domain combination frequency for {folder}")
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_FOLDER}{folder}_binding_domain_all_combinations_plddt_window_{PLDDT_SLIDING_WINDOW}_MIN_{threshold}.svg", bbox_inches="tight")

    plt.close("all")


def main():
    if not os.path.exists(OUTPUT_FOLDER):
        os.makedirs(OUTPUT_FOLDER)

    for superfolder, folders in SUPERFOLDER_TO_FOLDER.items():
        print(f"Processing {superfolder}", flush=True)
        merged_df = pd.DataFrame()
        for folder in folders:
            print(f"Processing {folder}", flush=True)

            merged_file = INPUT_FOLDER + folder + f"_merged_plddt_window_{PLDDT_SLIDING_WINDOW}.csv"
            output_file = OUTPUT_FOLDER + folder + f"_analysis_plddt_window_{PLDDT_SLIDING_WINDOW}.tsv"
            partial_merged_df = pd.read_csv(merged_file)
            merged_df = pd.concat([merged_df, partial_merged_df], ignore_index=True)

        # TODO: Adapt to not be hard-coded
        merged_df["binding domain"] = merged_df[["ROC_contacts", "COR-A_contacts", "COR-B_contacts"]].idxmax(axis=1)
        merged_df["binding domain"] = merged_df["binding domain"].map({"ROC_contacts": "ROC", "COR-A_contacts": "COR-A", "COR-B_contacts": "COR-B"})

        best_by_metric_idx = merged_df.groupby("fasta")[PREDICTION_THRESHOLD_METRIC].idxmax()
        top_by_metric_df = merged_df.loc[best_by_metric_idx].reset_index(drop=True)

        filtered_df = merged_df[(merged_df[PREDICTION_THRESHOLD_METRIC] >= PREDICTION_THRESHOLD_METRIC_VALUE) & (merged_df["between clashes"] <= MAX_CLASHES_THRESHOLD)]
        unique_uniprot_link = filtered_df["uniprot link"].unique()

        # print & write out filtered results
        # print(f"\n\n\n\n=================== {superfolder} ===================")
        # print(f"Filtered predictions for {superfolder}: {len(filtered_df)} models out of {len(merged_df)} models; {len(unique_uniprot_link)} unique proteins")
        # with open(output_file, "w") as f:
        #     f.write("uniprot\tprotein\torganism\thits\n")
        #     for uniprot_link in unique_uniprot_link:
        #         print()
        #         uniprot_details = print_uniprot_details(uniprot_link.split("/")[-2])
        #         hits = filtered_df[filtered_df["uniprot link"] == uniprot_link].shape[0]
        #         f.write(f"{uniprot_link}\t{uniprot_details['full_name']}\t{uniprot_details['organism']}\t{hits}\n")

        # plot distribution of prediction metric and clashes vs. RCKW
        g = sns.JointGrid(x=PREDICTION_THRESHOLD_METRIC, y="between clashes", data=top_by_metric_df, height=8)
        g.plot_marginals(sns.kdeplot)
        g.plot_joint(sns.scatterplot, alpha=0.5, edgecolor=None, s=16)
        # g.plot_joint(sns.scatterplot, alpha=0.5, edgecolor=None, s=10, hue=merged_df["binding domain"], palette="Set2")
        g.figure.suptitle(f"Scatter plot of {PREDICTION_THRESHOLD_METRIC} vs between clashes for {superfolder}", y=1.03)
        g.set_axis_labels(PREDICTION_THRESHOLD_METRIC, "between clashes")
        plt.savefig(f"{OUTPUT_FOLDER}{superfolder}_jointplot_{PREDICTION_THRESHOLD_METRIC}_vs_clashes_plddt_window_{PLDDT_SLIDING_WINDOW}.svg", bbox_inches="tight")

        # plot bar chart of binding domain frequency (calculated by the max number of contacts between the three binding domains)
        plt.figure(figsize=(10, 6))
        # TODO: Adapt to not be hard-coded
        binding_domain_order = ["ROC", "COR-A", "COR-B"]
        sns.countplot(data=merged_df, x="binding domain", alpha=0.5, label="All Data", order=binding_domain_order)
        sns.countplot(data=filtered_df, x="binding domain", alpha=0.8, label="Filtered", order=binding_domain_order)
        plt.title(f"Binding domain distribution for {superfolder}")
        plt.xlabel("Binding domain")
        plt.ylabel("Count")
        plt.yscale("log")
        ticks = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000]
        plt.yticks(ticks, [str(tick) for tick in ticks])
        plt.tight_layout()
        plt.savefig(f"{OUTPUT_FOLDER}{superfolder}_binding_domain_distribution_plddt_window_{PLDDT_SLIDING_WINDOW}.svg", bbox_inches="tight")

        for threshold in MIN_CONTACTS_THRESHOLDS:
            plot_binding_domain_combinations(superfolder, merged_df, filtered_df, threshold)


if __name__ == "__main__":
    main()
