import os
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from process_results import print_uniprot_details
from constants import FOLDERS

output_folder = "output_analyze_results/"
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

WINDOW_SIZE = 11

for folder in FOLDERS:
    merged_file = folder + f"_merged_plddt_window_{WINDOW_SIZE}.csv"
    output_file = output_folder + folder + f"_analysis_plddt_window_{WINDOW_SIZE}.tsv"
    merged_df = pd.read_csv(merged_file)
    merged_df["binding area"] = merged_df[["ROC_contacts", "COR-A_contacts", "COR-B_contacts"]].idxmax(axis=1)
    merged_df["binding area"] = merged_df["binding area"].map({
        "ROC_contacts": "ROC",
        "COR-A_contacts": "COR-A",
        "COR-B_contacts": "COR-B"
    })

    best_iptm_idx = merged_df.groupby("fasta")["ipTM"].idxmax()
    top_ipTM_df = merged_df.loc[best_iptm_idx].reset_index(drop=True)

    filtered_df = merged_df[(merged_df["ipTM"] >= 0.4) & (merged_df["between clashes"] < 1000)]
    unique_uniprot_link = filtered_df["uniprot link"].unique()

    # print & write out filtered results
    print(f"\n\n\n\n=================== {folder} ===================")
    print(f"Filtered predictions for {folder}: {len(filtered_df)} models out of {len(merged_df)} models; {len(unique_uniprot_link)} unique proteins")
    with open(output_file, "w") as f:
        f.write("uniprot\tprotein\torganism\thits\n")
        for uniprot_link in unique_uniprot_link:
            print()
            uniprot_details = print_uniprot_details(uniprot_link.split("/")[-2])
            hits = filtered_df[filtered_df["uniprot link"] == uniprot_link].shape[0]
            f.write(f"{uniprot_link}\t{uniprot_details['full_name']}\t{uniprot_details['organism']}\t{hits}\n")

    # plot distribution of ipTM and clashes vs. RCKW
    g = sns.JointGrid(x="ipTM", y="between clashes", data=top_ipTM_df, height=8)
    g.plot_marginals(sns.kdeplot)
    g.plot_joint(sns.scatterplot, alpha=0.5, edgecolor=None, s=16)
    # g.plot_joint(sns.scatterplot, alpha=0.5, edgecolor=None, s=10, hue=merged_df["binding area"], palette="Set2")
    g.figure.suptitle(f"Scatter plot of ipTM vs between clashes for {folder}", y=1.03)
    g.set_axis_labels("ipTM", "between clashes")
    plt.savefig(f"{output_folder}{folder}_jointplot_ipTM_vs_clashes_plddt_window_{WINDOW_SIZE}.svg", bbox_inches='tight')

    # plot bar chart of binding area frequency (calculated by the max number of contacts between the three binding areas)
    plt.figure(figsize=(10, 6))
    binding_area_order = ["ROC", "COR-A", "COR-B"]
    sns.countplot(data=merged_df, x="binding area", alpha=0.5, label="All Data", order=binding_area_order)
    sns.countplot(data=filtered_df, x="binding area", alpha=0.8, label="Filtered", order=binding_area_order)
    plt.title(f"Binding area distribution for {folder}")
    plt.xlabel("Binding area")
    plt.ylabel("Count")
    plt.yscale("log")
    ticks = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000]
    plt.yticks(ticks, [str(tick) for tick in ticks])
    plt.tight_layout()
    plt.savefig(f"{output_folder}{folder}_binding_area_distribution_plddt_window_{WINDOW_SIZE}.svg", bbox_inches='tight')