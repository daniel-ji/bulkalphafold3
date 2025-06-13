"""
A script to compare the filtered proteins from two different pipeline runs (AF2 vs. AF3 or different downstream filtering after AF3). 
This script reads the protein names (or uniprot links) from two files (one for AF2 and one for AF3) and compares them to find unique and common proteins.
"""
import pandas as pd

from constants import FOLDERS

protein_list_1 = f"output_filtered_protein_list/GAP_RCKW.tsv"
protein_list_2 = f"output_filtered_protein_list/GAP_ROC_COR.tsv"

protein_list_1_df = pd.read_csv(protein_list_1, sep="\t", header=None)
protein_list_2_df = pd.read_csv(protein_list_2, sep="\t", header=None)

proteins_1 = protein_list_1_df[0].values
proteins_2 = protein_list_2_df[0].values

only_1 = set(proteins_1) - set(proteins_2)
only_2 = set(proteins_2) - set(proteins_1)
both = set(proteins_1) & set(proteins_2)

print(f"{protein_list_1} only: {len(only_1)}")
for protein in only_1:
    print(protein)
print(f"{protein_list_2} only: {len(only_2)}")
for protein in only_2:
    print(protein)
print(f"Both: {len(both)}")
for protein in both:
    print(protein)
