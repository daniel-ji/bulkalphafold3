"""
A script to compare the filtered proteins from two different pipeline runs (AF2 vs. AF3 or different downstream filtering after AF3). 
This script reads the protein names (or uniprot links) from two files (one for AF2 and one for AF3) and compares them to find unique and common proteins.
"""

from constants import FOLDERS

PROTEIN_LIST_1_SUFFIX = "_AF2_proteins.txt"
PROTEIN_LIST_2_SUFFIX = "_AF3_proteins.txt"

for folder in FOLDERS:
    protein_list_1 = f"output_filtered_protein_list/{folder}{PROTEIN_LIST_1_SUFFIX}"
    protein_list_2 = f"output_filtered_protein_list/{folder}{PROTEIN_LIST_2_SUFFIX}"
    with open(protein_list_1, "r") as f:
        proteins_1 = set([line.strip().lower() for line in f.readlines()])
    with open(protein_list_2, "r") as f:
        proteins_2 = set([line.strip().lower() for line in f.readlines()])

    only_1 = proteins_1 - proteins_2
    only_2 = proteins_2 - proteins_1
    both = proteins_1 & proteins_2

    print(f"\n\n\n\n=================== {folder} ===================")
    print(f"{folder}{PROTEIN_LIST_1_SUFFIX} only: {len(only_1)}")
    for protein in only_2:
        print(protein)
    print(f"{folder}{PROTEIN_LIST_2_SUFFIX} only: {len(only_2)}")
    for protein in only_2:
        print(protein)
    print(f"Both: {len(both)}")
    for protein in both:
        print(protein)
