"""
A script to compare the filtered proteins from AF2 and AF3 predictions.
This script reads the protein names (or uniprot links) from two files (one for AF2 and one for AF3) and compares them to find unique and common proteins.
"""

from constants import FOLDERS

for folder in FOLDERS:
    af2_file = f"output_filtered_protein_list/{folder}_AF2_proteins.txt"
    af3_file = f"output_filtered_protein_list/{folder}_AF3_proteins.txt"
    with open(af2_file, "r") as f:
        af2_proteins = set([line.strip().lower() for line in f.readlines()])
    with open(af3_file, "r") as f:
        af3_proteins = set([line.strip().lower() for line in f.readlines()])
    
    only_af2 = af2_proteins - af3_proteins
    only_af3 = af3_proteins - af2_proteins
    both = af2_proteins & af3_proteins

    print(f"\n\n\n\n=================== {folder} ===================")
    print(f"AF2 only: {len(only_af2)}")
    for protein in only_af2:
        print(protein)
    print(f"AF3 only: {len(only_af3)}")
    for protein in only_af3:
        print(protein)
    print(f"Both: {len(both)}")
    for protein in both:
        print(protein)