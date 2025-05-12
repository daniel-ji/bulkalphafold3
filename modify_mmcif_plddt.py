# modify the pLDDT values in a mmCIF file (bfactor column)

import os
import numpy as np
from Bio.PDB import MMCIFParser, MMCIFIO

def get_plddt_sliding_window_mmcif(input_file, residue_sliding_window=1, target_chain_id="PREY"):
    """
    Gets the corresponding mmCIF file with modified pLDDT (B-factor) values based on a sliding window. If it does not exist, it creates it:
        Modify the B-factor column in a mmCIF file based on the pLDDT values.
        Changes **only** the atom B-factors in the mmCIF file.
        The B-factor is set to the pLDDT value for each residue in the sliding window (i.e., 1 is the residue itself, 3 is the residue and its neighbors, etc.).
        Writes to a new mmCIF file with the same name as the input file but with the sliding window as a suffix.
    """

    output_file = input_file.replace(".cif", f"_plddt_window_{residue_sliding_window}.cif")
    
    if os.path.exists(output_file):
        # print(f"Output file {output_file} already exists. Skipping modification.")
        return output_file
    
    parser = MMCIFParser(QUIET=True)

    structure = parser.get_structure("structure", input_file)
    for model in structure:
        # print(f" Processing Model: {model.id}")
        for chain in model:
            if chain.id != target_chain_id:
                continue
            
            # print(f"  Found target chain: {chain.id}")
            residue_bfactors = [np.mean([atom.get_bfactor() for atom in residue]) for residue in chain]
            new_residue_bfactors = []
            for i in range(len(residue_bfactors)):
                start = max(0, i - residue_sliding_window // 2)
                end = min(len(residue_bfactors), i + residue_sliding_window // 2 + 1)
                new_bfactor_value = np.mean(residue_bfactors[start:end])
                new_residue_bfactors.append(new_bfactor_value)

            # write the new B-factors back to the structure
            for i, residue in enumerate(chain):
                for atom in residue:
                    atom.set_bfactor(new_residue_bfactors[i])

    print(f"Writing structure with new pLDDT scores to {output_file}", flush=True)
    io = MMCIFIO()
    io.set_structure(structure)
    io.save(output_file)

    return output_file